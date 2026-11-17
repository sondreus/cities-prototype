"""Extract native WorldPop cells; generate a portable, offline map."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
import rasterio
from rasterio.windows import Window, from_bounds
from pyproj import Geod, Transformer

CITIES = [
    ('London', 'United Kingdom', 51.5074, -0.1278),
    ('Paris', 'France', 48.8566, 2.3522),
    ('New York', 'United States', 40.7128, -74.0060),
    ('Mexico City', 'Mexico', 19.4326, -99.1332),
    ('São Paulo', 'Brazil', -23.5505, -46.6333),
    ('Buenos Aires', 'Argentina', -34.6037, -58.3816),
    ('Lagos', 'Nigeria', 6.5244, 3.3792),
    ('Cairo', 'Egypt', 30.0444, 31.2357),
    ('Kinshasa', 'DR Congo', -4.3250, 15.3222),
    ('Johannesburg', 'South Africa', -26.2041, 28.0473),
    ('Istanbul', 'Türkiye', 41.0082, 28.9784),
    ('Karachi', 'Pakistan', 24.8607, 67.0011),
    ('Mumbai', 'India', 19.0760, 72.8777),
    ('Delhi', 'India', 28.6139, 77.2090),
    ('Dhaka', 'Bangladesh', 23.8103, 90.4125),
    ('Beijing', 'China', 39.9042, 116.4074),
    ('Shanghai', 'China', 31.2304, 121.4737),
    ('Tokyo', 'Japan', 35.6895, 139.6917),
    ('Jakarta', 'Indonesia', -6.2088, 106.8456),
    ('Manila', 'Philippines', 14.5995, 120.9842),
]
SOURCE_URL = 'https://data.worldpop.org/repo/prj/Global_2015_2030/R2024B/doc/Global2_Release_Statement_R2024B_v1.pdf'

def extract(src, city, samples):
    name, country, lat, lon = city
    geod = Geod(ellps='WGS84')
    project = Transformer.from_crs('EPSG:4326', f'+proj=aeqd +lat_0={lat} +lon_0={lon} +datum=WGS84 +units=m', always_xy=True)
    az = np.arange(360)
    blon, blat, _ = geod.fwd(np.full(360, lon), np.full(360, lat), az, np.full(360, 38000))
    win = from_bounds(blon.min(), blat.min(), blon.max(), blat.max(), src.transform)
    col, row = int(np.floor(win.col_off)), int(np.floor(win.row_off))
    win = Window(col, row, int(np.ceil(win.col_off+win.width))-col, int(np.ceil(win.row_off+win.height))-row)
    values = src.read(1, window=win, masked=True)
    valid = ~np.ma.getmaskarray(values)
    pop = values.filled(0).astype(float)
    assert np.all(pop >= 0) and np.all(np.isfinite(pop))
    h, w = pop.shape
    tr = src.window_transform(win)
    lons = tr.c + np.arange(w+1)*tr.a
    lats = tr.f + np.arange(h+1)*tr.e
    lng, lt = np.meshgrid(lons, lats)
    x, y = project.transform(lng, lt)
    row_areas = []
    for r in range(h):
        area, _ = geod.polygon_area_perimeter([lons[0],lons[1],lons[1],lons[0]], [lats[r],lats[r],lats[r+1],lats[r+1]])
        row_areas.append(abs(area)/1e6)
    density = pop / np.array(row_areas)[:,None]
    clng, clt = np.meshgrid((lons[:-1]+lons[1:])/2, (lats[:-1]+lats[1:])/2)
    _, _, distance = geod.inv(np.full(pop.shape, lon), np.full(pop.shape, lat), clng, clt)
    # A 2 km margin exceeds the diagonal of any native 30-arcsecond cell.
    fraction = (distance < 8000).astype(float)
    boundary = np.argwhere((distance >= 8000) & (distance <= 12000))
    offsets = (np.arange(samples)+.5)/samples
    for r, c in boundary:
        slon, slat = np.meshgrid(lons[c]+offsets*tr.a, lats[r]+offsets*tr.e)
        _, _, d = geod.inv(np.full(slon.shape, lon), np.full(slat.shape, lat), slon, slat)
        # Geodesic subcell areas vary with latitude; use those as weights.
        areas = []
        for t in range(samples):
            south = lats[r]+(t+1)*tr.e/samples
            north = lats[r]+t*tr.e/samples
            area, _ = geod.polygon_area_perimeter([lons[c], lons[c]+tr.a/samples, lons[c]+tr.a/samples, lons[c]], [north,north,south,south])
            areas.append(abs(area))
        weights = np.repeat(np.array(areas)[:,None], samples, axis=1)
        fraction[r,c] = np.sum(weights*(d <= 10000))/weights.sum()
    total = float(np.sum(pop*fraction))
    assert total > 0 and np.all((fraction >= 0) & (fraction <= 1))
    assert not np.any(fraction[[0,-1],:]) and not np.any(fraction[:,[0,-1]])
    print(f'{name}: {total:,.0f} ({h}×{w} cells)', flush=True)
    return dict(name=name, country=country, lat=lat, lon=lon, total=total,
                centroidTotal=float(np.sum(pop[distance<=10000])), w=w,h=h,
                x=np.round(x,1).flatten().tolist(),y=np.round(y,1).flatten().tolist(),
                pop=np.round(pop,3).flatten().tolist(),valid=valid.flatten().astype(int).tolist(),
                density=np.round(density,3).flatten().tolist(),
                fraction=np.round(fraction,6).flatten().tolist())

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--raster', type=Path, default=Path('/Users/sondresolstad/Github/population-density/source-data/global_pop_2026_UC_1km_R2024B_UA_v1.tif'))
    parser.add_argument('--samples',type=int,default=32)
    args=parser.parse_args()
    if args.samples<1: parser.error('--samples must be positive')
    out=Path(__file__).parent
    with rasterio.open(args.raster) as src:
        assert src.crs.to_epsg()==4326 and src.count==1
        data=dict(year=2026,source=args.raster.name,sourceUrl=SOURCE_URL,samples=args.samples,
                  cities=[extract(src,c,args.samples) for c in CITIES])
    (out/'data.json').write_text(json.dumps(data,separators=(',',':')))
    template=(out/'template.html').read_text()
    (out/'index.html').write_text(template.replace('__DATA__',json.dumps(data,separators=(',',':'))))
    with (out/'totals.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=['name','country','lat','lon','total','centroidTotal'])
        writer.writeheader()
        writer.writerows({k:c[k] for k in writer.fieldnames} for c in data['cities'])

if __name__=='__main__': main()
