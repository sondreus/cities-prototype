# Ten kilometres from here

Open `index.html` directly in a browser. It includes all extracted data, works offline, and needs no server. Choose among twenty illustrative major cities; scroll, pinch, or use the buttons to zoom, and drag to pan. Hover or tap to inspect a cell. Keyboard: plus/minus, arrow keys, Home.

Red brightness reflects people per square kilometre inside the 10 km circle. Context uses grey. Both use the same nonlinear logarithmic density scale, capped at 50,000 people/km², for every city. Totals use native cell population counts, rather than the density values used for colouring.

## Data and calculation

Source: `population-density/source-data/global_pop_2026_UC_1km_R2024B_UA_v1.tif`, WorldPop 2026 R2024B v1, unconstrained, UN-adjusted. [Release statement](https://data.worldpop.org/repo/prj/Global_2015_2030/R2024B/doc/Global2_Release_Statement_R2024B_v1.pdf). The source is a modelled population estimate, with 30-arcsecond cells (approximately 1 km at the equator). It remains untouched.

`build.py` lists the exact reference coordinates. Cities are a geographically varied selection, not the twenty largest. The radius is from these points and includes land across administrative borders. It does not represent the city's official total.

Native cells are read in local windows, with no resampling of population counts. Distances use WGS84 geodesics. Boundary-cell counts are multiplied by estimated overlap fractions from a 32 × 32 subcell grid, weighted by geodesic subcell area. This assumes uniform population within each native cell. NoData contributes no population and is retained as NoData in tooltips. Density is native count divided by geodesic cell area in km². Maps display native-cell polygons in a local azimuthal equidistant projection.

`totals.csv` contains unrounded overlap-weighted totals and an alternative cell-centre sum for comparison. The UI rounds weighted totals to the nearest thousand. `data.json` stores counts, density, validity, overlap fractions, and projected cell corners; `index.html` embeds it. `template.html` is the editable UI source.

## Rebuild

Install Python packages `numpy`, `rasterio`, and `pyproj`, then run:

```sh
python build.py --raster /path/to/global_pop_2026_UC_1km_R2024B_UA_v1.tif
```

Pass `--samples 64` for a finer overlap estimate. Source year and description are fixed to this 2026 product; update those if adapting to another year or raster product.

## Verification

Checked all twenty dropdown choices, displayed totals, button and wheel zoom, reset, hover tooltips, drag, keyboard controls, and mobile layout in Chromium. No JavaScript errors. Compared 32 × 32 and 64 × 64 boundary sampling for London, Cairo, Mumbai, and Tokyo: total differences were below 0.004%. This tests numerical overlap stability, not accuracy of the underlying population model.
