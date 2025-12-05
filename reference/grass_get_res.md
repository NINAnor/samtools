# Get resolution from a raster map in a GRASS GIS database

This function gets the resolution of a raster layer in GRASS GIS.

## Usage

``` r
grass_get_res(map)
```

## Arguments

- map:

  `[character]`  
  Name of the raster layer, possibly including the mapset in the format
  `"map_name@mapset_name"` if in a mapset different from the current
  mapset.

## Value

A number representing the pixel size, based on the `nsres` parameter
from `r.info`, in units defined based on the project/location CRS.
Should generally be meters.

## Examples

``` r
if(FALSE) {
  # connect to GRASS
  NinaR::grassConnect()
  # get resolution
  grass_get_res(map = "master_grid_100m_norway@p_sam_tools")
  grass_get_res(map = "master_grid_1km_norway@p_sam_tools")
}
```
