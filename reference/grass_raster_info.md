# Get information from a raster

This function gets information from a raster layer in GRASS GIS and
returns it as a list. It is a wrapper for module `r.info()` in GRASS
GIS.

## Usage

``` r
grass_raster_info(map)
```

## Arguments

- map:

  `[character]`  
  Name of the raster layer, possibly including the mapset in the format
  `"map_name@mapset_name"` if in a mapset different from the current
  mapset.

## Value

A list with 9 elements, extracted from `r.info`: the four coorinates
defining the bounding box of the raster (n, s, e, w); the raster
resolution, taken as nsres; the number of rows, cols, and cells, and the
type of data (e.g. CELL, FCELL).

## Examples

``` r
if(FALSE) {
  # connect to GRASS
  NinaR::grassConnect()
  # get raster parameters
  grass_raster_info(map = "master_grid_100m_norway@p_sam_tools")
}
```
