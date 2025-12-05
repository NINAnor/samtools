# Extract values of raster using in parallel

Function to extract values of rasters along lines using parallel
computation. The function is actually not very optimized currently.

## Usage

``` r
extract_along(
  r,
  sl,
  ws,
  col_step_length = "dist",
  step_id = "use_ava_data_animals_id",
  prefix = "along_",
  summarize = TRUE
)

# S3 method for class 'SpatRaster'
extract_along(
  r,
  sl,
  ws,
  col_step_length = "dist",
  step_id = "use_ava_data_animals_id",
  prefix = "along_",
  summarize = TRUE
)

# S3 method for class 'data.frame'
extract_along(
  r,
  sl,
  ws,
  col_step_length = "dist",
  step_id = "use_ava_data_animals_id",
  prefix = "along_"
)
```

## Arguments

- r:

  `[SpatRast,data.frame]`  
  Either a set of rasters, in the `SpatRaster` format from terra; or a
  `data.frame` resulting of
  [`terra::extract`](https://rspatial.github.io/terra/reference/extract.html)
  for linestrings, for a set of environmental variables (defined here as
  the columns).

- sl:

  `[sf,SpatVector]`  
  Vector of step lines, with a column representing the step id (stratum
  or just step number) and a LINESTRING geometry for each step. Can be a
  `sf` or a `SpatVector` object.

- ws:

  List of which summary statistics are to be computed. Currently, only
  "mean", "max", "min", and "sum" are implemented.

- col_step_length:

  `[character="dist"]`  
  String with the name of the column representing step length.

- step_id:

  `[character="use_ava_data_animals_id"]`  
  String with the name of the column representing step ID.

- prefix:

  `[character="along_"]`  
  Prefix to be added to the extracted variable names.
