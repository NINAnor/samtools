# Summarize number, length, and area of features per area/region

This function summarizes the number of vector features (for all vector
types), the total length of features (for Linestring vectors), and the
total area of features (for Polygon vectors) within a set of areas. Both
the areas and the vector layers/variables to be summarized should be in
a PostGIS database.

## Usage

``` r
db_summarize(
  con,
  areas,
  variables,
  variable_names,
  area_geom_col = "geom",
  variable_geom_col = "geom",
  areas_cols = c("gid", "name_area"),
  areas_cols_groupby = areas_cols,
  length_unit = c("m", "km")[2],
  area_unit = c("m2", "ha", "km2")[3],
  return_spatial = TRUE,
  verbose = FALSE
)
```

## Arguments

- con:

  Connection to PostGIS.

- areas:

  `[character]`  
  Areas/polygons for which the summary will be made. It should be in the
  format `"schema_name.table_name"`, representing a layer located within
  the PostGIS database.

- variables:

  `[vector,character]`  
  Array of layer names to be summarized, such as roads, houses, etc.
  Layers should be presented in the format `"schema_name.layer_name"`,
  representing layers located within the PostGIS database.

- variable_names:

  `[vector,character]`  
  Names of the variables for the output summary table, corresponding to
  the `variables` parameter. This should be an array with the same
  length as `variables`.

- area_geom_col:

  `[character="geom"]`  
  Character representing the name of the geometry column to be
  considered for the `areas` layer.

- variable_geom_col:

  `[vector,character="geom"]`  
  Character representing the name of the geometry column to be
  considered for the `variables` layers. Alternatively, it can be an
  array the same length as `variables`, with the geometry column name
  for each of the variables.

- areas_cols:

  `[vector,character=c("gid", "name_area")]`  
  Array of strings with the name of the columns that should represent
  the area ID, to define each of the areas in the `areas` layer.

- areas_cols_groupby:

  `[vector,character]`  
  Array of strings with the name of the columns that should used for
  grouping by the statistics. Typically (and by default) this is the
  same as `area_cols`.

- length_unit:

  `[character="m"]{"m", "km"}`  
  Unit for summarizing the length of the Linestring vectors. By default,
  "m" (meters).

- area_unit:

  `[character="km2"]{"m2", "ha, "km2"}`  
  Unit for summarizing the area of the Polygon vectors. By default,
  "km2" (square kilometers).

- return_spatial:

  `[logical(1)=TRUE]`  
  If `TRUE` (default), the results is vector of `areas` with the
  summaries for amount, length, and area of layers in the attribute
  table. If `FALSE`, the function returns only a `data.frame` with area
  names and summary values, but no spatial component.

- verbose:

  `[logical(1)=FALSE]`⁠  
  Should messages of the computation steps be printed in the prompt
  along the computation? Default is `FALSE`.

## Value

If `return_spatial = TRUE` (default), the the function produces a vector
of `areas` with the summaries for amount, length, and area of layers in
the attribute table. The geometry column is the same geometry from the
`areas` layer. If `return_spatial = FALSE`, the function returns only a
`data.frame` with area names and summary values, but no spatial
component.

## Examples

``` r
areas <- "sam_wrein_ancillary.reindeer_areas_official_2023"
variables <- c("sam_env.cabins_private_no", "sam_env.roads_private_renrein_no",
               "sam_env.reservoirs_all_no")
variable_names <- c("cabins_private", "roads_private", "reservoirs")
variable_geom_col <- c("geom", "geom", "geom_etrs33")
(summ <- db_summarize(con, areas, variables, variables_names,
                      variable_geom_col = variable_geom_col,
                      return_spatial = FALSE,
                      verbose = TRUE))
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'conn' in selecting a method for function 'dbGetQuery': object 'con' not found
```
