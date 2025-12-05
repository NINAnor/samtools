# Helper function to read vector layers from PostGIS database

Helper function to easily read vector layers from PostGIS.

## Usage

``` r
db_read_vect(con, dsn, return_format = c("sf", "terra")[1])
```

## Arguments

- con:

  Connection to PostGIS.

- dsn:

  `[character(1)]`  
  Name of the layer in the PostGIS database. It might contain the name
  of the schema and the table, in the format `"schema_name.table_name"`.

- return_format:

  `[character="sf"]{"sf", "terra"}`  
  Format for the vector object retrieved. Either "sf" (default) or
  "terra" for a `SpatVector` object.

## Value

If `return_format = "sf"` (default), a `sf` vector object. If
`return_format = "terra"`, a `SpatVector` object.
