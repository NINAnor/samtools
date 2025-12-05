# Helper function to upload vector layers to PostGIS database

Helper function to easily upload vector layers to PostGIS, while adding
a primary/unique key and possibly a comment to the layer. This function
also ensures that the owner of the layer has automatically write to
modify the layer as wanted.

## Usage

``` r
db_write_vect(vect, con, dsn, pkey_column = "pkey", comment = "")
```

## Arguments

- vect:

  Vector object in R, in `sf` format.

- con:

  Connection to PostGIS.

- dsn:

  `[character(1)]`  
  Name of the layer after saved in the PostGIS database. It might
  contain the name of the schema and the table, in the format
  `"schema_name.table_name"`.

- pkey_column:

  `[character(1)="pkey"]`  
  Name of the primary key column to be created. Default is `"pkey"`.

- comment:

  `[character=""]`  
  Comment to be added to the layer. Typically a description of the
  layer, how it was processed, the source, and other relevant
  information.

## Value

The function does not return anything, but writes de vector into the
PostGIS database.
