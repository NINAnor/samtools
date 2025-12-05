# Helper function to prepare queries to be sent to PostgreSQL

Helper function to easily concatenate strings and prepare queries for
PostgreSQL.

## Usage

``` r
db_make_query(...)
```

## Arguments

- ...:

  `[character]`  
  Strings representing the parts of a query.

## Value

A string with a query, to be used, for instance, within
[`DBI::dbExecute()`](https://dbi.r-dbi.org/reference/dbExecute.html).

## Examples

``` r
conditions <- paste0("\"KKOD\" IN ('5011', '5012', '5016', '5017')")

query1 <- db_make_query(
"CREATE TABLE IF NOT EXISTS sam_env.roads_public_se AS
SELECT *
FROM
\"Topography\".\"Sweden_Vagkartan_Roads_lines\"
WHERE ", paste(conditions, collapse = "AND "), ";")
query1
#> [1] "CREATE TABLE IF NOT EXISTS sam_env.roads_public_se AS SELECT * FROM \"Topography\".\"Sweden_Vagkartan_Roads_lines\" WHERE \"KKOD\" IN ('5011', '5012', '5016', '5017');"
# to execute, first one needs a connection 'con' to PostgreSQL, then it can be run as
# DBI::dbExecute(con, query1)
```
