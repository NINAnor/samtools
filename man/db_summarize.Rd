% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_summarize.R
\name{db_summarize}
\alias{db_summarize}
\title{Summarize number, length, and area of features per area/region}
\usage{
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
}
\arguments{
\item{con}{Connection to PostGIS.}

\item{areas}{\verb{[character]} \cr Areas/polygons for which the summary will be made.
It should be in the format \code{"schema_name.table_name"}, representing a layer
located within the PostGIS database.}

\item{variables}{\verb{[vector,character]} \cr Array of layer names to be summarized, such as
roads, houses, etc. Layers should be presented in the format
\code{"schema_name.layer_name"}, representing layers located within the PostGIS database.}

\item{variable_names}{\verb{[vector,character]} \cr Names of the variables for the output
summary table, corresponding to the \code{variables} parameter. This should be an array
with the same length as \code{variables}.}

\item{area_geom_col}{\verb{[character="geom"]} \cr Character representing the name of
the geometry column to be considered for the \code{areas} layer.}

\item{variable_geom_col}{\verb{[vector,character="geom"]} \cr Character representing the name of
the geometry column to be considered for the \code{variables} layers. Alternatively,
it can be an array the same length as \code{variables}, with the geometry column name
for each of the variables.}

\item{areas_cols}{\verb{[vector,character=c("gid", "name_area")]} \cr Array of strings
with the name of the columns that should represent the area ID, to define each
of the areas in the \code{areas} layer.}

\item{areas_cols_groupby}{\verb{[vector,character]} \cr Array of strings
with the name of the columns that should used for grouping by the statistics.
Typically (and by default) this is the same as \code{area_cols}.}

\item{length_unit}{\verb{[character="m"]\{"m", "km"\}} \cr Unit for summarizing the length of the
Linestring vectors. By default, "m" (meters).}

\item{area_unit}{\verb{[character="km2"]\{"m2", "ha, "km2"\}} \cr Unit for summarizing
the area of the Polygon vectors. By default, "km2" (square kilometers).}

\item{return_spatial}{\verb{[logical(1)=TRUE]} \cr If \code{TRUE} (default), the results
is vector of \code{areas} with the summaries for amount, length, and area of layers
in the attribute table. If \code{FALSE}, the function returns only a \code{data.frame}
with area names and summary values, but no spatial component.}

\item{verbose}{\verb{[logical(1)=FALSE]}⁠\cr Should messages of the computation steps
be printed in the prompt along the computation? Default is \code{FALSE}.}
}
\value{
If \code{return_spatial = TRUE} (default), the the function produces a
vector of \code{areas} with the summaries for amount, length, and area of layers
in the attribute table. The geometry column is the same geometry from the \code{areas}
layer. If \code{return_spatial = FALSE}, the function returns only a \code{data.frame}
with area names and summary values, but no spatial component.
}
\description{
This function summarizes the number of vector features (for all vector types),
the total length of features (for Linestring vectors), and the total area of
features (for Polygon vectors) within a set of areas. Both the areas and the
vector layers/variables to be summarized should be in a PostGIS database.
}
\examples{
areas <- "sam_wrein_ancillary.reindeer_areas_official_2023"
variables <- c("sam_env.cabins_private_no", "sam_env.roads_private_renrein_no",
               "sam_env.reservoirs_all_no")
variable_names <- c("cabins_private", "roads_private", "reservoirs")
variable_geom_col <- c("geom", "geom", "geom_etrs33")
(summ <- db_summarize(con, areas, variables, variables_names,
                      variable_geom_col = variable_geom_col,
                      return_spatial = FALSE,
                      verbose = TRUE))

}
