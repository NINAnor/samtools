% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_functions.R
\name{db_read_vect}
\alias{db_read_vect}
\title{Helper function to read vector layers from PostGIS database}
\usage{
db_read_vect(con, dsn, return_format = c("sf", "terra")[1])
}
\arguments{
\item{con}{Connection to PostGIS.}

\item{dsn}{\verb{[character(1)]} \cr Name of the layer in the PostGIS
database. It might contain the name of the schema and the table, in the format
\code{"schema_name.table_name"}.}

\item{return_format}{\verb{[character="sf"]\{"sf", "terra"\}} \cr Format for the
vector object retrieved. Either "sf" (default) or "terra" for a \code{SpatVector}
object.}
}
\value{
If \code{return_format = "sf"} (default), a \code{sf} vector object. If
\code{return_format = "terra"}, a \code{SpatVector} object.
}
\description{
Helper function to easily read vector layers from PostGIS.
}
