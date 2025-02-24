#' Get information from a raster
#'
#' This function gets information from a raster layer in GRASS GIS and returns
#' it as a list. It is a wrapper for [rgrass::r.info()].
#'
#' @param map `[character]` \cr Name of the raster layer, possibly including the
#' mapset in the format `"map_name@mapset_name"` if in a mapset different from the
#' current mapset.
#'
#' @return A list with 9 elements, extracted from `r.info`:
#' the four coorinates defining the bounding box
#' of the raster (n, s, e, w); the raster resolution, taken as nsres; the number
#' of rows, cols, and cells, and the type of data (e.g. CELL, FCELL).
#'
#' @examples
#' if(FALSE) {
#'   # connect to GRASS
#'   NinaR::grassConnect()
#'   # get raster parameters
#'   grass_raster_info(map = "master_grid_100m_norway@p_sam_tools")
#' }
#'
#' @export
grass_raster_info <- function(map) {

  parms <- rgrass::execGRASS("r.info", map = map, flags = "g", intern = T)
  parms <- sapply(parms, function(x) strsplit(x, "=")[[1]][2])

  list(n = as.numeric(parms[1]), s = as.numeric(parms[2]),
       e = as.numeric(parms[3]), w = as.numeric(parms[4]),
       res = as.numeric(parms[5]), rows = as.numeric(parms[7]),
       cols = as.numeric(parms[8]), cells = as.numeric(parms[9]),
       datatype = as.character(parms[10]))
}
