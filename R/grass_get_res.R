#' Get resolution from a raster map in a GRASS GIS database
#'
#' @param map `[character]` \cr Name of the raster layer, possibly including the
#' mapset in the format `"map_name@mapset_name"` if in a mapset different from the
#' current mapset.
#'
#' @return A number representing the pixel size, based on the `nsres` parameter from
#' `r.info`, in units defined based on the project/location CRS. Should generally be meters.
#'
#' @examples
#' # connect to GRASS
#' NinaR::grassConnect()
#' # get resolution
#' grass_get_res(map = "master_grid_100m_norway@p_sam_tools")
#' grass_get_res(map = "master_grid_1km_norway@p_sam_tools")
#'
#' @export
grass_get_res <- function(map) {

  # r.info
  resol <- rgrass::execGRASS("r.info", map = map, flags = "g", intern = TRUE)

  # get resolution from the nsres
  resol <- grep("nsres", resol, value = TRUE) |>
    strsplit("=") |> sapply(function(x) x[2]) |>
    as.numeric()

  resol
}
