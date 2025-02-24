#' Set the region in a GRASS GIS session
#'
#' Helper function to set the computation region of an open GRASS GIS
#' session connected to the R session, e.g. through [rgrass::initGRASS()].
#' To know more about the arguments, check: https://grass.osgeo.org/grass82/manuals/g.region.html
#'
#' @examples
#' if(FALSE) {
#'   grass_set_region(res = 30)
#' }
#'
#' @export
grass_set_region <- function(flags = c("p"), ...) {
  rgrass::execGRASS("g.region", parameters = list(...), flags = flags)
}
