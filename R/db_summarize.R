#' Summarize number, length, and area of features per area/region
#'
#' This function summarizes the number of vector features (for all vector types),
#' the total length of features (for Linestring vectors), and the total area of
#' features (for Polygon vectors) within a set of areas. Both the areas and the
#' vector layers/variables to be summarized should be in a PostGIS database.
#'
#' @param con Connection to PostGIS.
#' @param areas `[character]` \cr Areas/polygons for which the summary will be made.
#' It should be in the format `"schema_name.table_name"`, representing a layer
#' located within the PostGIS database.
#' @param variables `[vector,character]` \cr Array of layer names to be summarized, such as
#' roads, houses, etc. Layers should be presented in the format
#' `"schema_name.layer_name"`, representing layers located within the PostGIS database.
#' @param variable_names `[vector,character]` \cr Names of the variables for the output
#' summary table, corresponding to the `variables` parameter. This should be an array
#' with the same length as `variables`.
#' @param area_geom_col `[character="geom"]` \cr Character representing the name of
#' the geometry column to be considered for the `areas` layer.
#' @param variable_geom_col `[vector,character="geom"]` \cr Character representing the name of
#' the geometry column to be considered for the `variables` layers. Alternatively,
#' it can be an array the same length as `variables`, with the geometry column name
#' for each of the variables.
#' @param areas_cols `[vector,character=c("gid", "name_area")]` \cr Array of strings
#' with the name of the columns that should represent the area ID, to define each
#' of the areas in the `areas` layer.
#' @param areas_cols_groupby `[vector,character]` \cr Array of strings
#' with the name of the columns that should used for grouping by the statistics.
#' Typically (and by default) this is the same as `area_cols`.
#' @param length_unit `[character="m"]{"m", "km"}` \cr Unit for summarizing the length of the
#' Linestring vectors. By default, "m" (meters).
#' @param area_unit `[character="km2"]{"m2", "ha, "km2"}` \cr Unit for summarizing
#' the area of the Polygon vectors. By default, "km2" (square kilometers).
#' @param return_spatial `[logical(1)=TRUE]` \cr If `TRUE` (default), the results
#' is vector of `areas` with the summaries for amount, length, and area of layers
#' in the attribute table. If `FALSE`, the function returns only a `data.frame`
#' with area names and summary values, but no spatial component.
#' @param verbose `[logical(1)=FALSE]`⁠\cr Should messages of the computation steps
#' be printed in the prompt along the computation? Default is `FALSE`.
#'
#' @return If `return_spatial = TRUE` (default), the the function produces a
#' vector of `areas` with the summaries for amount, length, and area of layers
#' in the attribute table. The geometry column is the same geometry from the `areas`
#' layer. If `return_spatial = FALSE`, the function returns only a `data.frame`
#' with area names and summary values, but no spatial component.
#'
#' @examples
#' areas <- "sam_wrein_ancillary.reindeer_areas_official_2023"
#' variables <- c("sam_env.cabins_private_no", "sam_env.roads_private_renrein_no",
#'                "sam_env.reservoirs_all_no")
#' variable_names <- c("cabins_private", "roads_private", "reservoirs")
#' variable_geom_col <- c("geom", "geom", "geom_etrs33")
#' (summ <- db_summarize(con, areas, variables, variables_names,
#'                       variable_geom_col = variable_geom_col,
#'                       return_spatial = FALSE,
#'                       verbose = TRUE))
#'
#' @export
db_summarize <- function(con, areas, variables, variable_names,
                         area_geom_col = "geom", variable_geom_col = "geom",
                         areas_cols = c("gid", "name_area"),
                         areas_cols_groupby = areas_cols,
                         length_unit = c("m", "km")[2],
                         area_unit = c("m2", "ha", "km2")[3],
                         return_spatial = TRUE,
                         verbose = FALSE) {

  # set input, comment later
  # areas <- "sam_wrein_ancillary.reindeer_areas_official_2023"
  # variables <- c("sam_env.cabins_private_no", "sam_env.roads_private_renrein_no",
  #                "sam_env.reservoirs_all_no")
  # variable_names <- c("cabins_private", "roads_private", "reservoirs")
  ### check that variables and variable_names have similar lengths

  # check number of geom values vs number of variables
  # variable_geom_col <- "geom"
  if(length(variable_geom_col) == 1 & length(variables) > 1)
    variable_geom_col  <- rep(variable_geom_col, length(variables))
  # variable_geom_col <- c("geom", "geom", "geom_etrs33")

  # check geometry names first
  for(i in seq_along(variables)) {
    # get geometry column name
    schema_tab <- strsplit(variables[i], ".", fixed = TRUE)[[1]]
    q0 <- db_make_query("SELECT * from geometry_columns
                         WHERE f_table_schema = '", schema_tab[1],
                         "' AND f_table_name = '", schema_tab[2], "';")
    geoms <- DBI::dbGetQuery(con, q0)

    # check if table is empty
    if(nrow(geoms) < 1) stop(paste0("Check the name of the schema and table for the variable '", variables[i], "'.")) else
      # check if geometry exists in the table
      if(!(variable_geom_col[i] %in% geoms$f_geometry_column))
        stop(paste0("The geometry '", variable_geom_col[i], "' does not exist in the table '", variables[i], "' The geometry columns in this table are: '", paste(geoms$f_geometry_column, collapse = ","), "'."))
  }

  #---
  # set units
  # divisor for lengths
  div_length <- 1
  if(length_unit == "km") div_length <- 1000
  # divisor for areas
  div_area <- 1
  if(area_unit == "ha") div_area <- 1e4
  if(area_unit == "km2") div_area <- 1e6

  # initialize list of tables
  tbs <- list()

  #---
  # get total area
  q3 <- db_make_query("SELECT ", paste0("a.", areas_cols, collapse = ", "), ", ",
                      "'", areas, "' AS source_area, ",
                      "ST_Area(a.", area_geom_col, ")/", div_area, " AS area_total_", area_unit, " ",
                      "FROM ", areas, " a ",
                      "GROUP by ", paste0("a.", areas_cols_groupby, collapse = ", "), ";")
  print(q3)
  tb_tot <- DBI::dbGetQuery(con, q3)
  # tbs[[1]] <- tb_tot

  #---
  # loop over variables
  i <- 1
  for(i in seq_along(variables)) {

    # select layer
    v <- variables[i]
    v_name <- variable_names[i]
    v_geom <- variable_geom_col[i]

    if(verbose) print(paste0("Computing summaries for ", v_name, "..."))

    # get geometry type
    q1 <- db_make_query("SELECT ST_GeometryType(v.", v_geom, ") ",
                        "FROM ", v, " v ",
                        "LIMIT 1")
    geom_type <- DBI::dbGetQuery(con, q1)[1,1]

    # for points
    if(grepl("POINT", geom_type, ignore.case = TRUE)) {

      # compute number of points
      q2 <- db_make_query("SELECT ", paste0("a.", areas_cols, collapse = ", "), ", ",
                          "'", areas, "' AS source_area, ",
                          "'", v, "' AS source_", v_name, ", ",
                          "COUNT(v.*) AS count_", v_name, #", ",
                          # "count_", v_name, "/(ST_Area(a.", area_geom_col, ")/", div_area, ") ",
                          #"AS density_", v_name, "_per_", area_unit, " ",
                          ifelse(return_spatial, paste0(", a.", area_geom_col, " "), " "),
                          "FROM ", areas, " a ",
                          "LEFT JOIN ", v, " v ON ", "ST_Within(v.", v_geom, ", a.", area_geom_col, ") ",
                          "GROUP by ", paste0("a.", areas_cols_groupby, collapse = ", "), ";")

      # run query
      tb <- DBI::dbGetQuery(con, q2)

      # compute density
      tb <- tb |>
        dplyr::left_join(tb_tot) |>
        mutate(!!rlang::sym(paste0("density_", v_name, "_per_", area_unit)) := !!rlang::sym(paste0("count_", v_name))/!!rlang::sym(paste0("area_total_", area_unit)))
    }

    # for lines
    if(grepl("LINESTRING", geom_type, ignore.case = TRUE)) {

      # compute number and length of lines after cropped
      q2 <- db_make_query("SELECT ", paste0("a.", areas_cols, collapse = ", "), ", ",
                          "'", areas, "' AS source_area, ",
                          "'", v, "' AS source_", v_name, ", ",
                          "COUNT(ST_Intersection(v.", v_geom, ", a.", area_geom_col, ")) AS count_", v_name, ", ",
                          "SUM(ST_Length(ST_Intersection(v.", v_geom, ", a.", area_geom_col, ")))/", div_length,
                          " AS length_", v_name, "_", length_unit,
                          ifelse(return_spatial, paste0(", a.", area_geom_col, " "), " "),
                          "FROM ", areas, " a ",
                          "LEFT JOIN ", v, " v ON ",
                          "ST_Intersects(v.", v_geom, ", a.", area_geom_col, ") ",
                          "GROUP by ", paste0("a.", areas_cols_groupby, collapse = ", "), ";")

      # run query
      tb <- DBI::dbGetQuery(con, q2)

      # compute density
      tb <- tb |>
        dplyr::left_join(tb_tot) |>
        mutate(!!rlang::sym(paste0("length_", v_name, "_", length_unit)) := ifelse(is.na(!!rlang::sym(paste0("length_", v_name, "_", length_unit))), 0, !!rlang::sym(paste0("length_", v_name, "_", length_unit))),
               !!rlang::sym(paste0("density_", v_name, "_", length_unit, "_per_", area_unit)) := !!rlang::sym(paste0("length_", v_name, "_", length_unit))/!!rlang::sym(paste0("area_total_", area_unit)))

    }

    # for polygons
    if(grepl("POLYGON", geom_type, ignore.case = TRUE)) {

      # compute number and area of polygons after cropped
      q2 <- db_make_query("SELECT ", paste0("a.", areas_cols, collapse = ", "), ", ",
                          "'", areas, "' AS source_area, ",
                          "'", v, "' AS source_", v_name, ", ",
                          "COUNT(ST_Intersection(v.", v_geom, ", a.", area_geom_col, ")) AS count_", v_name, ", ",
                          "SUM(ST_Area(ST_Intersection(v.", v_geom, ", a.", area_geom_col, ")))/", div_area,
                          " AS area_", v_name, "_", area_unit,
                          ifelse(return_spatial, paste0(", a.", area_geom_col, " "), " "),
                          "FROM ", areas, " a ",
                          "LEFT JOIN ", v, " v ON ",
                          "ST_Intersects(v.", v_geom, ", a.", area_geom_col, ") ",
                          "GROUP by ", paste0("a.", areas_cols_groupby, collapse = ", "), ";")

      # run query
      tb <- DBI::dbGetQuery(con, q2)

      # compute density
      tb <- tb |>
        dplyr::left_join(tb_tot) |>
        mutate(!!rlang::sym(paste0("area_", v_name, "_", area_unit)) := ifelse(is.na(!!rlang::sym(paste0("area_", v_name, "_", area_unit))), 0, !!rlang::sym(paste0("area_", v_name, "_", area_unit))),
               !!rlang::sym(paste0("density_", v_name, "_", area_unit, "_per_", area_unit)) := !!rlang::sym(paste0("area_", v_name, "_", area_unit))/!!rlang::sym(paste0("area_total_", area_unit)))

    }

    tbs[[i]] <- tb
  }

  # columns to aggregate
  if(!return_spatial) areas_cols_groupby <- areas_cols_groupby[!areas_cols_groupby %in% area_geom_col]
  joinby <- c(areas_cols_groupby, "source_area", "area_total_km2")
  # bind cols
  tb_out <- tbs[[1]]
  if(length(tbs) > 1) {
    for(i in 2:length(tbs)) {
      tb_out <- tb_out |>
        dplyr::left_join(tbs[[i]], by = setNames(nm = joinby))
    }
  }

  # return table
  tb_out
}

# ----------------------
#   -- Get type of feature
# -- pts
# SELECT ST_GeometryType(v.geom)
# FROM sam_env.cabins_private_no v
# -- lin
# SELECT ST_GeometryType(v.geom)
# FROM sam_env.roads_private_renrein_no v
# --- pol
# SELECT ST_GeometryType(v.geom_etrs33)
# FROM sam_env.reservoirs_all_no v
#
# -----------------------
#   -- Number of points
# SELECT a.gid, a.name_area,
# 'sam_wrein_ancillary.reindeer_areas_official_2023' AS source_area,
# 'sam_env.cabins_no' AS source_cabins,
# Count(v.*) AS count_cabins,
# a.geom
# --Sum(ST_Area(ST_Within(v.geom, a.geom)))
# FROM sam_wrein_ancillary.reindeer_areas_official_2023 a
# LEFT JOIN sam_env.cabins_private_no v ON ST_Within(v.geom, a.geom)
# GROUP by a.gid, a.name_area
#
# -- number and length of lines
# SELECT a.gid, a.name_area,
# 'sam_wrein_ancillary.reindeer_areas_official_2023' AS source_area,
# 'sam_env.roads_minor_no' AS source_roads_minor,
# --Count(a.geom)--,
# COUNT(ST_Intersection(v.geom, a.geom)) AS count_roads_minor,
# SUM(ST_Length(ST_Intersection(v.geom, a.geom)))/1000 AS length_roads_minor
# --a.geom
# --Sum(ST_Area(ST_Within(v.geom, a.geom)))
# FROM sam_wrein_ancillary.reindeer_areas_official_2023 a
# LEFT JOIN sam_env.roads_minor_no v ON ST_Intersects(v.geom, a.geom)
# GROUP by a.gid, a.name_area
#
# -- number and area of polygons
# SELECT a.gid, a.name_area,
# 'sam_wrein_ancillary.reindeer_areas_official_2023' AS source_area,
# 'sam_env.reservoirs_all_no' AS source_reservoirs,
# COUNT(ST_Intersection(v.geom_etrs33, a.geom)) AS count_reservoirs,
# SUM(ST_Area(ST_Intersection(v.geom_etrs33, a.geom)))/1e6 AS length_reservoirs
# --a.geom
# FROM sam_wrein_ancillary.reindeer_areas_official_2023 a
# LEFT JOIN sam_env.reservoirs_all_no v ON ST_Intersects(v.geom_etrs33, a.geom)
# GROUP by a.gid, a.name_area
#
# -- identify type of geometry
# -- set type of summary stats - count, length, area
# -- set unit of output (m, km, m2, ha, km2)
# -- compute total area and proportion/relative value
# -- cbind all infrastructure types
# -- options - return or not a geometry, group by again or not
