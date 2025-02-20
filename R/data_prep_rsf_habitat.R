#' Prepare data for RSF analyses for wild and semi-domesticated reindeer within SAM
#'
#' This function takes in the already annotated data for RSF with either wild or semi-domesticated
#' reindeer and further prepare columns and variables for the model fitting procedure.
#' The function is suited for point resource selection functions (RSF) only.
#'
#' @param dat `[data.frame]` \cr Data set annotated with environmental covariates, almost
#' ready for analysis.
#' @param season `[string]{"sum", "cal", "win"}` \cr Season of interest.
#' One of "sum", "cal", or "win".
#' @param prediction `[logical(1)=FALSE]` \cr Additional changes that should
#' only be done in the grid for prediction.
#' @param land_cover_factor `[logical(1)=FALSE]` \cr Logical variable stating whether
#' the land cover variable should be treated as a factor (one single column with
#' land cover type names) or not (as multiple columns as dummy variables).
#' @param land_cover `[string(1)="norut_smd"]{"norut_smd", "norut", "smd", "nmd"}` \cr
#' Which land cover map should be used for analysis. It prepares the corresponding
#' classes, using the `ref_landcover` as the reference class.
#' @param ref_landcover `[string="heathland"]` \cr Reference class for land cover
#' variables.
#' @param include_zoi_nearest `[logical(1)=FALSE]` \cr If `TRUE`, variables representing
#' the zone of influence (ZOI) of the nearest feature are also computed, based on
#' the distance to the nearest features. Only relevant if `prediction = TRUE`.
#' @param species `[string="wrein"]{"wrein", "trein"}` \cr The species/management system
#' for which we should prepare the data for. Either `"wrein"` for wild reindeer or
#' `"trein"` for semi-domesticated reindeer. It only implies some differences in data
#' which are particular to each of the species/management systems.
#' @param prefix `[string=""]` \cr Prefix to be added to the variable/column names
#' of the dataset. Default is `""`, but it could be e.g. `"startpt_"`, `"endpt_"`,
#' or `"along_"` for variables extracted at the starting or ending point or along
#' steps.
#' @param use_binary_vars_as_categorical `[logical(1)=FALSE]` \cr If `TRUE`, binary
#' variables (e.g. lakes, reservoirs) should be treated as categorical. Otherwise,
#' they are treated as numerical.
#'
#' @export
data_prep_rsf_habitat_rein <- function(dat, season,
                                       prediction = FALSE,
                                       fixwind = TRUE,
                                       land_cover_factor = FALSE,
                                       land_cover = c("norut_smd", "norut", "smd", "nmd")[1],
                                       ref_landcover = "heathland",
                                       include_zoi_nearest = FALSE,
                                       species = c("wrein", "trein")[1],
                                       prefix = "",
                                       use_binary_vars_as_categorical = FALSE) {

  #---
  # rename variables and compute new, derived variables

  #--
  # all infrastructure
  # compute ZOI of the nearest feature for prediction
  if(prediction & include_zoi_nearest) {

    # add ZOI nearest
    cols_euclidean <- colnames(dat) |>
      grep(pattern = "euclidean", value = TRUE)
    radii <- c(100, 250, 500, 1000, 2500, 5000, 10000)
    # test
    # dd <- oneimpact::dist_decay(ldat[[dat]][[seas]][["houses_euclidean"]], radius = radii[7],
    #                       type = "bartlett")
    # plot(ldat[[dat]][[seas]][["houses_euclidean"]], dd)
    for(r in radii) {
      for(cc in cols_euclidean) {
        new_col <- base::sub("euclidean", paste0("nearest", r), cc)
        dat[[new_col]] <- oneimpact::dist_decay(dat[[cc]],
                                                radius = r,
                                                type = "bartlett")
      }
    }

  }

  #--
  # transport and urban

  # houses -- nothing here
  string <- "houses"
  cols_rhou_n <- grep(string, names(dat))
  names(dat)[cols_rhou_n]
  # roads major
  string <- "roads_major"
  cols_rmaj_n <- grep(string, names(dat))
  names(dat)[cols_rmaj_n]
  roads_seas <- ifelse(season == "sum", "sum", "win")
  names(dat)[cols_rmaj_n] <- sub(paste0(string, "_", roads_seas), string, names(dat)[cols_rmaj_n])
  # roads minor
  string <- "roads_minor"
  cols_rmin_n <- grep(string, names(dat))
  names(dat)[cols_rmin_n]
  roads_seas <- ifelse(season == "sum", "sum", "win")
  names(dat)[cols_rmin_n] <- sub(paste0(string, "_", roads_seas), string, names(dat)[cols_rmin_n])
  # roads public
  string <- "roads_public"
  cols_rpub_n <- grep(string, names(dat))
  names(dat)[cols_rpub_n]
  roads_seas <- ifelse(season == "sum", "sum", "win")
  names(dat)[cols_rpub_n] <- sub(paste0(string, "_", roads_seas), string, names(dat)[cols_rpub_n])
  # roads private
  string <- "roads_private"
  cols_rpriv_n <- grep(string, names(dat))
  names(dat)[cols_rpriv_n]
  roads_seas <- ifelse(season == "sum", "sum", "win")
  names(dat)[cols_rpriv_n] <- sub(paste0(string, "_", roads_seas), string, names(dat)[cols_rpriv_n])
  # railways
  string <- "railways"
  cols_rrail_n <- grep(string, names(dat))
  names(dat)[cols_rrail_n]
  rail_seas <- ifelse(season == "sum", "sum", "win") ## check here with Anna and Per
  names(dat)[cols_rrail_n] <- sub(paste0(string, "_", rail_seas), string, names(dat)[cols_rrail_n])

  #--
  # industry

  # industrial buildings
  string <- "industrial_buildings"
  cols_rind_n <- grep(string, names(dat))
  names(dat)[cols_rind_n]
  names(dat)[cols_rind_n] <- sub(paste0(string), "ind_build", names(dat)[cols_rind_n])
  # mining - nothing here
  string <- "mining"
  cols_min_n <- grep(string, names(dat))
  names(dat)[cols_min_n]
  # power lines
  string <- "power_lines"
  cols_rpl_n <- grep(string, names(dat))
  names(dat)[cols_rpl_n]
  names(dat)[cols_rpl_n] <- sub(paste0(string), "powerlines", names(dat)[cols_rpl_n])
  # Create radii (zone of) of influence for distance powerline
  # radii <- c(100, 250, 500, 1000, 2500, 5000, 10000, 20000)
  # for (i in radii){
  #   tmp <- ifelse(dat$dist_powerlines>i, i, dat$dist_powerlines)
  #   dat <- cbind(dat, tmp)
  #   names(dat)[ ncol(dat)] <- paste0("dist_powerlines_", i)#
  # }
  # reservoirs - transform into factor
  if(use_binary_vars_as_categorical) {
    string <- "reservoirs_hydro"
    cols_hyd_n <- grep(string, names(dat))
    names(dat)[cols_hyd_n]
    dat[names(dat)[cols_hyd_n]] <- factor(dat[[names(dat)[cols_hyd_n]]])
  }

  # if (fixwind){
  #   ## Create several windpark var  NB: DO NOT CONSIDER "CONSTRUCTION" - USE ACTIVE OR NOT (IN SUMMER THERE ARE NO CONSTRUCTION WORKS, AND WE DO NOT HAVE MONHTS SPECIFIED)
  #   #add year and windpark_activity for the available points
  #   if (prediction){
  #     dat$windpark_dist <- ifelse(is.na(dat$windpark_dist), 20000, dat$windpark_dist)
  #     dat$windpark_activity <- ifelse(dat$windpark_dist>=20000, "none", "active")
  #   }
  #
  #   # distance only from windfarms that are active or under construction
  #   dat$windpark_dist_act <- ifelse(dat$windpark_activity=='none',20000, dat$windpark_dist) # if there is not yet a windpark active or under construciton, set dist to 20 km
  #   dat$log_windpark_dist_act <- log10(dat$windpark_dist_act+1) # same, log transformed
  #   dat$log_windpark_dist <- log10(dat$windpark_dist+1) # log transformed distance to all windparks (also none active ones)
  #
  #   # in-outside windpark polygon
  #   dat$inside_windpark <- as.numeric(dat$windpark_dist==0) # 0=inside windpark (no matter if windpark exists or not!), 1=outside
  #   dat$inside_windpark_act <- as.numeric(dat$windpark_dist_act==0) # 0=inside windpark active/under constr , 1=outside
  #
  #   # Create radii (zone of) of influence for distance variables - 3 variables:
  #   ## windpark_100 = distance to nearest windpower only active or under construction. Before activity/construciotn equal > 20 km from wind. So we can use this var alone , no need for interation with activity
  #   ## windpark_all_100 = distance to all windpower - before, active or under construction
  #   ## windpark_activity_100 = what is the status of the nearest windpark: under construction, active or non active/too far
  #   radii <- c(100, 250, 500, 1000, 2500, 5000, 10000, 20000)
  #   for (i in radii){
  #     tmp <- ifelse(dat$windpark_dist>i, i, dat$windpark_dist)
  #     tmp2 <- ifelse(dat$windpark_activity=="none", i, tmp)
  #     dat <- cbind(dat, tmp, tmp2)
  #     names(dat)[c((ncol(dat)-1):ncol(dat))] <- c(paste0("windpark_all_", i), paste0("windpark_", i))#
  #   }
  #   for (i in radii){
  #     tmp <- ifelse(dat$windpark_dist>i, i, dat$windpark_dist)
  #     tmp3 <- ifelse(tmp==i, "none", dat$windpark_activity)
  #     dat <- cbind(dat, tmp3)
  #     names(dat)[ncol(dat)] <- paste0("windpark_activity_", i)#
  #   }
  # }

  #--
  # tourism

  # cabins private - nothing here
  string <- "cabins_private"
  cols_cpriv_n <- grep(string, names(dat))
  names(dat)[cols_cpriv_n]
  # cabins public
  string <- "cabins_public"
  cols_cpub_n <- grep(string, names(dat))
  names(dat)[cols_cpub_n]
  # cabins low cumulative = sum cabins small and medium
  # cabins low nearest = minimum cabins small and medium
  radii <- c(100, 250, 500, 1000, 2500, 5000, 10000)
  for (i in radii){
    # cumulative
    tmp <- dat[,paste0("cabins_public_medium_bartlett", i)] +
      dat[,paste0("cabins_public_small_bartlett", i)]
    dat <- cbind(dat, tmp)
    names(dat)[ncol(dat)] <- paste0("cabins_public_low_bartlett", i)
    # nearest
    if(!prediction | (prediction & include_zoi_nearest)) {
      tmp <- min(dat[,paste0("cabins_public_medium_nearest", i)],
                 dat[,paste0("cabins_public_small_nearest", i)], na.rm = TRUE)
      dat <- cbind(dat, tmp)
      names(dat)[ncol(dat)] <- paste0("cabins_public_low_nearest", i)
    }
  }
  # # merge WIN public cabins high-low - very few of each
  # dat$pub_cabins_win_100 <- (dat$pub_cabins_winter_high_100 + dat$pub_cabins_winter_low_100)
  # dat$pub_cabins_win_250 <- (dat$pub_cabins_winter_high_250 + dat$pub_cabins_winter_low_250)
  # dat$pub_cabins_win_500 <- (dat$pub_cabins_winter_high_500 + dat$pub_cabins_winter_low_500)
  # dat$pub_cabins_win_1000 <- (dat$pub_cabins_winter_high_1000 + dat$pub_cabins_winter_low_1000)
  # dat$pub_cabins_win_2500 <- (dat$pub_cabins_winter_high_2500 + dat$pub_cabins_winter_low_2500)
  # dat$pub_cabins_win_5000 <- (dat$pub_cabins_winter_high_5000 + dat$pub_cabins_winter_low_5000)
  # dat$pub_cabins_win_10000 <-(dat$pub_cabins_winter_high_10000 + dat$pub_cabins_winter_low_10000)
  # hotels
  string <- "hotels"
  cols_chot_n <- grep(string, names(dat))
  names(dat)[cols_chot_n]
  # cabins low cumulative = sum cabins large and hotels norway
  # cabins low nearest = minimum cabins large and hotels in norway
  radii <- c(100, 250, 500, 1000, 2500, 5000, 10000)
  for (i in radii){
    # cumulative
    tmp <- dat[,paste0("cabins_public_large_bartlett", i)] +
      dat[,paste0("hotels_bartlett", i)]
    dat <- cbind(dat, tmp)
    names(dat)[ncol(dat)] <- paste0("cabins_public_high_bartlett", i)
    # nearest
    tmp <- min(dat[,paste0("cabins_public_large_bartlett", i)],
               dat[,paste0("hotels_bartlett", i)], na.rm = TRUE)
    dat <- cbind(dat, tmp)
    names(dat)[ncol(dat)] <- paste0("cabins_public_high_nearest", i)
  }
  # summer trails
  string <- "trails"
  cols_tr_n <- grep(string, names(dat))
  names(dat)[cols_tr_n]
  names(dat)[cols_tr_n] <- sub("trails_summer", "trails", names(dat)[cols_tr_n])
  # log pseudotui
  radii <- c(100, 250, 500, 1000, 2500, 5000, 10000)
  for (i in radii){
    tmp <- log(dat[,paste0("trails_pseudotui_bartlett", i)]+1)
    dat <- cbind(dat, tmp)
    names(dat)[ncol(dat)] <- paste0("trails_log_pseudotui_bartlett", i)#
  }
  # ski tracks
  string <- "skitracks_high"
  cols_ski_n <- grep(string, names(dat))
  names(dat)[cols_ski_n]
  ski_seas <- ifelse(season == "win", "win", "cal")
  names(dat)[cols_ski_n] <- sub(paste0(string, "_", ski_seas), string, names(dat)[cols_ski_n])
  string <- "skitracks_low"
  cols_skil_n <- grep(string, names(dat))
  names(dat)[cols_skil_n]
  ski_seas <- ifelse(season == "win", "win", "cal")
  names(dat)[cols_skil_n] <- sub(paste0(string, "_", ski_seas), string, names(dat)[cols_skil_n])

  ## merge skitracks_win_      high-low - very few of each
  # dat$ski_win_100 <- (dat$skitracks_win_high_100 + dat$skitracks_win_low_100)
  # dat$ski_win_250 <- (dat$skitracks_win_high_250 + dat$skitracks_win_low_250)
  # dat$ski_win_500 <- (dat$skitracks_win_high_500 + dat$skitracks_win_low_500)
  # dat$ski_win_1000 <- (dat$skitracks_win_high_1000 + dat$skitracks_win_low_1000)
  # dat$ski_win_2500 <- (dat$skitracks_win_high_2500 + dat$skitracks_win_low_2500)
  # dat$ski_win_5000 <- (dat$skitracks_win_high_5000 + dat$skitracks_win_low_5000)
  # dat$ski_win_10000 <- (dat$skitracks_win_high_10000 + dat$skitracks_win_low_10000)
  #
  # ## merge skitracks_cal_      high-low - very few of each
  # dat$ski_cal_100 <- (dat$skitracks_cal_high_100 + dat$skitracks_cal_low_100)
  # dat$ski_cal_250 <- (dat$skitracks_cal_high_250 + dat$skitracks_cal_low_250)
  # dat$ski_cal_500 <- (dat$skitracks_cal_high_500 + dat$skitracks_cal_low_500)
  # dat$ski_cal_1000 <- (dat$skitracks_cal_high_1000 + dat$skitracks_cal_low_1000)
  # dat$ski_cal_2500 <- (dat$skitracks_cal_high_2500 + dat$skitracks_cal_low_2500)
  # dat$ski_cal_5000 <- (dat$skitracks_cal_high_5000 + dat$skitracks_cal_low_5000)
  # dat$ski_cal_10000 <- (dat$skitracks_cal_high_10000 + dat$skitracks_cal_low_10000)

  #--
  # landscape

  # PCAs
  names(dat)[names(dat)=="Norway_PCA_klima_axis1"] <- "norway_pca1"
  names(dat)[names(dat)=="Norway_PCA_klima_axis2"] <- "norway_pca2"
  names(dat)[names(dat)=="Norway_PCA_klima_axis3"] <- "norway_pca3"
  names(dat)[names(dat)=="Norway_PCA_klima_axis4"] <- "norway_pca4"
  # dem
  names(dat)[names(dat)=="dem_slope"] <- "slope"
  names(dat)[names(dat)=="dem_aspect"] <- "aspect"
  # TPI
  string <- "dem_tpi"
  cols_tpi_n <- grep(string, names(dat))
  names(dat)[cols_tpi_n]
  names(dat)[cols_tpi_n] <- sub(string, "tpi", names(dat)[cols_tpi_n])
  # for(xx in cols_tpi_n) {
  #
  # }
  # valley depth
  names(dat)[names(dat)=="valley_depth_100m"] <- "valley_depth"
  # solar radiation
  string <- "solar_radiation"
  cols_sol_n <- grep(string, names(dat))
  names(dat)[cols_sol_n]
  solar_seas <- ifelse(season == "win", "winter", ifelse(season == "sum", "summer", "spring"))
  names(dat)[cols_sol_n] <- sub(paste0(string, "_", solar_seas), string, names(dat)[cols_sol_n])
  # lakes - nothing here
  if(use_binary_vars_as_categorical) {
    string <- "lakes"
    cols_lak_n <- grep(string, names(dat))
    names(dat)[cols_lak_n]
    dat[names(dat)[cols_lak_n]] <- factor(dat[[names(dat)[cols_lak_n]]])
  }
  # set biomass to zero for missing data:
  dat$digestible_biomass_summer[is.na(dat$digestible_biomass_summer)] <- 0
  dat$digestible_biomass_winter[is.na(dat$digestible_biomass_winter)] <- 0
  # lichen
  names(dat)[grep(paste0(prefix, "lichen"), names(dat))]
  names(dat)[grep(paste0(prefix, "lichen_wrein"), names(dat))] <- paste0(prefix, "lichen_nina")
  names(dat)[grep(paste0(prefix, "lichen_model_se"), names(dat))] <- paste0(prefix, "lichen_slu")
  names(dat)[grep(paste0(prefix, "lichen_model_975"), names(dat))] <- paste0(prefix, "lichen_nose")

  # get name of classes for land cover layers
  if(prediction) {
    # read classes
    classes_norut <- read.csv("data/classes_norut_prodchange.csv", header = FALSE)
    classes_norut_smd <- read.csv("data/classes_norut_smd_oneimpact.csv", header = FALSE)

    #-- NORUT
    # replace by names, as a factor
    dat$landcover_NORUT_2009 <- ifelse(dat$landcover_NORUT_2009 == 0, NA, dat$landcover_NORUT_2009)
    dat$landcover_NORUT_2009 <- classes_norut$V3[factor(dat$landcover_NORUT_2009)] |>
      factor(levels = unique(classes_norut$V3))

    dat <- dat |>
      dplyr::mutate(lc_no_norut = landcover_NORUT_2009) |>
      fastDummies::dummy_cols(select_columns = "lc_no_norut", ignore_na = TRUE, remove_selected_columns = TRUE)

    #-- NORUT-SMD
    # replace by names, as a factor
    dat$landcover_norut_smd <- ifelse(dat$landcover_norut_smd == 0, NA, dat$landcover_norut_smd)
    dat$landcover_norut_smd <- classes_norut_smd$V3[factor(dat$landcover_norut_smd)] |>
      factor(levels = unique(classes_norut_smd$V3))

    dat <- dat |>
      dplyr::mutate(lc_nose_norut_smd = landcover_norut_smd) |>
      fastDummies::dummy_cols(select_columns = "lc_nose_norut_smd", ignore_na = TRUE, remove_selected_columns = TRUE)

    if(species == "trein") {
      levels_norut_smd <- levels(dat[[paste0(prefix, "landcover_norut_smd")]])
      clear_cuts_lines <- which(!is.na(dat[[paste0(prefix, "clear_cuts_se")]]))
      dat[[paste0(prefix, "landcover_norut_smd")]] <- as.character(dat[[paste0(prefix, "landcover_norut_smd")]])
      dat[[paste0(prefix, "landcover_norut_smd")]][clear_cuts_lines] <-
        ifelse(dat[["year"]][clear_cuts_lines] > dat[[paste0(prefix, "clear_cuts_se")]][clear_cuts_lines],
               "clear_cut", dat[[paste0(prefix, "landcover_norut_smd")]][clear_cuts_lines])
      dat[[paste0(prefix, "landcover_norut_smd")]] <- factor(dat[[paste0(prefix, "landcover_norut_smd")]],
                                                             levels = levels_norut_smd)
    }
  }

  # land cover
  if(land_cover_factor) {

    if(land_cover == "norut_smd") {
      string <- "lc_nose_norut_smd"
      cols_lc_n <- grep(string, names(dat))
      names(dat)[cols_lc_n]
      names(dat)[cols_lc_n] <- sub(string, "lc", names(dat)[cols_lc_n])

      dat$landcover_norut_smd <- relevel(dat$landcover_norut_smd, ref = ref_landcover)
      # remove lichen to the reference category
      # another possibility would be to move it to "other"
      #------- NEED TO BE CHANGED in case the ref category is not heathland
      dat$landcover_norut_smd <- forcats::fct_collapse(dat$landcover_norut_smd, heathland = c("heathland", "lichen"))
    } else {

      if(land_cover == "norut") {
        string <- "lc_no_norut"
        cols_lc_n <- grep(string, names(dat))
        names(dat)[cols_lc_n]
        names(dat)[cols_lc_n] <- sub(string, "lc", names(dat)[cols_lc_n])

        dat$landcover_NORUT_2009 <- relevel(dat$landcover_NORUT_2009, ref = ref_landcover)
        # remove lichen to the reference category
        # another possibility would be to move it to "other"
        #------- NEED TO BE CHANGED in case the ref category is not heathland
        dat$landcover_NORUT_2009 <- forcats::fct_collapse(dat$landcover_NORUT_2009, heathland = c("heathland", "lichen"))
      }

      #--- add options for SMD and NMD later for Sweden
    }

  } else {

    string <- "lc_nose_norut_smd"
    cols_lc_n <- grep(string, names(dat))
    names(dat)[cols_lc_n]
    names(dat)[cols_lc_n] <- sub(string, "lc", names(dat)[cols_lc_n])

    # remove lakes and reservoirs from water
    # dat[["lc_water_no_lak_res"]] <- ifelse(dat[["lc_water"]] == dat[["reservoirs_hydro"]] |
    #                                          dat[["lc_water"]] == dat[["lakes"]], 0, dat[["lc_water"]])

  }


  #--------------------------
  # add option to use LC as factors instead (is there any difference?)
  # using the tables for reclassifying from numbers to strings

  #--
  # climate phenology

  # onset of spring
  names(dat)[names(dat)=="sprAverageAllYrs"] <- "onset_spring"
  # snow cover duration (days)
  names(dat)[grep(paste0(prefix, "snow_cover_duration"), names(dat))] <- paste0(prefix, "snow_cover_days")
  # bio12 annual prec
  names(dat)[names(dat)=="BIOCLIM_bio12_historic"] <- "bio12_annual_prec"
  # bio18 warmest prec
  names(dat)[names(dat)=="BIOCLIM_bio18_historic"] <- "bio18_warmest_prec"
  # bio19 coldest prec
  names(dat)[names(dat)=="BIOCLIM_bio19_historic"] <- "bio19_coldest_prec"
  # ENVIREM continentality
  names(dat)[names(dat)=="ENVIREM_continentality_historic"] <- "continentality"
  # ENVIREM growing degree days 5C
  names(dat)[names(dat)=="ENVIREM_growingDegDays5_historic"] <- "growing_deg_days5"
  # CHELSA snow cover days
  names(dat)[names(dat)=="CHELSA_scd_historic"] <- "chelsa_scd"

  dat <- dat |>
    dplyr::mutate(
      tpi_500_small = tpi_500,
      tpi_500_small_2 = tpi_500**2,
      tpi_1000_2 = tpi_1000**2,
      tpi_5000_large = tpi_5000,
      tpi_5000_large_2 = tpi_5000**2,
      tpi_250_small = tpi_250,
      tpi_250_small_2 = tpi_250**2,
      tpi_2500_large = tpi_2500,
      tpi_2500_large_2 = tpi_2500**2,
      onset_spring_2 = onset_spring**2,
      snow_cover_days_2 = snow_cover_days**2,
      slope_2 = slope**2,
      solar_radiation_2 = solar_radiation**2,
      # bio12_annual_prec_2 = bio12_annual_prec**2,
      bio18_warmest_prec_2 = bio18_warmest_prec**2,
      bio19_coldest_prec_2 = bio19_coldest_prec**2,
      continentality_2 = continentality**2,
      growing_deg_days5_2 = growing_deg_days5**2,
      chelsa_scd_2 = chelsa_scd**2)

  #--
  # species
  #--------------- CHECK HERE!!
  # for prediction
  # compute average density of grazing animals
  # average grazing animals
  # if(prediction) {
  #   string <- "grazing_animals_\\d{4}_bartlett"
  #   cols_graz_n <- grep(string, names(dat))
  #   names(dat)[cols_graz_n]
  #   radii <- c(100, 250, 500, 1000, 2500, 5000, 10000)
  #   rr <- radii[1]
  #   for(rr in radii) {
  #     vars <- grep(paste0("bartlett", rr, "$"), names(dat)[cols_graz_n], value = TRUE)
  #     new_var <- sub("\\d{4}", "avg_yrs", vars[1])
  #     dat[[new_var]] <- apply(dat[, match(vars, colnames(dat))], 1, mean, na.rm = TRUE)
  #     dat <- dat[, -match(vars, colnames(dat))]
  #   }
  # }
  # rename vars
  string <- "grazing_animals_avg_yrs_bartlett"
  cols_graz_avg_n <- grep(string, names(dat))
  names(dat)[cols_graz_avg_n]
  names(dat)[cols_graz_avg_n] <- sub(string, "grazing_animals_avg", names(dat)[cols_graz_avg_n])


  # # merge roads  high-low - very few of each
  # dat$roads_100 <- (dat$roads_summer_high_100 + dat$roads_summer_low_100)
  # dat$roads_250 <- (dat$roads_summer_high_250 + dat$roads_summer_low_250)
  # dat$roads_500 <- (dat$roads_summer_high_500 + dat$roads_summer_low_500)
  # dat$roads_1000 <- (dat$roads_summer_high_1000 + dat$roads_summer_low_1000)
  # dat$roads_2500 <- (dat$roads_summer_high_2500 + dat$roads_summer_low_2500)
  # dat$roads_5000 <- (dat$roads_summer_high_5000 + dat$roads_summer_low_5000)
  # dat$roads_10000 <- (dat$roads_summer_high_10000 + dat$roads_summer_low_10000)
  #
  # # merge WINTER roads  high-low - very few of each
  # dat$roads_win_100 <- (dat$roads_winter_high_100 + dat$roads_winter_low_100)
  # dat$roads_win_250 <- (dat$roads_winter_high_250 + dat$roads_winter_low_250)
  # dat$roads_win_500 <- (dat$roads_winter_high_500 + dat$roads_winter_low_500)
  # dat$roads_win_1000 <- (dat$roads_winter_high_1000 + dat$roads_winter_low_1000)
  # dat$roads_win_2500 <- (dat$roads_winter_high_2500 + dat$roads_winter_low_2500)
  # dat$roads_win_5000 <- (dat$roads_winter_high_5000 + dat$roads_winter_low_5000)
  # dat$roads_win_10000 <- (dat$roads_winter_high_10000 + dat$roads_winter_low_10000)
  #
  # # Powerlines without roads (variable name: powerlines_noroad_XXXX)
  # radii <- c(100, 250, 500, 1000, 2500, 5000, 10000)
  # for (rad in radii){
  #   dat[,c(paste0("powerlines_noroad_", rad))] <- dat[,c(paste0("powerlines_", rad))]-dat[,c(paste0("roads_win_", rad))]
  #   dat[,c(paste0("powerlines_noroad_", rad))] <- ifelse(dat[,c(paste0("powerlines_noroad_", rad))]<0, 0, dat[,c(paste0("powerlines_noroad_", rad))])
  # }
  #
  # # Powerlines without roads AND railway (variable name: powerlines_noroadrail_XXXX)
  # radii <- c(100, 250, 500, 1000, 2500, 5000, 10000)
  # for (rad in radii){
  #   dat[,c(paste0("powerlines_noroadrail_", rad))] <- dat[,c(paste0("powerlines_noroad_", rad))]-dat[,c(paste0("railway_", rad))]
  #   dat[,c(paste0("powerlines_noroadrail_", rad))] <- ifelse(dat[,c(paste0("powerlines_noroadrail_", rad))]<0, 0, dat[,c(paste0("powerlines_noroadrail_", rad))])
  # }

  #-------------
  # Add dummy variables for prediction, in the same way we did not the input data
  # add options for norut, norut_smd, smd, and nmd

  # add NORUT columns (dummy variable) - before we had density at 100 m - now we go for simple pixel-scale var (#Separate NORUT recklass into separate var)
  # if (prediction){
  #   reclass_norut <- function(x){
  #     y <- NA
  #     y <- ifelse((x>=1 & x<=2)|(x>=4 & x<=7), "skog", y)
  #     y <- ifelse(x==3 | x==8, "lavskog",y)
  #     y <- ifelse(x>=9 & x<=11, "myr",y)
  #     y <- ifelse(x==12, "ridges",y)
  #     y <- ifelse(x==13, "grasses",y)
  #     y <- ifelse(x==14, "heather_ridges",y)
  #     y <- ifelse(x==15, "lichen",y)
  #     y <- ifelse(x==16, "heather_lowland",y)
  #     y <- ifelse(x==17, "heathland",y)
  #     y <- ifelse(x==18, "meadows",y)
  #     y <- ifelse(x==19, "snowbed",y)
  #     y <- ifelse(x==20, "snow",y)
  #     y <- ifelse(x==21, "glacier",y)
  #     y <- ifelse(x==23, "dyrka",y)
  #     y <- ifelse(x==22, "water",y)
  #     y <- ifelse(x==24, "urban",y)
  #     y <- ifelse(x==25, "noclass",y)
  #     return(y)
  #   }
  #
  #   dat$norut_reclass <- reclass_norut(dat$norut)
  # }
  #
  # norut_columns <- function(x){
  #   y <- data.frame(dyrka=as.numeric(x=="dyrka"),
  #                   glacier=as.numeric(x=="glacier"),
  #                   grasses=as.numeric(x=="grasses"),
  #                   heather_lowland=as.numeric(x=="heather_lowland"),
  #                   heather_ridges=as.numeric(x=="heather_ridges"),
  #                   heathland=as.numeric(x=="heathland"),
  #                   lavskog=as.numeric(x=="lavskog"),
  #                   lichen=as.numeric(x=="lichen"),
  #                   meadows=as.numeric(x=="meadows"),
  #                   myr=as.numeric(x=="myr"),
  #                   ridges=as.numeric(x=="ridges"),
  #                   skog=as.numeric(x=="skog"),
  #                   snow=as.numeric(x=="snow"),
  #                   snowbed=as.numeric(x=="snowbed"),
  #                   noclass=as.numeric(x=="noclass"),
  #                   urban=as.numeric(x=="urban"),
  #                   water=as.numeric(x=="water"))
  #   return(y)
  # }
  # dat <- cbind(dat, norut_columns(dat$norut_reclass))
  # dat$norut_reclass <- NULL
  #
  # if (use_norut_prop){
  #   dat$dyrka <- dat$norut_dyrka_100
  #   dat$glacier <- dat$norut_glacier_100
  #   dat$grasses <- dat$norut_grasses_100
  #   dat$heather_lowland <- dat$norut_heather_lowland_100
  #   dat$heather_ridges <- dat$norut_heather_ridges_100
  #   dat$heathland <- dat$norut_heathland_100
  #   dat$lavskog <- dat$norut_lavskog_100
  #   dat$lichen <- dat$norut_lichen_100
  #   dat$meadows <- dat$norut_meadows_100
  #   dat$myr <- dat$norut_myr_100
  #   dat$ridges <- dat$norut_ridges_100
  #   dat$skog <- dat$norut_skog_100
  #   dat$snow <- dat$norut_snow_100
  #   dat$snowbed <- dat$norut_snowbed_100
  # }

  if(prediction) {

    # fill NAs
    dat[[paste0(prefix, "onset_spring")]][is.na(dat[[paste0(prefix, "onset_spring")]])] <- max(dat[[paste0(prefix, "onset_spring")]], na.rm = TRUE)

    if(sum(grepl(paste0(prefix, "lichen_nina"), names(dat))) > 1)
      dat[[paste0(prefix, "lichen_nina")]][is.na(dat[[paste0(prefix, "lichen_nina")]])] <- mean(dat[[paste0(prefix, "lichen_nina")]], na.rm = TRUE)
    if(sum(grepl(paste0(prefix, "lichen_slu"), names(dat))) > 1)
      dat[[paste0(prefix, "lichen_slu")]][is.na(dat[[paste0(prefix, "lichen_slu")]])] <- mean(dat[[paste0(prefix, "lichen_slu")]], na.rm = TRUE)
    # dat[[paste0(prefix, "lichen_nose")]][is.na(dat[[paste0(prefix, "lichen_nose")]])] <- mean(dat[[paste0(prefix, "lichen_nose")]], na.rm = TRUE)
  }

  # identify non-numeric columns,  that shoud not be standardized, and move non-numeric to the front - so it will be easier to standardize in the ModSel files
  # cls <- NA
  # for (i in c(1:ncol(dat))){cls[i] <- class(dat[,i])}
  # nums <- which((cls %in% c("numeric","integer")))
  # nums <- nums[nums>6]
  # nonums <- which(!(cls %in% c("numeric","integer")))
  # nonums <- nonums[nonums>6]
  # dat <- cbind(dat[,c(1:6)], dat[,nonums], dat[,nums])# columns that are not numeric or integer => should NOT standardize these, in the ModSel files

  return(dat)
}
