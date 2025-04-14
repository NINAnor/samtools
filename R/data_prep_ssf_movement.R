#' Prepare data for SSF movement/permeability analyses for wild and semi-domesticated reindeer within SAM
#'
#' This function takes in the already annotated data for SSF with either wild or semi-domesticated
#' reindeer and further prepare columns and variables for the model fitting procedure.
#' The function is suited for step resource selection functions (SSF) only.
#'
#' @param dat `[data.frame]` \cr Annotated data set for analysis.
#' @param season `[string]{"sum", "cal", "win"}` \cr Season of interest.
#' One of "sum", "cal", or "win".
#' @param prediction `[logical(1)=FALSE]` \cr Additional changes that should
#' only be done in the grid for prediction.
#' @param land_cover `[string(1)="norut_smd"]{"norut_smd", "norut", "smd", "nmd"}` \cr
#' Which land cover map should be used for analysis. It prepares the corresponding
#' classes as
#' @param prefix `[string(1)="endpt_]{"endpt_", "along_", ""}` \cr Prefix for the variable,
#' whether `"endpt_"` or `"along_"` for variables measures in the end point or along
#' the step, and `""` for normal RSF or for the prediction using the grid.
#' @param include_zoi_nearest `[logical(1)=FALSE]` \cr By default, `FALSE`. It should be
#' set to `TRUE` if there are variables representing the ZOI of the nearest feature in
#' the formula. Relevant only for prediction (using the grid).
#'
#' @export
data_prep_ssf_movement_rein <- function(dat, season,
                                        prediction = FALSE,
                                        fixwind = TRUE,
                                        land_cover_factor = FALSE,
                                        land_cover = c("norut_smd", "norut", "smd", "nmd")[1],
                                        ref_landcover = "heathland",
                                        tpi_factor = FALSE,
                                        include_zoi_nearest = FALSE,
                                        species = c("wrein", "trein")[1],
                                        prefix = c("along_", "endpt_", "")[1],
                                        prefix_cross = "cross",
                                        reference_year = lubridate::year(lubridate::now()),
                                        use_binary_vars_as_categorical = FALSE,
                                        formula = NULL){

  #---
  # rename variables and compute new, derived variables

  # movement variables
  if(!prediction) {
    # we need to compute that before hand!!! relative angle
    # dat$rel_angle[is.na(dat$rel_angle)] <-
    dat$step_length <- dat$step_length + 50
    dat$log_step_length <- log10(dat$step_length)

    cross_vars <- grep("cross", names(dat))
    names(dat)[cross_vars]
    for(i in cross_vars) {
      dat[[i]][dat[[i]] < 0] <- 0
    }

  } else {
    dat$step_length <- 100
    dat$log_step_length <- log10(dat$step_length)

    if(!is.null(formula)) {
      all_vars <- all.vars(f)[-c(1,2)]
    }
  }

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
  # log road volume in roads major
  radii <- c(100, 250, 500, 1000, 2500, 5000, 10000)
  i <- radii[1]
  for (i in radii){
    cols <- grep(pattern = "bartlett", names(dat)[cols_rmaj_n], value = TRUE) |>
      grep(pattern = paste0("bartlett", i, "_"), fixed = TRUE, value = TRUE) |>
      grep(pattern = "win|sum|cal", invert = TRUE, value = TRUE)
    cc <- cols[1]
    for(cc in cols) {
      tmp <- log10(dat[,cc]+1)
      dat <- cbind(dat, tmp)
      names(dat)[ncol(dat)] <- paste0(cc, "_log")
    }
    # tail(names(dat), 20)
  }

  # roads major - cross
  string <- "cross_roads_major"
  cols_rmaj_n <- grep(string, names(dat))
  names(dat)[cols_rmaj_n]
  roads_seas <- season
  names(dat)[cols_rmaj_n] <- sub(paste0(string, "_bin_", roads_seas), string, names(dat)[cols_rmaj_n])
  if(prediction) {
    v_cross_roads <- grep(string, all_vars, value = TRUE)
    string_env <- "roads_major"
    grep(string_env, names(dat), value = TRUE)
    roads_seas <- ifelse(season == "sum", "sum", "win")
    dat[[v_cross_roads]] <- dat[[string_env]]
    #### Add a condition to turn into 0, 1 if the variable is binary
  }
  # log road volume in roads major
  cols <- grep(pattern = "win_|sum_|cal_", names(dat)[cols_rmaj_n], invert = TRUE, value = TRUE)
  cc <- cols[1]
  for(cc in cols) {
    tmp <- log10(dat[,cc]+1)
    dat <- cbind(dat, tmp)
    names(dat)[ncol(dat)] <- paste0(cc, "_log")
  }

  # roads minor
  string <- "roads_minor"
  cols_rmin_n <- grep(string, names(dat))
  names(dat)[cols_rmin_n]
  roads_seas <- ifelse(season == "sum", "sum", "win")
  names(dat)[cols_rmin_n] <- sub(paste0(string, "_", roads_seas), string, names(dat)[cols_rmin_n])

  # roads minor - cross
  string <- "cross_roads_minor"
  cols_rmin_n <- grep(string, names(dat))
  names(dat)[cols_rmin_n]
  roads_seas <- season
  names(dat)[cols_rmin_n] <- sub(paste0(string, "_bin_", roads_seas), string, names(dat)[cols_rmin_n])
  if(prediction) {
    v_cross_roads <- grep(string, all_vars, value = TRUE)
    string_env <- "roads_minor"
    grep(string_env, names(dat), value = TRUE)
    roads_seas <- ifelse(season == "sum", "sum", "win")
    dat[[v_cross_roads]] <- dat[[string_env]]
    #### Add a condition to turn into 0, 1 if the variable is binary
  }

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

  # roads minor - cross
  string <- "cross_railways"
  cols_rmin_n <- grep(string, names(dat))
  names(dat)[cols_rmin_n]
  rail_seas <- ifelse(season == "sum", "sum", "win")
  names(dat)[cols_rmin_n] <- sub(paste0(string, "_bin_", rail_seas), string, names(dat)[cols_rmin_n])
  if(prediction) {
    v_cross_rails <- grep(string, all_vars, value = TRUE)
    string_env <- "railways"
    grep(string_env, names(dat), value = TRUE)
    dat[[v_cross_rails]] <- dat[[paste0(string_env)]]
    #### Add a condition to turn into 0, 1 if the variable is binary
  }

  # fences cross
  if(prediction & species == "trein") {
    string <- "cross_fences"
    v_cross_roads <- grep(string, all_vars, value = TRUE)
    string_env <- "fences"
    grep(string_env, names(dat), value = TRUE)

    dat[[v_cross_roads]] <- dat[[string_env]]
    #### Add a condition to turn into 0, 1 if the variable is binary
  }

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

  # wind turbines -- add to trein
  string <- "wind_turbines"
  cols_rpl_n <- grep(string, names(dat))
  names(dat)[cols_rpl_n]

  # reservoirs - transform into factor
  string <- "reservoirs_hydro"
  cols_hyd_n <- grep(string, names(dat))
  nam_hyd <- names(dat)[cols_hyd_n] |>
    grep(pattern = "cross", invert = TRUE, value = TRUE) |>
    grep(pattern = "along", invert = TRUE, value = TRUE)
  if(length(nam_hyd) > 0) {
    if(prefix == "endpt_") {
      dat[[nam_hyd]] <- factor(dat[[nam_hyd]]) ######## NOT FOR ALONG???
    }
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
    if(prediction) {

      var_cab_low <- grep("cabins_public_low_bartlett", all.vars(f), value = T)
      summ <- strsplit(var_cab_low[1], split = "_")[[1]] |> dplyr::last()

      grep("cabins_public_medium_bartlett", names(dat), value = T)
      # if cumulative
      if(grepl("bartlett", var_cab_low[1])) {
        if(summ == "mean") {
          tmp <- dat[,paste0("cabins_public_medium_bartlett", i)] +
            dat[,paste0("cabins_public_small_bartlett", i)]
          dat <- cbind(dat, tmp)
          names(dat)[ncol(dat)] <- paste0(prefix, "cabins_public_low_bartlett", i, "_", summ)
        } else {
          if(summ == "max") {
            tmp <- max(dat[,paste0("cabins_public_medium_bartlett", i)],
                       dat[,paste0("cabins_public_small_bartlett", i)], na.rm = TRUE)
            dat <- cbind(dat, tmp)
            names(dat)[ncol(dat)] <- paste0(prefix, "cabins_public_low_bartlett", i, "_", summ)
          }
        }


      } else {

        if(grepl("nearest", var_cab_low[1])) {
          if(summ == "mean") {
            tmp <- dat[,paste0("cabins_public_medium_nearest", i)] +
              dat[,paste0("cabins_public_small_nearest", i)]
            dat <- cbind(dat, tmp)
            names(dat)[ncol(dat)] <- paste0(prefix, "cabins_public_low_nearest", i, "_", summ)
          } else {
            if(summ == "max") {
              tmp <- min(dat[,paste0("cabins_public_medium_nearest", i)],
                         dat[,paste0("cabins_public_small_nearest", i)], na.rm = TRUE)
              dat <- cbind(dat, tmp)
              names(dat)[ncol(dat)] <- paste0(prefix, "cabins_public_low_nearest", i, "_", summ)
            }
          }
        }


      }

    } else {

      for(summ in c("_max", "_mean")) {
        if(summ == "_mean") {
          tmp <- dat[,paste0(prefix, "cabins_public_medium_bartlett", i, summ)] +
            dat[,paste0(prefix, "cabins_public_small_bartlett", i, summ)]
          dat <- cbind(dat, tmp)
          names(dat)[ncol(dat)] <- paste0(prefix, "cabins_public_low_bartlett", i, summ)
        } else {
          if(summ == "_max") {
            tmp <- max(dat[,paste0(prefix, "cabins_public_medium_bartlett", i, summ)],
                       dat[,paste0(prefix, "cabins_public_small_bartlett", i, summ)], na.rm = TRUE)
            dat <- cbind(dat, tmp)
            names(dat)[ncol(dat)] <- paste0(prefix, "cabins_public_low_bartlett", i, summ)
          }
        }
      }

      # nearest
      if(!prediction) {
        for(summ in c("_max", "_mean")) {
          if(summ == "_mean") {
            tmp <- dat[,paste0(prefix, "cabins_public_medium_nearest", i, summ)] +
              dat[,paste0(prefix, "cabins_public_small_nearest", i, summ)]
            dat <- cbind(dat, tmp)
            names(dat)[ncol(dat)] <- paste0(prefix, "cabins_public_low_nearest", i, summ)
          } else {
            if(summ == "_max") {
              tmp <- min(dat[,paste0(prefix, "cabins_public_medium_nearest", i, summ)],
                         dat[,paste0(prefix, "cabins_public_small_nearest", i, summ)], na.rm = TRUE)
              dat <- cbind(dat, tmp)
              names(dat)[ncol(dat)] <- paste0(prefix, "cabins_public_low_nearest", i, summ)
            }
          }
        }
      }

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

    if(prediction) {

      var_cab_low <- grep("cabins_public_high_bartlett", all.vars(f), value = T)
      summ <- strsplit(var_cab_low[1], split = "_")[[1]] |> dplyr::last()
      # if cumulative
      if(grepl("bartlett", var_cab_low[1])) {
        if(summ == "mean") {
          tmp <- dat[,paste0("cabins_public_large_bartlett", i)] +
            dat[,paste0("hotels_bartlett", i)]
          dat <- cbind(dat, tmp)
          names(dat)[ncol(dat)] <- paste0(prefix, "cabins_public_high_bartlett", i, "_", summ)
        } else {
          if(summ == "max") {
            tmp <- max(dat[,paste0("cabins_public_large_bartlett", i)],
                       dat[,paste0("hotels_bartlett", i)], na.rm = TRUE)
            dat <- cbind(dat, tmp)
            names(dat)[ncol(dat)] <- paste0(prefix, "cabins_public_high_bartlett", i, "_", summ)
          }
        }


      } else {

        if(grepl("nearest", var_cab_low[1])) {
          if(summ == "mean") {
            tmp <- dat[,paste0("cabins_public_large_nearest", i)] +
              dat[,paste0("hotels_nearest", i)]
            dat <- cbind(dat, tmp)
            names(dat)[ncol(dat)] <- paste0(prefix, "cabins_public_high_nearest", i, "_", summ)
          } else {
            if(summ == "max") {
              tmp <- min(dat[,paste0("cabins_public_large_nearest", i)],
                         dat[,paste0("hotels_nearest", i)], na.rm = TRUE)
              dat <- cbind(dat, tmp)
              names(dat)[ncol(dat)] <- paste0(prefix, "cabins_public_high_nearest", i, "_", summ)
            }
          }
        }


      }

    } else {

      # cumulative
      for(summ in c("_max", "_mean")) {
        if(summ == "_mean") {
          tmp <- dat[,paste0(prefix, "cabins_public_large_bartlett", i, summ)] +
            dat[,paste0(prefix, "hotels_bartlett", i, summ)]
          dat <- cbind(dat, tmp)
          names(dat)[ncol(dat)] <- paste0(prefix, "cabins_public_high_bartlett", i, summ)
        } else {
          if(summ == "_max") {
            tmp <- max(dat[,paste0(prefix, "cabins_public_large_bartlett", i, summ)],
                       dat[,paste0(prefix, "hotels_bartlett", i, summ)], na.rm = TRUE)
            dat <- cbind(dat, tmp)
            names(dat)[ncol(dat)] <- paste0(prefix, "cabins_public_high_bartlett", i, summ)
          }
        }
      }

      # nearest
      if(!prediction & include_zoi_nearest) {
        for(summ in c("_max", "_mean")) {
          if(summ == "_mean") {
            tmp <- dat[,paste0(prefix, "cabins_public_large_nearest", i, summ)] +
              dat[,paste0(prefix, "hotels_nearest", i, summ)]
            dat <- cbind(dat, tmp)
            names(dat)[ncol(dat)] <- paste0(prefix, "cabins_public_high_nearest", i, summ)
          } else {
            if(summ == "_max") {
              tmp <- min(dat[,paste0(prefix, "cabins_public_large_nearest", i, summ)],
                         dat[,paste0(prefix, "hotels_nearest", i, summ)], na.rm = TRUE)
              dat <- cbind(dat, tmp)
              names(dat)[ncol(dat)] <- paste0(prefix, "cabins_public_high_nearest", i, summ)
            }
          }
        }
      }

    }

  }

  # summer trails
  if(season == "sum") {
    string <- "trails"
    cols_tr_n <- grep(string, names(dat))
    names(dat)[cols_tr_n]
    names(dat)[cols_tr_n] <- sub("trails_summer", "trails", names(dat)[cols_tr_n])
    # log pseudotui
    # radii <- c(100, 250, 500, 1000, 2500, 5000, 10000)
    # for (i in radii){
    #   tmp <- log(dat[,paste0(prefix, "trails_pseudotui_bartlett", i)]+1) ## correct here - check below
    #   dat <- cbind(dat, tmp)
    #   names(dat)[ncol(dat)] <- paste0(prefix, "trails_log_pseudotui_bartlett", i)#
    # }
    cols_trails <- grep("cross|euclidean|nearest", names(dat)[cols_tr_n], value = TRUE, invert = TRUE)
    cols_tr_n <- grep(paste0(cols_trails, collapse = "|"), names(dat))
    names(dat)[cols_tr_n]
    for(i in cols_tr_n) {
      tmp <- log10(dat[, i]+1)
      dat <- cbind(dat, tmp)
      names(dat)[ncol(dat)] <- sub("pseudotui", "log_pseudotui", names(dat)[i])
    }

    if(prediction) {
      string <- "cross_trails"
      v_cross_trails <- grep(string, all_vars, value = TRUE)
      string_env <- "trails"
      grep(string_env, names(dat), value = TRUE)
      dat[[v_cross_trails]] <- dat[["trails_log_pseudotui_bartlett100"]]
      #### chnge later when we annotate the pseudotui correctly here, compute log - like here below

    } else {
      # summer trails - log cross
      string <- "cross_trails"
      cols_tr_n <- grep(string, names(dat))
      names(dat)[cols_tr_n]
      for(i in cols_tr_n) {
        tmp <- log10(dat[, i]+1)
        dat <- cbind(dat, tmp)
        names(dat)[ncol(dat)] <- sub("pseudotui", "log_pseudotui", names(dat)[i])
      }
    }


  } else {

    # ski tracks
    string <- "skitracks_high"
    cols_ski_n <- grep(string, names(dat))
    names(dat)[cols_ski_n]
    ski_seas <- ifelse(season == "win", "win", "cal")
    names(dat)[cols_ski_n] <- sub(paste0(string, "_", ski_seas), string, names(dat)[cols_ski_n])
    names(dat)[cols_ski_n] <- sub(paste0(string, "_bin_", ski_seas), string, names(dat)[cols_ski_n])
    string <- "skitracks_low"
    cols_skil_n <- grep(string, names(dat))
    names(dat)[cols_skil_n]
    ski_seas <- ifelse(season == "win", "win", "cal")
    names(dat)[cols_skil_n] <- sub(paste0(string, "_", ski_seas), string, names(dat)[cols_skil_n])
    names(dat)[cols_skil_n] <- sub(paste0(string, "_bin_", ski_seas), string, names(dat)[cols_skil_n])

    # crossing ski tracks - add HERE
    if(prediction) {
      # skitracks high
      string <- "cross_skitracks_high"
      v_cross_ski <- grep(string, all_vars, value = TRUE)
      string_env <- "skitracks_high"
      grep(string_env, names(dat), value = TRUE)
      dat[[v_cross_ski]] <- dat[["skitracks_high"]]

      # skitracks low
      string <- "cross_skitracks_low"
      v_cross_skil <- grep(string, all_vars, value = TRUE)
      string_env <- "skitracks_low"
      grep(string_env, names(dat), value = TRUE)
      dat[[v_cross_skil]] <- dat[["skitracks_low"]]
      #### chnge later when we annotate the pseudotui correctly here, compute log - like here below

    }

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

  }



  #--
  # landscape

  # PCAs
  names(dat)[grep("Norway_PCA_klima_axis1", names(dat))] <- sub("Norway_PCA_klima_axis1", "norway_pca1", names(dat)[grep("Norway_PCA_klima_axis1", names(dat))])
  names(dat)[grep("Norway_PCA_klima_axis2", names(dat))] <- sub("Norway_PCA_klima_axis2", "norway_pca2", names(dat)[grep("Norway_PCA_klima_axis2", names(dat))])
  names(dat)[grep("Norway_PCA_klima_axis3", names(dat))] <- sub("Norway_PCA_klima_axis3", "norway_pca3", names(dat)[grep("Norway_PCA_klima_axis3", names(dat))])
  names(dat)[grep("Norway_PCA_klima_axis4", names(dat))] <- sub("Norway_PCA_klima_axis4", "norway_pca4", names(dat)[grep("Norway_PCA_klima_axis4", names(dat))])

  # dem
  if(!prediction) {
    if(grepl("along", prefix)) {
      dat$along_dem_range <- abs(dat$along_dem_max - dat$along_dem_min)
      dat$dem_diff_start_end <- abs(dat$endpt_dem_10m - dat$startpt_dem)
      dat$dem_diff_range_startend <- dat$along_dem_range - dat$dem_diff_start_end
      # dat[c("startpt_dem", "endpt_dem_10m", "dem_diff_start_end", "along_dem_range", "along_dem_cum_diff",  "dem_diff_range_startend")]
      # dat$
      # dat$along_dem_cum_diff
    }
    dat$log_along_dem_cum_diff <- log10(dat$along_dem_cum_diff + 1)
  }


  names(dat)[grep("dem_slope", names(dat))] <- sub("dem_slope", "slope", names(dat)[grep("dem_slope", names(dat))])
  names(dat)[grep("dem_aspect", names(dat))] <- sub("dem_aspect", "aspect", names(dat)[grep("dem_aspect", names(dat))])
  # TPI
  string <- "dem_tpi"
  cols_tpi_n <- grep(string, names(dat))
  names(dat)[cols_tpi_n]
  names(dat)[cols_tpi_n] <- sub(string, "tpi", names(dat)[cols_tpi_n])
  # TPI - categories
  if(tpi_factor) {

    tpi_scales <- c(250, 2500)
    string <- "tpi"
    cols_tpi_n <- grep(string, names(dat))
    tpi_sc <- tpi_scales[1]
    for(tpi_sc in tpi_scales) {

      if(!prediction) {
        cols_tpi_n <- grep(paste0("tpi_", tpi_sc, "_mean"), names(dat))
      } else {
        cols_tpi_n <- grep(paste0("tpi_", tpi_sc, "$"), names(dat))
      }

      # this should be the same SD for all the data!! we need to store that in the bag
      SD <- sd(dat[[names(dat)[cols_tpi_n]]], na.rm = TRUE)
      new_col_name <- paste0(names(dat)[cols_tpi_n], "_cat")
      dat[[new_col_name]] <- dplyr::case_when(
        dat[[names(dat)[cols_tpi_n]]] <= -SD[1] ~ "valley",
        dat[[names(dat)[cols_tpi_n]]] > -SD[1] & dat[[names(dat)[cols_tpi_n]]] <= -SD[1]/2  ~ "lower_slope",
        dat[[names(dat)[cols_tpi_n]]] > -SD[1]/2 & dat[[names(dat)[cols_tpi_n]]] <= 0  ~ "flat_area",
        dat[[names(dat)[cols_tpi_n]]] > 0 & dat[[names(dat)[cols_tpi_n]]] <= SD[1]/2  ~ "medium_slope",
        dat[[names(dat)[cols_tpi_n]]] > SD[1]/2 & dat[[names(dat)[cols_tpi_n]]] <= SD[1]  ~ "upper_slope",
        TRUE ~ "ridge"
      )
      rm(SD)
      dat <- dat |>

        fastDummies::dummy_cols(select_columns = paste0(new_col_name), ignore_na = TRUE, remove_selected_columns = FALSE)

    }
    # dat$along_dem_tpi_2500_mean

  }

  # valley depth
  names(dat)[grep("valley_depth_100m", names(dat))] <- sub("valley_depth_100m", "valley_depth", names(dat)[grep("valley_depth_100m", names(dat))])
  # solar radiation
  string <- "solar_radiation"
  cols_sol_n <- grep(string, names(dat))
  names(dat)[cols_sol_n]
  solar_seas <- ifelse(season == "win", "winter", ifelse(season == "sum", "summer", "spring"))
  names(dat)[cols_sol_n] <- sub(paste0(string, "_", solar_seas), string, names(dat)[cols_sol_n])
  # lakes - nothing here
  string <- "lakes"
  cols_lak_n <- grep(string, names(dat))
  names(dat)[cols_lak_n]
  nam_lak <- names(dat)[cols_lak_n] |>
    grep(pattern = "cross", invert = TRUE, value = TRUE) |>
    grep(pattern = "along", invert = TRUE, value = TRUE)
  if(length(nam_lak) > 0) {
    if(prefix == "endpt_") {
      dat[[nam_lak]] <- factor(dat[[nam_lak]]) ########## IF ENDPT, not ALONG?
    }
  }

  # set biomass to zero for missing data:
  string <- "digestible_biomass"
  cols_biomass_n <- grep(string, names(dat))
  for(i in cols_biomass_n) {
    dat[,i][is.na(dat[,i])] <- 0
  }

  # lichen
  #### add condition for a different layer if tamrein in norway or sweden
  names(dat)[grep(paste0(prefix, "lichen"), names(dat))]
  names(dat)[grep(paste0(prefix, "lichen_wrein"), names(dat))] <- sub("lichen_wrein_nina", "lichen_nina", names(dat)[grep(paste0(prefix, "lichen_wrein"), names(dat))])
  names(dat)[grep(paste0(prefix, "lichen_nina_no"), names(dat))] <- sub("lichen_nina_no", "lichen_nina", names(dat)[grep(paste0(prefix, "lichen_nina_no"), names(dat))])#lichen_nina_no
  names(dat)[grep(paste0(prefix, "lichen_model_se"), names(dat))] <- sub("lichen_model_se", "lichen_slu", names(dat)[grep(paste0(prefix, "lichen_model_se"), names(dat))])
  names(dat)[grep(paste0(prefix, "lichen_model_975"), names(dat))] <- sub("lichen_model_975", "lichen_nose", names(dat)[grep(paste0(prefix, "lichen_model_975"), names(dat))])
  names(dat)[grep(paste0(prefix, "lichen_nose_01"), names(dat))] <- sub("lichen_nose_01", "lichen_nose", names(dat)[grep(paste0(prefix, "lichen_nose_01"), names(dat))])#lichen_nina_no

  # dat[[grep(paste0(prefix, "lichen_nina"), names(dat))]][is.na(dat[grep(paste0(prefix, "lichen_nina"), names(dat))])] <- mean(dat[[grep(paste0(prefix, "lichen_nina"), names(dat))]], na.rm = TRUE)
  # get name of classes for land cover layers
  ############## CHANGE HERE ADD SMD, NMD FOR TAMREIN

  if(prediction) {

    # read classes
    classes_norut <- read.csv("data/classes_norut_prodchange.csv", header = FALSE)
    classes_norut_smd <- read.csv("data/classes_norut_smd_oneimpact.csv", header = FALSE)

    if(land_cover == "norut") {

      #-- NORUT
      # replace by names, as a factor
      dat[[paste0(prefix, "landcover_NORUT_2009")]] <- ifelse(dat[[paste0(prefix, "landcover_NORUT_2009")]] == 0, NA, dat[[paste0(prefix, "landcover_NORUT_2009")]])
      dat[[paste0(prefix, "landcover_NORUT_2009")]] <- classes_norut$V3[factor(dat[[paste0(prefix, "landcover_NORUT_2009")]])] |>
        factor(levels = unique(classes_norut$V3))

      dat[[paste0(prefix, "lc")]] <- dat[[paste0(prefix, "landcover_NORUT_2009")]]
      dat <- dat |>
        fastDummies::dummy_cols(select_columns = paste0(prefix, "lc"), ignore_na = TRUE, remove_selected_columns = TRUE)
    } else {

      if(land_cover == "norut_smd") {

        #-- NORUT-SMD
        # replace by names, as a factor
        dat[[paste0(prefix, "landcover_norut_smd")]] <- ifelse(dat[[paste0(prefix, "landcover_norut_smd")]] == 0, NA,
                                                               dat[[paste0(prefix, "landcover_norut_smd")]])
        dat[[paste0(prefix, "landcover_norut_smd")]] <- classes_norut_smd$V3[factor(dat[[paste0(prefix, "landcover_norut_smd")]])] |>
          factor(levels = unique(classes_norut_smd$V3))

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

        dat[[paste0(prefix, "lc")]] <- dat[[paste0(prefix, "landcover_norut_smd")]]
        dat <- dat |>
          fastDummies::dummy_cols(select_columns = paste0(prefix, "lc"), ignore_na = TRUE, remove_selected_columns = TRUE)

      }

    }

  }

  # land cover
  if(land_cover == "norut_smd") {

    if(land_cover_factor) {

      string <- "lc_nose_norut_smd"
      cols_lc_n <- grep(string, names(dat))
      names(dat)[cols_lc_n]
      names(dat)[cols_lc_n] <- sub(string, "lc", names(dat)[cols_lc_n])

      dat[[paste0(prefix, "landcover_norut_smd")]] <- relevel(dat[[paste0(prefix, "landcover_norut_smd")]], ref = ref_landcover)
      # remove lichen to the reference category
      # another possibility would be to move it to "other"
      #------- NEED TO BE CHANGED in case the ref category is not heathland
      dat[[paste0(prefix, "landcover_norut_smd")]] <- forcats::fct_collapse(dat[[paste0(prefix, "landcover_norut_smd")]], heathland = c("heathland", "lichen"))

    } else {

      string <- "landcover_norut_smd"
      cols_lc_n <- grep(string, names(dat))
      names(dat)[cols_lc_n]
      names(dat)[cols_lc_n] <- sub(string, ifelse(prediction, "lc_", "lc"),
                                   names(dat)[cols_lc_n])

      string <- "headland"
      cols_lc_n <- grep(string, names(dat))
      names(dat)[cols_lc_n]
      names(dat)[cols_lc_n] <- sub(string, "heathland", names(dat)[cols_lc_n])

    }

  } else {

    if(land_cover == "norut") {

      if(land_cover_factor) {

        string <- "lc_no_norut"
        cols_lc_n <- grep(string, names(dat))
        names(dat)[cols_lc_n]
        names(dat)[cols_lc_n] <- sub(string, "lc", names(dat)[cols_lc_n])

        dat[[paste0(prefix, "landcover_NORUT_2009")]] <- relevel(dat[[paste0(prefix, "landcover_NORUT_2009")]], ref = ref_landcover)
        # remove lichen to the reference category
        # another possibility would be to move it to "other"
        #------- NEED TO BE CHANGED in case the ref category is not heathland
        dat[[paste0(prefix, "landcover_NORUT_2009")]] <- forcats::fct_collapse(dat[[paste0(prefix, "landcover_NORUT_2009")]], heathland = c("heathland", "lichen"))
        # dat[[paste0(prefix, "landcover_NORUT_2009")]] |> levels()

      } else {

        string <- "landcover_NORUT_2009"
        cols_lc_n <- grep(string, names(dat))
        names(dat)[cols_lc_n]
        names(dat)[cols_lc_n] <- sub(string, "lc", names(dat)[cols_lc_n])

        string <- "headland"
        cols_lc_n <- grep(string, names(dat))
        names(dat)[cols_lc_n]
        names(dat)[cols_lc_n] <- sub(string, "heathland", names(dat)[cols_lc_n])

      }

    }

    #--- add options for SMD and NMD later for Sweden
  }

  #--------------------------
  # add option to use LC as factors instead (is there any difference?)
  # using the tables for reclassifying from numbers to strings

  #--
  # climate phenology

  # onset of spring
  names(dat)[grep(paste0(prefix, "sprAverageAllYrs"), names(dat))] <- paste0(prefix, "onset_spring")
  # bio12 annual prec
  names(dat)[grep(paste0(prefix, "BIOCLIM_bio12_historic"), names(dat))] <- paste0(prefix, "bio12_annual_prec")
  # bio18 warmest prec
  names(dat)[grep(paste0(prefix, "BIOCLIM_bio18_historic"), names(dat))] <- paste0(prefix, "bio18_warmest_prec")
  # bio19 coldest prec
  names(dat)[grep(paste0(prefix, "BIOCLIM_bio19_historic"), names(dat))] <- paste0(prefix, "bio19_coldest_prec")
  # ENVIREM continentality
  names(dat)[grep(paste0(prefix, "ENVIREM_continentality_historic"), names(dat))] <- paste0(prefix, "continentality")
  # ENVIREM growing degree days 5C
  names(dat)[grep(paste0(prefix, "ENVIREM_growingDegDays5_historic"), names(dat))] <- paste0(prefix, "growing_deg_days5")
  # snow coverdays, CHELSA snow cover days
  names(dat)[grep(paste0("snow_cover_duration"), names(dat))] <- sub("snow_cover_duration", "snow_cover_days", names(dat)[grep(paste0("snow_cover_duration"), names(dat))])
  names(dat)[grep(paste0(prefix, "CHELSA_scd"), names(dat))] <- paste0(prefix, "chelsa_scd")

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
  for(i in cols_graz_avg_n) {
    dat[,i][is.na(dat[,i])] <- 0
  }

  if(prediction) {

    grep("river", names(dat), value = T)
    grep("reservoir", names(dat), value = T)
    # grep("lake", names(sddat), value = T)

    # adding crossing variables
    # river
    river_var_formula <- grep("river", all_vars, value = T)
    if(length(river_var_formula) != 1) river_var_formula <- "cross_rivers_large_bin_sum"
    dat[[river_var_formula]] <- dat[["rivers_large"]]
    # lakes
    lakes_var_formula <- grep("river", all_vars, value = T)
    if(length(lakes_var_formula) != 1) lakes_var_formula <- "cross_lakes_bin_sum"
    dat[[lakes_var_formula]] <- dat[["lakes"]]
    # lakes
    reserv_var_formula <- grep("reservoir", all_vars, value = T) |>
      grep("cross", all_vars, value = T)
    if(length(reserv_var_formula) != 1) reserv_var_formula <- "cross_reservoirs_hydro_bin_sum"
    dat[[reserv_var_formula]] <- dat[["reservoirs_hydro"]]

    # fill NAs
    dat[[paste0(prefix, "onset_spring")]][is.na(dat[[paste0(prefix, "onset_spring")]])] <- max(dat[[paste0(prefix, "onset_spring")]], na.rm = TRUE)

    if(sum(grepl(paste0(prefix, "lichen_nina"), names(dat))) > 0)
      dat[[paste0(prefix, "lichen_nina")]][is.na(dat[[paste0(prefix, "lichen_nina")]])] <- 0#mean(dat[[paste0(prefix, "lichen_nina")]], na.rm = TRUE)
    if(sum(grepl(paste0(prefix, "lichen_slu"), names(dat))) > 0)
      dat[[paste0(prefix, "lichen_slu")]][is.na(dat[[paste0(prefix, "lichen_slu")]])] <- 0#mean(dat[[paste0(prefix, "lichen_slu")]], na.rm = TRUE)
    # dat[[paste0(prefix, "lichen_nose")]][is.na(dat[[paste0(prefix, "lichen_nose")]])] <- mean(dat[[paste0(prefix, "lichen_nose")]], na.rm = TRUE)

    # compute length of crossing all water
    dat[["cross_water_length"]] <- as.numeric(dat[[river_var_formula]]) + as.numeric(dat[[lakes_var_formula]]) + as.numeric(dat[[reserv_var_formula]])
  }



  if(prediction) {

    zoi_shape <- "bartlett"

    # vars_table <- tibble::tribble(
    #   ~name_formula, ~name_envdata,
    #   "houses", "houses",
    #   "roads_major", "roads_major_bartlett",
    #   "roads_minor", "roads_minor_bartlett",
    #   "railways", "railways_bartlett",
    #   "powerlines", "powerlines_bartlett",
    #   "cabins_private", "cabins_private_bartlett",
    #
    # )

    vars_table <- tibble::tibble(
      name_formula = c("houses", "roads_major", "roads_minor", "railways", "powerlines", "wind_turbines",
                       "cabins_private", "cabins_public_high", "cabins_public_low", "trails_log_pseudotui",
                       "grazing_animals_2018", "skitracks_high", "skitracks_low"),
      name_envdata = paste0(name_formula, "_", zoi_shape),
      suffix = ifelse(grepl("cabins_public|skitracks", name_formula), "", "_")) |>
      dplyr::bind_rows(
        tibble::tibble(
          name_formula = c("tpi", "onset", "snow_cover", "slope", "solar_radiation", "lakes", "reservoirs", "lc", "lichen_nina"),
          name_envdata = paste0(name_formula),
          suffix = ifelse(grepl("onset|snow_cover", name_formula), "", "_")
        )
      )

    ii = 1
    for(ii in seq_len(nrow(vars_table))) {

      string <- vars_table$name_formula[ii]
      cols_rhou_n <- grep(string, all_vars)
      (rmin_vars <- all_vars[cols_rhou_n])
      string_env <- vars_table$name_envdata[ii]
      (rmin_vars_envdat <- grep(string_env, names(dat), value = TRUE))

      i <- rmin_vars_envdat[2]
      for(i in rmin_vars_envdat) {
        if(length(nn <- grep(paste0(i, vars_table$suffix[ii]), rmin_vars)) > 0 & i != "lc_lichen") {
          # print(i)
          if(length(nn) > 1) nn <- nn[1]
          dat[[rmin_vars[nn]]] <- dat[[i]]
        }
      }

    }


    ######
    # add these variables!!
    # "along_tpi_250_mean"                "along_tpi_2500_mean"               "along_onset_spring"                "along_snow_cover_days"             "dem_diff_range_startend"           "along_slope_max"
    # [7] "along_solar_radiation_mean"        "along_lichen_nina_mean"            "along_lakes_sum"                   "along_reservoirs_hydro_sum"        "along_lc_conif_forest_mean"        "along_lc_mix_forest_mean"
    # [13] "along_lc_lichen_conif_forest_mean" "along_lc_decid_forest_mean"        "along_lc_lichen_decid_forest_mean" "along_lc_mires_mean"               "along_lc_wet_mires_mean"           "along_lc_exp_ridges_mean"
    # [19] "along_lc_grasses_mean"             "along_lc_meadows_mean"             "along_lc_snowbed_mean"             "along_lc_snow_mean"                "along_lc_agriculture_mean"         "along_lc_builtup_mean"
    # [25] "along_lc_clear_cut_mean"           "along_lc_other_mean"

  }

  if(prediction) {
    dat <- dat |>
      dplyr::mutate(
        # along_tpi_500_mean_small = along_tpi_500_mean,
        # along_tpi_500_mean_small_2 = along_tpi_500_mean**2,
        # along_tpi_1000_mean_2 = along_tpi_1000_mean**2,
        # along_tpi_5000_mean_large = along_tpi_1000_mean,
        # along_tpi_1000_mean_large_2 = along_tpi_5000_mean**2,
        # along_tpi_250_mean_small = along_tpi_250_mean,
        along_tpi_250_mean_small_2 = along_tpi_250_mean_small**2,
        # along_tpi_250_mean_small = along_tpi_2500_mean,
        along_tpi_2500_mean_large_2 = along_tpi_2500_mean_large**2,
        # onset_spring_2 = onset_spring**2,
        along_snow_cover_days_mean_2 = along_snow_cover_days_mean**2,
        # along_slope_mean_2 = along_slope_mean**2,
        along_slope_max_2 = along_slope_max**2,
        # solar_radiation_2 = solar_radiation**2,
        # bio12_annual_prec_2 = bio12_annual_prec**2,
        # bio18_warmest_prec_2 = bio18_warmest_prec**2,
        # bio19_coldest_prec_2 = bio19_coldest_prec**2,
        # continentality_2 = continentality**2,
        # growing_deg_days5_2 = growing_deg_days5**2,
        # chelsa_scd_2 = chelsa_scd**2
      )
  } else {
    dat <- dat |>
      dplyr::mutate(
        along_tpi_500_mean_small = along_tpi_500_mean,
        along_tpi_500_mean_small_2 = along_tpi_500_mean_small**2,
        along_tpi_1000_mean_2 = along_tpi_1000_mean**2,
        along_tpi_5000_mean_large = along_tpi_5000_mean,
        along_tpi_5000_mean_large_2 = along_tpi_5000_mean_large**2,
        along_tpi_250_mean_small = along_tpi_250_mean,
        along_tpi_250_mean_small_2 = along_tpi_250_mean_small**2,
        along_tpi_2500_mean_large = along_tpi_2500_mean,
        along_tpi_2500_mean_large_2 = along_tpi_2500_mean_large**2,
        onset_spring_2 = onset_spring**2,
        along_snow_cover_days_mean_2 = along_snow_cover_days_mean**2,
        along_slope_mean_2 = along_slope_mean**2,
        along_slope_max_2 = along_slope_max**2,
        # solar_radiation_2 = solar_radiation**2,
        # bio12_annual_prec_2 = bio12_annual_prec**2,
        # bio18_warmest_prec_2 = bio18_warmest_prec**2,
        # bio19_coldest_prec_2 = bio19_coldest_prec**2,
        # continentality_2 = continentality**2,
        # growing_deg_days5_2 = growing_deg_days5**2,
        # chelsa_scd_2 = chelsa_scd**2
      )
  }




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

  # ADD OTHER OPERATIONS PRESENT IN THE DATA PREPARATION HERE AS WELL!?
  # FOR THE PREDICTIONS BASED ON THE GRID

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
