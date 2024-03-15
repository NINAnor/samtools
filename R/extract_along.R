#' Extract values of raster using in parallel
#'
#' Function to extract values of rasters along lines using parallel computation.
#' The function is actually not very optimized currently.
#'
#' @param r `[SpatRast,data.frame]` \cr Either a set of rasters, in the `SpatRaster`
#' format from [terra]; or a `data.frame` resulting of `terra::extract` for
#' linestrings, for a set of environmental variables (defined here as the columns).
#' @param sl Vector of step lines, with a column representing the step id
#' (stratum or just step number) and a LINESTRING geometry for each step.
#' Can be a `sf` or a `SpatVector` object.
#' @param ws List of which summary statistics are to be computed. Currently, only "mean", "max",
#' "min", and "sum" are implemented.
#' @param col_step_length `[character="dist"]` \cr String with the name of the
#' column representing step length.
#' @param step_id `[character="use_ava_data_animals_id"]` \cr String with the name
#' of the column representing step ID.
#' @param prefix `[character="along_"]` \cr Prefix to be added to the extracted
#' variable names.
#'
#' @name extract_along
#'
#' @export
extract_along <- function(r, sl, ws,
                          col_step_length = "dist",
                          step_id = "use_ava_data_animals_id",
                          prefix = "along_") {
  UseMethod("extract_along")
}


#' @rdname extract_along
#' @export
extract_along.SpatRaster <- function(r, sl, ws,
                                     col_step_length = "dist",
                                     step_id = "use_ava_data_animals_id",
                                     prefix = "along_"){

  # check number of summaries
  n_summaries <- sapply(ws, length)
  # repeat raster indices
  # not needed, it would require performing the extraction multiple times
  # rasts <- rep(r, n_summaries)
  # summaries
  # we add "mean" to account for the stepID
  ws <- c("mean", unlist(ws))

  # remove NAs
  whichNA <- which(is.na(sl[[col_step_length]]))
  whichkeep <- which(!is.na(sl[[col_step_length]]))

  # extract values
  extracted_values <- terra::extract(r, sl[whichkeep,])
  # extracted_values <- terra::extract(r, terra::vect(sl), fun=mean)
  # ext_v <- raster::extract(raster::raster(r), sl)

  # repeat column for which there is more than one summary statistic
  n_summaries <- c(1, n_summaries) # add 1 for ID
  indices <- rep(1:ncol(extracted_values), times = n_summaries)
  ext_val_rep <- extracted_values[,indices]
  # split extracted values by step (column ID automatically created)
  ext_val_split <- split(ext_val_rep, extracted_values$ID)

  # #max slope
  # step_values <- data.frame(step_id = sl$use_ava_data_animals_id)
  #
  # tmp <- unlist(lapply(extracted_values, function(x){ get_summary(matrix(x[,1]), "max")}))
  # step_values$max_slope <- tmp

  #all other variables
  summaries <- data.frame(do.call("rbind", lapply(ext_val_split, get_summary, y = ws)))

  colnames(summaries) <- c(step_id,
                           paste0(prefix,
                                  sapply(strsplit(colnames(ext_val_rep)[-1], split="@"), function(x) x[1]), "_", ws[-1]))

  # re-add NAs
  if(length(whichNA) > 0) {
    summariesNA <- summaries[1:length(whichNA),]
    summariesNA[[step_id]] <- sl[[step_id]][whichNA]
    summariesNA[,!grepl(step_id, colnames(summariesNA))] <- NA
    summaries <- rbind(summaries, summariesNA)
    summaries <- summaries[order(c(whichkeep, whichNA)),]
  }

  summaries[,1] <- sl[[step_id]]

  return(summaries)
}

#' @rdname extract_along
#' @export
extract_along.data.frame <- function(r, sl, ws,
                                     col_step_length = "dist",
                                     step_id = "use_ava_data_animals_id",
                                     prefix = "along_"){

  # check number of summaries
  n_summaries <- sapply(ws, length)
  # repeat raster indices
  # not needed, it would require performing the extraction multiple times
  # rasts <- rep(r, n_summaries)
  # summaries
  # we add "mean" to account for the stepID
  ws <- c("mean", unlist(ws))

  # remove NAs
  # whichNA <- which(is.na(sl[[col_step_length]]))
  # whichkeep <- which(!is.na(sl[[col_step_length]]))

  # extract values
  extracted_values <- r
  # extracted_values <- terra::extract(r, terra::vect(sl), fun=mean)
  # ext_v <- raster::extract(raster::raster(r), sl)

  # repeat column for which there is more than one summary statistic
  n_summaries <- c(1, n_summaries) # add 1 for stepID
  indices <- rep(1:ncol(extracted_values), times = n_summaries)
  ext_val_rep <- extracted_values[,indices]
  # split extracted values by step (column ID automatically created)
  step_id_factor <- factor(as.character(extracted_values[[step_id]]),
                           levels = unique(as.character(extracted_values[[step_id]])))
  ext_val_split <- split(ext_val_rep, step_id_factor)

  # #max slope
  # step_values <- data.frame(step_id = sl$use_ava_data_animals_id)
  #
  # tmp <- unlist(lapply(extracted_values, function(x){ get_summary(matrix(x[,1]), "max")}))
  # step_values$max_slope <- tmp

  #all other variables
  summaries <- data.frame(do.call("rbind", lapply(ext_val_split, get_summary, y = ws)))

  colnames(summaries) <- c(step_id,
                           paste0(prefix,
                                  sapply(strsplit(colnames(ext_val_rep)[-1], split="@"), function(x) x[1]), "_", ws[-1]))

  # re-add NAs
  if(length(whichNA) > 0) {
    summaries[whichNA, !grepl(step_id, colnames(summariesNA))] <- NA
  }

  # summaries[,1] == sl[[step_id]]

  return(summaries)
}

# function to get summary statistics for a list of
# matrices with values extracted from a raster
get_summary <- function(x, y){

  toto <- function(a,b,c){
    if (c[a]=="mean"){ res <- base::mean(as.numeric(as.character(b[,a])), na.rm=T) }
    if (c[a]=="max") { res <- base::max(as.numeric(as.character(b[,a])), na.rm=T) }
    if (c[a]=="min") { res <- base::min(as.numeric(as.character(b[,a])), na.rm=T) }
    if (c[a]=="sum") { res <- base::sum(as.numeric(as.character(b[,a])), na.rm=T) }
    return(res)
  }
  res <- unlist(lapply(c(1:ncol(x)), toto, b = x, c = y))
  # for(i in c(1:ncol(x))) toto(a = i, b = x, c = y)
  return(res)
}
