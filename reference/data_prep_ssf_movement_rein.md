# Prepare data for SSF movement/permeability analyses for wild and semi-domesticated reindeer within SAM

This function takes in the already annotated data for SSF with either
wild or semi-domesticated reindeer and further prepare columns and
variables for the model fitting procedure. The function is suited for
step resource selection functions (SSF) only.

## Usage

``` r
data_prep_ssf_movement_rein(
  dat,
  season,
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
  formula = NULL
)
```

## Arguments

- dat:

  `[data.frame]`  
  Annotated data set for analysis.

- season:

  `[string]{"sum", "cal", "win"}`  
  Season of interest. One of "sum", "cal", or "win".

- prediction:

  `[logical(1)=FALSE]`  
  Additional changes that should only be done in the grid for
  prediction.

- land_cover:

  `[string(1)="norut_smd"]{"norut_smd", "norut", "smd", "nmd"}`  
  Which land cover map should be used for analysis. It prepares the
  corresponding classes as

- include_zoi_nearest:

  `[logical(1)=FALSE]`  
  By default, `FALSE`. It should be set to `TRUE` if there are variables
  representing the ZOI of the nearest feature in the formula. Relevant
  only for prediction (using the grid).

- prefix:

  `[string(1)="endpt_]{"endpt_", "along_", ""}`  
  Prefix for the variable, whether `"endpt_"` or `"along_"` for
  variables measures in the end point or along the step, and `""` for
  normal RSF or for the prediction using the grid.
