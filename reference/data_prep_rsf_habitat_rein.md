# Prepare data for RSF analyses for wild and semi-domesticated reindeer within SAM

This function takes in the already annotated data for RSF with either
wild or semi-domesticated reindeer and further prepare columns and
variables for the model fitting procedure. The function is suited for
point resource selection functions (RSF) only.

## Usage

``` r
data_prep_rsf_habitat_rein(
  dat,
  season,
  prediction = FALSE,
  fixwind = TRUE,
  land_cover_factor = FALSE,
  land_cover = c("norut_smd", "norut", "smd", "nmd")[1],
  ref_landcover = "heathland",
  include_zoi_nearest = FALSE,
  species = c("wrein", "trein")[1],
  prefix = "",
  reference_year = lubridate::year(lubridate::now()),
  use_binary_vars_as_categorical = FALSE
)
```

## Arguments

- dat:

  `[data.frame]`  
  Data set annotated with environmental covariates, almost ready for
  analysis.

- season:

  `[string]{"sum", "cal", "win"}`  
  Season of interest. One of "sum", "cal", or "win".

- prediction:

  `[logical(1)=FALSE]`  
  Additional changes that should only be done in the grid for
  prediction.

- land_cover_factor:

  `[logical(1)=FALSE]`  
  Logical variable stating whether the land cover variable should be
  treated as a factor (one single column with land cover type names) or
  not (as multiple columns as dummy variables).

- land_cover:

  `[string(1)="norut_smd"]{"norut_smd", "norut", "smd", "nmd"}`  
  Which land cover map should be used for analysis. It prepares the
  corresponding classes, using the `ref_landcover` as the reference
  class.

- ref_landcover:

  `[string="heathland"]`  
  Reference class for land cover variables.

- include_zoi_nearest:

  `[logical(1)=FALSE]`  
  If `TRUE`, variables representing the zone of influence (ZOI) of the
  nearest feature are also computed, based on the distance to the
  nearest features. Only relevant if `prediction = TRUE`.

- species:

  `[string="wrein"]{"wrein", "trein"}`  
  The species/management system for which we should prepare the data
  for. Either `"wrein"` for wild reindeer or `"trein"` for
  semi-domesticated reindeer. It only implies some differences in data
  which are particular to each of the species/management systems.

- prefix:

  `[string=""]`  
  Prefix to be added to the variable/column names of the dataset.
  Default is `""`, but it could be e.g. `"startpt_"`, `"endpt_"`, or
  `"along_"` for variables extracted at the starting or ending point or
  along steps.

- use_binary_vars_as_categorical:

  `[logical(1)=FALSE]`  
  If `TRUE`, binary variables (e.g. lakes, reservoirs) should be treated
  as categorical. Otherwise, they are treated as numerical.
