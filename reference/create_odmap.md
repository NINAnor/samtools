# Create folder with ODMAP structure

This function creates an empty folder following the ODMAP structure
(Overview, Data, Modeling, Assessment, and Prediction). It is intended
to be used for habitat, movement, and connectivity modeling and
suggested to be used by species (i.e. one folder per species), but this
can be used in other ways depending of the user preferences.

## Usage

``` r
create_odmap(path = ".", overwrite = FALSE)
```

## Arguments

- path:

  `[character]`  
  Path where the folder will be created, including intended folder name.
  If no argument is provided, the structure is created in the current
  working directory,

- overwrite:

  `[logical(1)=FALSE]`  
  Should the folder be overwritten if' it already exists? Default is
  `FALSE`.

## Details

The function will probably throw some warnings because R project files
are created in nested folders where a R proj file is already existing.
The user needs to mannually accept that to allow the creation of these
files, which makes it easier for one to set relative paths for each of
the modeling steps.

## Examples

``` r
create_odmap(path = "./test")
#> ✔ Creating ./test/.
#> ✔ Creating ./test/01_overview/.
#> ! New project ./test/02_data would be nested inside an existing project
#>   /home/runner/work/samtools/samtools/, which is rarely a good idea.
#> ℹ If this is unexpected, the here package has a function, `here::dr_here()`
#>   that reveals why a particular path is regarded as a project. To learn more,
#>   run `here::dr_here()` in a fresh R session that has
#>   /home/runner/work/samtools/samtools/ as working directory.
#> Error in ui_yep(x = x, yes = yes, no = no, n_yes = n_yes, n_no = n_no,     shuffle = shuffle, .envir = .envir): ✖ User input required, but session is not interactive.
#> ℹ Query: "Do you want to create anyway?"
```
