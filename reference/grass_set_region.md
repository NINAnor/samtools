# Set the region in a GRASS GIS session

Helper function to set the computation region of an open GRASS GIS
session connected to the R session, e.g. through
[`rgrass::initGRASS()`](https://osgeo.github.io/rgrass/reference/initGRASS.html).
To know more about the arguments, check:
https://grass.osgeo.org/grass82/manuals/g.region.html

## Usage

``` r
grass_set_region(flags = c("p"), ...)
```

## Examples

``` r
if(FALSE) {
  grass_set_region(res = 30)
}
```
