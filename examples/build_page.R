# Build page
library(pkgdown)
library(usethis)

# Run once to configure package to use pkgdown
usethis::use_pkgdown()
# Run to build the website
pkgdown::build_site()

# CI
usethis::use_pkgdown_github_pages()
