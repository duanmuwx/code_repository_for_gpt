
library(devtools)
library(roxygen2)
library(usethis)

usethis::use_package_doc()
#devtools::use_rcpp()
devtools::document()


build_manual(pkg = ".", path = NULL)
