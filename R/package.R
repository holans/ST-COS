# We needed some workarounds to get the package to work. See:
# - https://stackoverflow.com/questions/36605955/c-function-not-available
# - https://cran.r-project.org/web/packages/roxygen2/vignettes/namespace.html

#' @useDynLib stcos, .registration = TRUE
#' @import sf
#' @import Rcpp
#' @import Matrix
#' @importFrom MASS ginv
#' @importFrom stats dgamma dnorm logLik optim rgamma rnorm runif quantile sd
#' @importFrom R6 R6Class
#' @importFrom Rcpp evalCpp
NULL

