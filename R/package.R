#' stcos: Space-Time Change of Support
#'
#' An R Package for Space-Time Change of Support (STCOS) modeling.
#' 
#' @useDynLib stcos, .registration = TRUE
#' @import sf
#' @import Rcpp
#' @import Matrix
#' @importFrom stats dgamma dnorm logLik optim rgamma rnorm runif quantile sd
#' @importFrom R6 R6Class
#' @importFrom Rcpp evalCpp
#' @importFrom dplyr mutate row_number
#' 
#' @details Supports building and running STCOS and related models. A guide
#' on package use is given in <arXiv:1904.12092>.
#' 
#' @references
#' Jonathan R. Bradley, Christopher K. Wikle, and Scott H. Holan (2015).
#' Spatio-temporal change of support with application to American Community
#' Survey multi-year period estimates. STAT 4 pp.255-270.
#' \url{https://doi.org/10.1002/sta4.94}.
#'
#' Andrew M. Raim, Scott H. Holan, Jonathan R. Bradley, and Christopher K.
#' Wikle (2017). A model selection study for spatio-temporal change of
#' support. In JSM Proceedings, Government Statistics Section. Alexandria,
#' VA: American Statistical Association, pp.1524-1540.
#' 
#' Andrew M. Raim, Scott H. Holan, Jonathan R. Bradley, and Christopher K.
#' Wikle (2019+). An R package for spatio-temporal change of support.
#' \url{https://arxiv.org/abs/1904.12092}.
#' 
#' @docType package
#' @name stcos
NULL
