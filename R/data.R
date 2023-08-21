#' City of Columbia neighborhoods.
#'
#' An \code{sf} object containing the geometry of four neighborhoods in the
#' City of Columbia, Boone County, Missouri. Based on shapefiles provided by
#' the Office of Information Technology / GIS, City of Columbia, Missouri.
#'
#' @format An \code{sf} object with 4 features (neighborhoods).
"columbia_neighbs"

#' Shapes and ACS estimates for Boone County, MO.
#'
#' An \code{sf} object with ACS estimates for:
#' \itemize{
#' \item Boone County, Missouri
#' \item Table B19013
#' \item Block group level geography
#' \item Years 2013 - 2017
#' }
#' 
#' This dataset was assembled using shapefiles from the \code{tigris} package
#' and ACS estimates from the American FactFinder website on 2/28/2019.
#' See \code{data-prep-aff.R} in the Columbia example code for details.
#' American FactFinder has since been deprecated, and similar data are
#' available at \url{http://data.census.gov}.
#'
#' @format \code{sf} objects.
#' @name acs_sf
NULL

#' @name acs_sf
"acs5_2013"

#' @name acs_sf
"acs5_2014"

#' @name acs_sf
"acs5_2015"

#' @name acs_sf
"acs5_2016"

#' @name acs_sf
"acs5_2017"

