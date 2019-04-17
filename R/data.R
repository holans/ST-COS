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
#' Shapefiles were gathered via the \code{tigris} package, and ACS estimates were
#' downloaded from the American FactFinder \url{http://factfinder.census.gov}.
#' Data was assembled on 2/28/2019. See \code{data-prep-aff.R} in the Columbia
#' example code for details.
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

