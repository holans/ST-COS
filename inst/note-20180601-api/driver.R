#' ---
#' title: How to Retrieve Data from the Census API
#' author: Andrew Raim
#' date: June 1, 2018
#' output: pdf_document
#' ---

#' Here is an example of using the `rjson` package to download ACS estimates.
library(rjson)

# How to call Census JSON API
# https://census.gov/content/dam/Census/programs-surveys/acs/guidance/training-presentations/06212017_ACS_Census_API.pdf
json_file <- "https://api.census.gov/data/2013/acs1?get=NAME,B17001_002E&for=state:*"
json_data <- fromJSON(file = json_file)

json_file <- "https://api.census.gov/data/2015/acs5?get=NAME,B17001_002E&for=county:*"
json_data <- fromJSON(file = json_file)

#' Here is an example of using the `jsonlite` package to download ACS estimates.
#' It returns the data in a more convenient form. Note that ACS data available
#' through the API only seems to go back to 2012.
library(jsonlite)

json_file <- "https://api.census.gov/data/2015/acs5?get=NAME,B17001_002E&for=county:*"
json_data <- fromJSON(json_file)
dat <- data.frame(json_data[-1,])
colnames(dat) <- json_data[1,]

json_file <- "https://api.census.gov/data/2012/acs1?get=NAME,B17001_002E&for=county:*"
json_data <- fromJSON(json_file)
dat <- data.frame(json_data[-1,])
colnames(dat) <- json_data[1,]

#' # Now what about getting shapefiles??

#' Tigris is a convenient way to grab shapefiles in R
library(rgeos)
library(tigris)
library(sf)

# dfw <- tracts(state = 'TX', county = c('Dallas', 'Tarrant'))
dfw <- tigris::counties(state = 'PA', year = 2012)
df <- st_as_sf(dfw)
plot(df[,3])

dfw <- tigris::counties(year = 2011, cb = TRUE)
df <- st_as_sf(dfw)
plot(df[,3])
