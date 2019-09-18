library(jsonlite)
library(sf)
library(tigris)
library(dplyr)
library(ggplot2)

# Assemble Source Supports
# The following loop pulls ACS direct estimates and associated MOEs from the
# Census Bureau's Data API, and merges them into a single data frame. For more
# information about the API, see
# <https://www.census.gov/data/developers/guidance/api-user-guide.html>.
# Examples of URLs to pull direct estimates and associated MOEs:
# 
# - <https://api.census.gov/data/2015/acs/acs5?get=NAME,B19013_001E&for=block%20group:*&in=state:29+county:019>
# - <https://api.census.gov/data/2015/acs/acs5?get=NAME,B19013_001M&for=block%20group:*&in=state:29+county:019>
# 
# Note: Block group data for years earlier than 2015 appear to require a
# slightly different call through the API - there we need to specify which
# tract(s) we are interested in, as well.

year_levels = 2015:2017
sf_list = list()
# dat.missing = list()

for (idx in 1:length(year_levels)) {
	year = year_levels[idx]

	est_url = paste('https://api.census.gov/data/', year,
		'/acs/acs5?get=NAME,B19013_001E&for=block%20group:*&in=state:29+county:019',
		sep = '')
	json_data = fromJSON(est_url)
	est_dat = data.frame(json_data[-1,])
	colnames(est_dat) = json_data[1,]

	moe_url = paste('https://api.census.gov/data/', year,
		'/acs/acs5?get=NAME,B19013_001M&for=block%20group:*&in=state:29+county:019',
		sep = '')
	json_data = fromJSON(moe_url)
	moe_dat = data.frame(json_data[-1,])
	colnames(moe_dat) = json_data[1,]

	my_dat = est_dat %>%
		inner_join(moe_dat, by = c('state' = 'state', 'county' = 'county',
			'tract' = 'tract', 'block group' = 'block group')) %>%
		select(state, county, tract, blockgroup = `block group`,
			DirectEst = B19013_001E, DirectMOE = B19013_001M) %>%
		mutate(state = as.character(state)) %>%
		mutate(county = as.character(county)) %>%
		mutate(tract = as.character(tract)) %>%
		mutate(blockgroup = as.character(blockgroup)) %>%
		mutate(DirectEst = as.numeric(as.character(DirectEst))) %>%
		mutate(DirectMOE = as.numeric(as.character(DirectMOE))) %>%
		mutate(DirectEst = replace(DirectEst, DirectEst < 0, NA)) %>%
		mutate(DirectMOE = replace(DirectMOE, DirectMOE < 0, NA)) %>%
		mutate(DirectVar = (DirectMOE / qnorm(0.95))^2) %>%
		arrange(tract, blockgroup)

	my_shp = block_groups(state = '29', county = '019', year = year) %>%
		st_as_sf() %>%
		st_transform(crs = 3857)

	my_sf = my_shp %>%
		inner_join(my_dat, by = c('STATEFP' = 'state', 'COUNTYFP' = 'county',
			'TRACTCE' = 'tract', 'BLKGRPCE' = 'blockgroup')) %>%
		select(geoid = GEOID, state = STATEFP, county = COUNTYFP,
			tract = TRACTCE, blockgroup = BLKGRPCE,
			DirectEst, DirectMOE, DirectVar)

	sf_list[[idx]] = my_sf
}

# Our assembled source supports are now `acs5_2015`, ..., `acs5_2017`
acs5_2015 = sf_list[[1]]
acs5_2016 = sf_list[[2]]
acs5_2017 = sf_list[[3]]
rm(sf_list)
