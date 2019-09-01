library(readr)
library(dplyr)

# Assemble Source Supports
# We downloaded ACS data files from the Census Bureau's American FactFinder website
# <https://factfinder.census.gov> for:
# - Boone County, Missouri
# - Table B19013
# - Block group level geography
# - Years 2013 - 2017
#
# The following loop reads those files and puts them in a form needed for our
# analysis. We merge data with shapefiles, which are requested from the Census Bureau's
# TIGER/Line service via the `tigris` package.

year_levels = 2013:2017
sf_list = list()

for (idx in 1:length(year_levels)) {
	year = year_levels[idx]
	acs_file = sprintf("ACS_%s_5YR_B19013_with_ann.csv", substr(year,3,4))

	my_dat = read_csv(acs_file) %>%
		slice(-1) %>%
		mutate(geoid = GEO.id2) %>%
		mutate(HD01_VD01 = replace(HD01_VD01, HD01_VD01 == '-', NA)) %>%
		mutate(HD02_VD01 = replace(HD02_VD01, HD02_VD01 == '**', NA)) %>%
		mutate(DirectEst = as.numeric(HD01_VD01)) %>%
		mutate(DirectMOE = as.numeric(HD02_VD01)) %>%
		mutate(DirectVar = (DirectMOE^2 / qnorm(0.95)^2)) %>%
		mutate(state = substr(geoid, 1, 2)) %>%
		mutate(county = substr(geoid, 3, 5)) %>%
		mutate(tract = substr(geoid, 6, 11)) %>%
		mutate(blockgroup = substr(geoid, 12, 12)) %>%
		select(geoid, state, county, tract, blockgroup, DirectEst, DirectMOE, DirectVar) %>%
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

# Our assembled source supports are now `acs5_2013`, ..., `acs5_2017`
acs5_2013 = sf_list[[1]]
acs5_2014 = sf_list[[2]]
acs5_2015 = sf_list[[3]]
acs5_2016 = sf_list[[4]]
acs5_2017 = sf_list[[5]]
save(acs5_2013, acs5_2014, acs5_2015, acs5_2016, acs5_2017, file="data/acs_sf.rda")

