library(sf)
library(dplyr)

options(width = 80)

neighbs = st_read("neighborhoods.shp") %>% st_transform(crs = 3857)

year = 2017
est_url = paste('https://api.census.gov/data/', year,
	'/acs/acs5?get=NAME,B19013_001E&for=block%20group:*&in=state:29+county:019',
	sep = '')
moe_url = paste('https://api.census.gov/data/', year,
	'/acs/acs5?get=NAME,B19013_001M&for=block%20group:*&in=state:29+county:019',
	sep = '')

json_data = jsonlite::fromJSON(est_url)
est_dat = data.frame(json_data[-1,])
colnames(est_dat) = json_data[1,]

json_data = jsonlite::fromJSON(moe_url )
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
head(my_dat)

my_shp = tigris::block_groups(state = '29', county = '019', year = 2017) %>%
	st_as_sf() %>%
	st_transform(crs = 3857)

acs5_2017 = my_shp %>%
	inner_join(my_dat, by = c('STATEFP' = 'state', 'COUNTYFP' = 'county',
		'TRACTCE' = 'tract', 'BLKGRPCE' = 'blockgroup')) %>%
	select(geoid = GEOID, state = STATEFP, county = COUNTYFP,
		tract = TRACTCE, blockgroup = BLKGRPCE,
		DirectEst, DirectMOE, DirectVar)

head(acs5_2017)

dom_fine = tigris::block_groups(state = '29', county = '019', year = 2017) %>% 
	st_as_sf() %>%
	st_transform(crs = 3857)
