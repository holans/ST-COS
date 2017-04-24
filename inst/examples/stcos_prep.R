library(R6)
library(stcos)
library(sf)

setwd("/home/araim/Documents/simulations/ST-COS/run-20170324-checkoverlaps-sf")

load.domain <- function(shpfile, datfile, layername, crs.tx = NULL, crs.orig = NULL)
{
	if (!is.null(crs.orig)) {
		area <- st_read(shpfile, layer = layername, crs = crs.orig)
	} else {
		area <- st_read(shpfile, layer = layername)
	}
	if (!is.null(crs.tx)) {
		area <- st_transform(area, crs = crs.tx)
	}

	states <- area$STATE
	# plot(area)

	# Extract the continental U.S.
	idx.HI <- which(states == '02')
	idx.AK <- which(states == '15')
	idx.PR <- which(states == '72')
	idx.keep <- setdiff(1:nrow(area), c(idx.HI, idx.AK, idx.PR))

	area <- area[idx.keep,]
	# plot(area)

	dat <- read.csv(datfile, header = FALSE)
	eZagg <- dat[idx.keep, 1]
	esigmavar <- (dat[idx.keep, 2] / 1.645)^2		# known survey variances
	Zagg <- log(eZagg)
	sigmavar <- esigmavar * (1 / eZagg)^2

	area$EST <- eZagg
	area$VAR <- sigmavar
	return(area)
}

res3 <- load.domain("shp/period3.shp", "shp/period3.csv", "period3")
sp <- STCOSPrep$new(res3)

# ACS 1-year estimates for 2013
res1 <- load.domain("shp/period1.shp", "shp/period1.csv", "period1", crs.tx = st_crs(res3$area))
sp$add_obs(res1, time = 2013, period = 2013, estimate_name = "EST", variance_name = "VAR")

# ACS 5-year estimates for 2012
res3.2012 <- load.domain("shp/period3_2012.shp", "shp/period3_2012.csv", "period3_2012", crs.tx = st_crs(res3$area))
sp$add_obs(res3.2012, time = 2012, period = 2009:2013, estimate_name = "EST", variance_name = "VAR")


