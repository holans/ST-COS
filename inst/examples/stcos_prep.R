library(stcos)
library(sf)

setwd("/home/araim/Documents/simulations/ST-COS")

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
	area$VAR <- esigmavar
	return(area)
}

# Fine-level domain comes from ACS 5-year estimates for 2013
acs5.2013 <- load.domain("shp/period3.shp", "shp/period3.csv", "period3")

# Set up knots for bisquare basis
knots.sp <- read.csv("dat/knots250_ACS_amr.csv", head = FALSE)
knots.t <- c(2012.5, 2011, 2011.5, 2011, 2010, 2010.5, 2010, 2009, 2009.5, 2009,
	2008, 2008.5, 2008, 2007, 2007.5, 2007, 2006.5, 2006, 2005.5)
knots <- merge(knots.sp, knots.t)
names(knots) <- c("x", "y", "t")
basis <- BisquareBasis$new(knots[,1], knots[,2], knots[,3], w.s = 0.5, w.t = 0.5)

# Construct a STCOSPrep object, then add space-time domains with observations
sp <- STCOSPrep$new(fine_domain = acs5.2013, basis = basis)

# ACS 1-year estimates for 2013
acs1.2013 <- load.domain("shp/period1.shp", "shp/period1.csv", "period1", crs.tx = st_crs(acs5.2013))
sp$add_obs(acs1.2013, time = 2013, period = 2013, estimate_name = "EST", variance_name = "VAR")

# ACS 5-year estimates for 2012
acs5.2012 <- load.domain("shp/period3_2012.shp", "shp/period3_2012.csv", "period3_2012", crs.tx = st_crs(acs5.2013))
sp$add_obs(acs5.2012, time = 2012, period = 2008:2012, estimate_name = "EST", variance_name = "VAR")


Z <- sp$get_Z()
V <- sp$get_V()
H <- sp$get_H()

# Here's where we would think about a reduction for S
S <- sp$get_S()
sp$set_basis_reduction(identity)
S.reduced <- sp$get_reduced_S()

