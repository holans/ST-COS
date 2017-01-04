#library(rgdal)
#library(rgeos)
#library(raster)

# From Matlab code:
# proj4string(res3$area)
# mstruct = defaultm('eqaconicstd');
# mstruct.geoid = almanac('earth','grs80','foot');
# mstruct.origin = [24.395833 -96.000000];
# mstruct.falseeasting = 0.0;
# mstruct.falsenorthing = 0.0;
# mstruct.mapparallels = [29.5 45.5];

# TBD: Looks like we need to project all datasets onto the same geography...
# I tried to do something below, but I don't know if there's a 'proper' way to do it.
# shpfile <- "shp/period2_2007.shp"
# datfile <- "shp/period2_2007.csv"
# layername <- "period2_2007"
# proj4string(area.orig)

# area.orig <- readOGR(shpfile, layer = layername, verbose = FALSE)

get.period.data <- function(shpfile, datfile, layername)
{
	area.orig <- readOGR(shpfile, layer = layername, verbose = FALSE)
	#if (is.na(proj4string(area.orig))) {
	#	area.orig <- readOGR(shpfile, layer = layername, verbose = FALSE,
	#		p4s = "+proj=aea +lat_0=24.395833 +lon_0=-96 +x_0=0 +y_0=0 +units=us-ft +ellps=GRS80")
	#} else {
	#	area.orig <- spTransform(area.orig, "+proj=aea +lat_0=24.395833 +lon_0=-96 +x_0=0 +y_0=0 +units=us-ft +ellps=GRS80")
	#}
	# area.tx <- spTransform(area.orig, "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=24.395833 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=us-ft +no_defs +ellps=GRS80 +towgs84=0,0,0")

	states <- area.orig@data$STATE
	# plot(area)

	# Extract the continental U.S.
	idx.HI <- which(states == '02')
	idx.AK <- which(states == '15')
	idx.PR <- which(states == '72')
	idx.keep <- setdiff(1:nrow(area.orig), c(idx.HI, idx.AK, idx.PR))

	area <- area.orig[idx.keep,]
	# plot(area)

	area.dat <- read.csv(datfile, header = FALSE)
	eZagg <- area.dat[idx.keep, 1]
	esigmavar <- (area.dat[idx.keep, 2] / 1.645)^2		# known survey variances
	Zagg <- log(eZagg)
	sigmavar <- esigmavar * (1 / eZagg)^2
	
	list(area = area, area.dat = area.dat, eZagg = eZagg,
		esigmavar = esigmavar, Zagg = Zagg, sigmavar = sigmavar)
}

# period3 is our target geography
res1 <- get.period.data("shp/period1.shp", "shp/period1.csv", "period1")
res2 <- get.period.data("shp/period2.shp", "shp/period2.csv", "period2")
res3 <- get.period.data("shp/period3.shp", "shp/period3.csv", "period3")

res1.2006 <- get.period.data("shp/period1_2006.shp", "shp/period1_2006.csv", "period1_2006")
res1.2007 <- get.period.data("shp/period1_2007.shp", "shp/period1_2007.csv", "period1_2007")
res1.2008 <- get.period.data("shp/period1_2008.shp", "shp/period1_2008.csv", "period1_2008")
res1.2009 <- get.period.data("shp/period1_2009.shp", "shp/period1_2009.csv", "period1_2009")
res1.2010 <- get.period.data("shp/period1_2010.shp", "shp/period1_2010.csv", "period1_2010")
res1.2011 <- get.period.data("shp/period1_2011.shp", "shp/period1_2011.csv", "period1_2011")
res1.2012 <- get.period.data("shp/period1_2012.shp", "shp/period1_2012.csv", "period1_2012")

res2.2007 <- get.period.data("shp/period2_2007.shp", "shp/period2_2007.csv", "period2_2007")
res2.2008 <- get.period.data("shp/period2_2008.shp", "shp/period2_2008.csv", "period2_2008")
res2.2009 <- get.period.data("shp/period2_2009.shp", "shp/period2_2009.csv", "period2_2009")
res2.2010 <- get.period.data("shp/period2_2010.shp", "shp/period2_2010.csv", "period2_2010")
res2.2011 <- get.period.data("shp/period2_2011.shp", "shp/period2_2011.csv", "period2_2011")
res2.2012 <- get.period.data("shp/period2_2012.shp", "shp/period2_2012.csv", "period2_2012")

res3.2009 <- get.period.data("shp/period3_2009.shp", "shp/period3_2009.csv", "period3_2009")
res3.2010 <- get.period.data("shp/period3_2010.shp", "shp/period3_2010.csv", "period3_2010")
res3.2011 <- get.period.data("shp/period3_2011.shp", "shp/period3_2011.csv", "period3_2011")
res3.2012 <- get.period.data("shp/period3_2012.shp", "shp/period3_2012.csv", "period3_2012")

H2.prime <- compute.overlap(D = res3$area, G = res2$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H1.prime <- compute.overlap(D = res3$area, G = res1$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)

H1_2006.prime <- compute.overlap(D = res3$area, G = res1.2006$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H1_2007.prime <- compute.overlap(D = res3$area, G = res1.2007$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H1_2008.prime <- compute.overlap(D = res3$area, G = res1.2008$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H1_2009.prime <- compute.overlap(D = res3$area, G = res1.2009$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H1_2010.prime <- compute.overlap(D = res3$area, G = res1.2010$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H1_2011.prime <- compute.overlap(D = res3$area, G = res1.2011$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H1_2012.prime <- compute.overlap(D = res3$area, G = res1.2012$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)

H2_2007.prime <- compute.overlap(D = res3$area, G = res2.2007$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H2_2008.prime <- compute.overlap(D = res3$area, G = res2.2008$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H2_2009.prime <- compute.overlap(D = res3$area, G = res2.2009$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H2_2010.prime <- compute.overlap(D = res3$area, G = res2.2010$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H2_2011.prime <- compute.overlap(D = res3$area, G = res2.2011$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H2_2012.prime <- compute.overlap(D = res3$area, G = res2.2012$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)

H3_2009.prime <- compute.overlap(D = res3$area, G = res3.2009$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H3_2010.prime <- compute.overlap(D = res3$area, G = res3.2010$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H3_2011.prime <- compute.overlap(D = res3$area, G = res3.2011$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H3_2012.prime <- compute.overlap(D = res3$area, G = res3.2012$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)

H2 <- apply(H2.prime, 2, normalize)
H1 <- apply(H1.prime, 2, normalize)

H1_2006 <- apply(H1_2006.prime, 2, normalize)
H1_2007 <- apply(H1_2007.prime, 2, normalize)
H1_2008 <- apply(H1_2008.prime, 2, normalize)
H1_2009 <- apply(H1_2009.prime, 2, normalize)
H1_2010 <- apply(H1_2010.prime, 2, normalize)
H1_2011 <- apply(H1_2011.prime, 2, normalize)
H1_2012 <- apply(H1_2012.prime, 2, normalize)

H2_2007 <- apply(H2_2007.prime, 2, normalize)
H2_2008 <- apply(H2_2008.prime, 2, normalize)
H2_2009 <- apply(H2_2009.prime, 2, normalize)
H2_2010 <- apply(H2_2010.prime, 2, normalize)
H2_2011 <- apply(H2_2011.prime, 2, normalize)
H2_2012 <- apply(H2_2012.prime, 2, normalize)

H3_2009 <- apply(H3_2009.prime, 2, normalize)
H3_2010 <- apply(H3_2010.prime, 2, normalize)
H3_2011 <- apply(H3_2011.prime, 2, normalize)
H3_2012 <- apply(H3_2012.prime, 2, normalize)
