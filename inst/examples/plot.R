library(rgdal)
library(ggplot2)
library(maptools)
library(dplyr)

get.period.data <- function(shpfile, datfile, layername)
{
	area.orig <- readOGR(shpfile, layer = layername, verbose = FALSE)
	states <- area.orig@data$STATE

	# Extract the continental U.S.
	idx.HI <- which(states == '02')
	idx.AK <- which(states == '15')
	idx.PR <- which(states == '72')
	idx.keep <- setdiff(1:nrow(area.orig), c(idx.HI, idx.AK, idx.PR))

	area <- area.orig[idx.keep,]
	area.dat <- read.csv(datfile, header = FALSE)
	eZagg <- area.dat[idx.keep, 1]
	esigmavar <- (area.dat[idx.keep, 2] / 1.645)^2		# known survey variances
	Zagg <- log(eZagg)
	sigmavar <- esigmavar * (1 / eZagg)^2

	list(area = area, area.dat = area.dat, eZagg = eZagg,
		 esigmavar = esigmavar, Zagg = Zagg, sigmavar = sigmavar)
}

setwd("~/Documents/simulations/ST-COS/")

res1 <- get.period.data("shp/period1.shp", "shp/period1.csv", "period1")
res2 <- get.period.data("shp/period2.shp", "shp/period2.csv", "period2")
res3 <- get.period.data("shp/period3.shp", "shp/period3.csv", "period3")

area <- res3$area
Zagg <- res3$Zagg
sigmavar <- res3$sigmavar
esigmavar <- res3$esigmavar

# http://stackoverflow.com/questions/16462290/obtaining-latitude-and-longitude-with-from-spatial-objects-in-r
coordinates(area)
proj4string(area)

# area <- spTransform(area, "+proj=aea +lat_0=24.395833 +lon_0=-96 +x_0=0 +y_0=0 +units=us-ft +ellps=GRS80")
crs <- CRS("+proj=longlat +lat_1=29.5 +lat_2=45.5 +lat_0=24.395833 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=us-ft +no_defs +ellps=GRS80 +towgs84=0,0,0")
area <- spTransform(area, crs)

map.dat <- fortify(area, region="GEO_ID")
Z.dat <- data.frame(Zagg = Zagg, sigmavar = sigmavar, id = area@data$GEO_ID)	# CHECK THIS!!!
map.dat$id <- as.character(map.dat$id)
Z.dat$id <- as.character(Z.dat$id)
dat <- left_join(map.dat, Z.dat)

# Plot Zagg
p <- ggplot() +
	geom_polygon(data = dat, aes(x = long, y = lat, group = group, fill = Zagg),
				 color = "black", size = 0.25) +
	theme_bw () +
	ggtitle("Plot of Zagg") +
	scale_fill_gradient(low = "white", high = "red") +
	theme(panel.grid.major = element_blank(),
		  panel.grid.minor = element_blank())
print(p)

# Plot sigmavar
p <- ggplot() +
	geom_polygon(data = dat, aes(x = long, y = lat, group = group, fill = sigmavar),
				 color = "black", size = 0.25) +
	theme_bw () +
	ggtitle("Plot of sigmavar") +
	scale_fill_gradient(low = "white", high = "red") +
	theme(panel.grid.major = element_blank(),
		  panel.grid.minor = element_blank())
print(p)
