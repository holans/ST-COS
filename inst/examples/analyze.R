library(rgdal)
library(ggplot2)
library(maptools)
library(dplyr)
library(coda)
library(Matrix)

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
	esigmavar <- (area.dat[idx.keep, 2] / 1.645)^2	  # known survey variances
	Zagg <- log(eZagg)
	sigmavar <- esigmavar * (1 / eZagg)^2

	list(area = area, area.dat = area.dat[idx.keep,], eZagg = eZagg,
		 esigmavar = esigmavar, Zagg = Zagg, sigmavar = sigmavar)
}

load("output.Rdata")

# Read data
Z <- as.matrix(read.csv("dat/Zagg.txt.gz", header = FALSE))
sig2eps <- as.matrix(read.csv("dat/sigmavar.txt.gz", header = FALSE))[,1]
C.inv <- as.matrix(read.csv("dat/Kinv.txt.gz", header = FALSE))

H.el <- read.csv("dat/H_sparse.txt.gz", header = FALSE)
H <- sparseMatrix(i = H.el[,1], j = H.el[,2], x = H.el[,3])
rm(H.el)

S.el <- read.csv("dat/S1_sparse.txt.gz", header = FALSE)
S <- sparseMatrix(i = S.el[,1], j = S.el[,2], x = S.el[,3])
rm(S.el)

# Rename variables to be like the Stat paper.
# eta is already named correctly, I think.
R <- length(sig2xi)
mu.mcmc <- xi; rm(xi)
xi.mcmc <- xi2; rm(xi2)
eta.mcmc <- eta; rm(eta)
sig2mu.mcmc <- sig2xi; rm(sig2xi)
sig2K.mcmc <- t(lambda_eta); rm(lambda_eta)
sig2xi.mcmc <- mcmc(matrix(1, R, 1))

# My computer doesn't have enough RAM for this, so just keep
# some of the draws for now
idx <- seq(1, R, 10)
temp <- mu.mcmc[idx,]; mu.mcmc <- temp
temp <- xi.mcmc[idx,]; xi.mcmc <- temp
temp <- eta.mcmc[idx,]; eta.mcmc <- temp
temp <- sig2mu.mcmc[idx,]; sig2mu.mcmc <- temp
temp <- sig2K.mcmc[idx,]; sig2K.mcmc <- temp
temp <- sig2xi.mcmc[idx,]; sig2xi.mcmc <- temp
R <- length(idx)

n <- ncol(mu.mcmc)
r <- ncol(eta.mcmc)
N <- ncol(xi.mcmc)

source.info <- read.table("source_info.dat.gz", sep = ",", head = FALSE)
colnames(source.info) <- c("year", "acs_type", "state", "county")
source.info$county <- sprintf("%03d", source.info$county)
source.info$state <- sprintf("%02d", source.info$state)
source.info$geo_id <- sprintf("0500000US%s%s", source.info$state, source.info$county)

# I think this is the relevant part for the predictions...
idx <- which(source.info$year == 2013 & source.info$acs_type == 3)
source.info.new <- source.info[idx,]

# Load S.new
S.new.el <- read.csv("basis_2013_1_sparse.dat.gz", header = FALSE)
S.new <- sparseMatrix(i = S.new.el[,1], j = S.new.el[,2], x = S.new.el[,3])
rm(S.new.el)
N.new <- nrow(S.new)

# Load H.new
# If we are predicting for 2013 geography, maybe H.new should just be the identity matrix??
H.new <- diag(1, 3109)

# Compute draws for E(Y | Data) for 2013 ACS 3 year
Yhat.mcmc <- matrix(NA, R, N.new)
for (r in 1:R) {
	if (r %% 100 == 0) logger("Computing Y for rep %d\n", r)
	eta <- eta.mcmc[r,]
	mu <- mu.mcmc[r,]
	Yhat.mcmc[r,] <- as.numeric(H.new %*% mu) + as.numeric(S.new %*% eta)
}
Y.post.mean <- apply(Yhat.mcmc, 2, mean)
Y.post.sd <- apply(Yhat.mcmc, 2, sd)

# Need to join posterior results back to source info so it can be plotted on map.
source <- cbind(source.info.new, Y.post.mean, Y.post.sd)

# Read shape files
res1 <- get.period.data("shp/period1.shp", "shp/period1.csv", "period1")
res2 <- get.period.data("shp/period2.shp", "shp/period2.csv", "period2")
res3 <- get.period.data("shp/period3.shp", "shp/period3.csv", "period3")

# Observed data
area <- res3$area
Zagg <- res3$Zagg
sigmavar <- res3$sigmavar
esigmavar <- res3$esigmavar

# http://stackoverflow.com/questions/16462290/obtaining-latitude-and-longitude-with-from-spatial-objects-in-r
# coordinates(area)
# proj4string(area)

# area <- spTransform(area, "+proj=aea +lat_0=24.395833 +lon_0=-96 +x_0=0 +y_0=0 +units=us-ft +ellps=GRS80")
crs <- CRS("+proj=longlat +lat_1=29.5 +lat_2=45.5 +lat_0=24.395833 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=us-ft +no_defs +ellps=GRS80 +towgs84=0,0,0")
area <- spTransform(area, crs)

# Join model output to area
# Make sure to maintain all the shape stuff. Krista and I did before...
# TBD: Assuming that Zagg, sigmavar, and esigmavar come in the same order
# as the counties in the shapefile. CHECK THIS!!
map.dat <- fortify(area, region="GEO_ID")
Z.dat <- data.frame(Zagg = Zagg, sigmavar = sigmavar, id = area@data$GEO_ID)
map.dat$id <- as.character(map.dat$id)
Z.dat$id <- as.character(Z.dat$id)
dat1 <- left_join(map.dat, Z.dat)
dat <- left_join(dat1, source, by = c("id" = "geo_id"))

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

# Plot Y.post.mean
p <- ggplot() +
	geom_polygon(data = dat, aes(x = long, y = lat, group = group,
				 fill = Y.post.mean),
				 color = "black", size = 0.25) +
	theme_bw () +
	ggtitle("Plot of Y.post.mean") +
	scale_fill_gradient(low = "white", high = "red") +
	theme(panel.grid.major = element_blank(),
		  panel.grid.minor = element_blank())
print(p)

# Plot Y.post.sd
p <- ggplot() +
	geom_polygon(data = dat, aes(x = long, y = lat, group = group,
				 fill = Y.post.sd),
				 color = "black", size = 0.25) +
	theme_bw () +
	ggtitle("Plot of Y.post.sd") +
	scale_fill_gradient(low = "white", high = "red") +
	theme(panel.grid.major = element_blank(),
		  panel.grid.minor = element_blank())
print(p)
