library(stcos)

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
	
	area@data$eZagg <- eZagg
	area@data$esigmavar <- esigmavar
	area@data$Zagg <- Zagg
	area@data$sigmavar <- sigmavar
	
	list(area = area, area.dat = area.dat[idx.keep,], eZagg = eZagg,
		 esigmavar = esigmavar, Zagg = Zagg, sigmavar = sigmavar)
}

setwd("~/Documents/simulations/ST-COS/run-20170109-scott/")

# Load Jon's S matrix (before applying licols!)
S.el <- read.csv("dat/S_sparse.txt.gz", header = FALSE)
S <- sparseMatrix(i = S.el[,1], j = S.el[,2], x = S.el[,3])
rm(S.el)

res1 <- get.period.data("shp/period1.shp", "shp/period1.csv", "period1")
res2 <- get.period.data("shp/period2.shp", "shp/period2.csv", "period2")
res3 <- get.period.data("shp/period3.shp", "shp/period3.csv", "period3")

crs <- CRS("+proj=longlat +lat_1=29.5 +lat_2=45.5 +lat_0=24.395833 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=us-ft +no_defs +ellps=GRS80 +towgs84=0,0,0")
area <- spTransform(res3$area, crs)
D <- area[2,]
plot(D)

cc <- read.csv("dat/knots250_ACS_amr.csv", head = FALSE)


## This is how Jon constructs the knots in Matlab
## knots at each time point
n.cc <- nrow(cc)
level.times = c(
	rep(2012.5, n.cc),
	rep(2011, n.cc),
	rep(2011.5, n.cc),
	rep(2011, n.cc),
	rep(2010, n.cc),
	rep(2010.5, n.cc),
	rep(2010, n.cc),
	rep(2009, n.cc),
	rep(2009.5, n.cc),
	rep(2009, n.cc),
	rep(2008, n.cc),
	rep(2008.5, n.cc),
	rep(2008, n.cc),
	rep(2007, n.cc),
	rep(2007.5, n.cc),
	rep(2007, n.cc),
	rep(2006.5, n.cc),
	rep(2006, n.cc),
	rep(2005.5, n.cc)
)
level <- cbind(cc = cc, gg = level.times)
level2 <- matrix(NA, 0, 0)
level3 <- matrix(NA, 0, 0)

# S1 <- local.bisquare.basis.domain(D, times = 2013, w.s = 0.5, w.t = 0.5, cc = level[,1:2], gg = level[,3], n = 100)
# S2 <- local.bisquare.basis.domain(D, times = 2009:2013, w.s = 0.5, w.t = 0.5, cc, gg = 2012.5, n = 500)
# dim(S1)

S1 <- ArealBi2(D, times = 2013, level, B = 100, srad = 0.5, trad = 0.5, report.period = 1)
