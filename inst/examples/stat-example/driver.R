library(stcos)
library(MASS)

get.period.data <- function(shpfile, datfile, layername)
{
	area.orig <- readOGR(shpfile, layer = layername, verbose = FALSE)
	states <- area.orig@data$STATE
	
	# Extract the continental U.S.
	idx.HI <- which(states == '02')
	idx.AK <- which(states == '15')
	idx.PR <- which(states == '72')
	idx.keep <- setdiff(1:nrow(area.orig), c(idx.HI, idx.AK, idx.PR))
	
	crs <- CRS("+proj=longlat +lat_1=29.5 +lat_2=45.5 +lat_0=24.395833 +lon_0=-96 +x_0=0 +y_0=0 +datum=WGS84 +units=us-ft +no_defs +ellps=GRS80 +towgs84=0,0,0")
	area.orig <- spTransform(area.orig, crs)

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

# ----- Read all the spatio-temporal domains -----
county1 <- get.period.data("shp/period1.shp", "shp/period1.csv", "period1")
county2 <- get.period.data("shp/period2.shp", "shp/period2.csv", "period2")
county3 <- get.period.data("shp/period3.shp", "shp/period3.csv", "period3")
county1_2012 <- get.period.data("shp/period1_2012.shp", "shp/period1_2012.csv", "period1_2012")
county2_2012 <- get.period.data("shp/period2_2012.shp", "shp/period2_2012.csv", "period2_2012")
county3_2012 <- get.period.data("shp/period3_2012.shp", "shp/period3_2012.csv", "period3_2012")
county1_2011 <- get.period.data("shp/period1_2011.shp", "shp/period1_2011.csv", "period1_2011")
county2_2011 <- get.period.data("shp/period2_2011.shp", "shp/period2_2011.csv", "period2_2011")
county3_2011 <- get.period.data("shp/period3_2011.shp", "shp/period3_2011.csv", "period3_2011")
county1_2010 <- get.period.data("shp/period1_2010.shp", "shp/period1_2010.csv", "period1_2010")
county2_2010 <- get.period.data("shp/period2_2010.shp", "shp/period2_2010.csv", "period2_2010")
county3_2010 <- get.period.data("shp/period3_2010.shp", "shp/period3_2010.csv", "period3_2010")
county1_2009 <- get.period.data("shp/period1_2009.shp", "shp/period1_2009.csv", "period1_2009")
county2_2009 <- get.period.data("shp/period2_2009.shp", "shp/period2_2009.csv", "period2_2009")
county3_2009 <- get.period.data("shp/period3_2009.shp", "shp/period3_2009.csv", "period3_2009")
county1_2008 <- get.period.data("shp/period1_2008.shp", "shp/period1_2008.csv", "period1_2008")
county2_2008 <- get.period.data("shp/period2_2008.shp", "shp/period2_2008.csv", "period2_2008")
county1_2007 <- get.period.data("shp/period1_2007.shp", "shp/period1_2007.csv", "period1_2007")
county2_2007 <- get.period.data("shp/period2_2007.shp", "shp/period2_2007.csv", "period2_2007")
county1_2006 <- get.period.data("shp/period1_2006.shp", "shp/period1_2006.csv", "period1_2006")

# ----- Set up spatio-temporal knots for computing bisquare basis
# Read spatial knots
cc <- read.csv("dat/knots250_ACS_amr.csv", head = FALSE)

# This is how Jon constructs the matrix of spatio-temporal knots in Matlab
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

# ----- Compute the full S matrix -----
# Do it as close as possible to Jon's code
SGBF1 <- ArealBi2(county1$area, times = 2013, level1 = level, B = 100, srad = 0.5, trad = 0.5)
SGBF3 <- ArealBi2(county3$area, times = 2009:2013, level1 = level, B = 100, srad = 0.5, trad = 0.5)
SGBF1_2012 <- ArealBi2(county1_2012$area, times = 2012, level1 = level, B = 100, srad = 0.5, trad = 0.5)
SGBF2_2012 <- ArealBi2(county2_2012$area, times = 2010:2012, level1 = level, B = 100, srad = 0.5, trad = 0.5)
SGBF3_2012 <- ArealBi2(county3_2012$area, times = 2008:2012, level1 = level, B = 100, srad = 0.5, trad = 0.5)
SGBF1_2011 <- ArealBi2(county1_2011$area, times = 2011, level1 = level, B = 100, srad = 0.5, trad = 0.5)
SGBF2_2011 <- ArealBi2(county2_2011$area, times = 2009:2011, level1 = level, B = 100, srad = 0.5, trad = 0.5)
SGBF3_2011 <- ArealBi2(county3_2011$area, times = 2007:2011, level1 = level, B = 100, srad = 0.5, trad = 0.5)
SGBF1_2010 <- ArealBi2(county1_2010$area, times = 2010, level1 = level, B = 100, srad = 0.5, trad = 0.5)
SGBF2_2010 <- ArealBi2(county2_2010$area, times = 2008:2010, level1 = level, B = 100, srad = 0.5, trad = 0.5)
SGBF3_2010 <- ArealBi2(county3_2010$area, times = 2006:2010, level1 = level, B = 100, srad = 0.5, trad = 0.5)
SGBF1_2009 <- ArealBi2(county1_2009$area, times = 2009, level1 = level, B = 100, srad = 0.5, trad = 0.5)
SGBF2_2009 <- ArealBi2(county2_2009$area, times = 2007:2009, level1 = level, B = 100, srad = 0.5, trad = 0.5)
SGBF3_2009 <- ArealBi2(county3_2009$area, times = 2005:2009, level1 = level, B = 100, srad = 0.5, trad = 0.5)
SGBF1_2008 <- ArealBi2(county1_2008$area, times = 2008, level1 = level, B = 100, srad = 0.5, trad = 0.5)
SGBF2_2008 <- ArealBi2(county2_2008$area, times = 2006:2008, level1 = level, B = 100, srad = 0.5, trad = 0.5)
SGBF1_2007 <- ArealBi2(county1_2007$area, times = 2007, level1 = level, B = 100, srad = 0.5, trad = 0.5)
SGBF2_2007 <- ArealBi2(county2_2007$area, times = 2005:2007, level1 = level,100, srad = 0.5, trad = 0.5)
SGBF1_2006 <- ArealBi2(county1_2006$area, times = 2006, level1 = level, B = 100, srad = 0.5, trad = 0.5)
S <- rbind(SGBF1, SGBF1_2012, SGBF1_2011, SGBF1_2010, SGBF1_2009, SGBF1_2008,
	SGBF1_2007, SGBF1_2006, SGBF2_2012, SGBF2_2011, SGBF2_2010, SGBF2_2009,
	SGBF2_2008, SGBF2_2007, SGBF3, SGBF3_2012, SGBF3_2011, SGBF3_2010, SGBF3_2009)

# TBD: Extract submatrix of interest and save.
# [S1,idx]=licols(S);


# ----- Compute and save the H matrix -----
H2.prime <- compute.overlap(D = county3$area, G = county2$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H1.prime <- compute.overlap(D = county3$area, G = county1$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H1_2006.prime <- compute.overlap(D = county3$area, G = county1_2006$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H1_2007.prime <- compute.overlap(D = county3$area, G = county1_2007$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H1_2008.prime <- compute.overlap(D = county3$area, G = county1_2008$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H1_2009.prime <- compute.overlap(D = county3$area, G = county1_2009$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H1_2010.prime <- compute.overlap(D = county3$area, G = county1_2010$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H1_2011.prime <- compute.overlap(D = county3$area, G = county1_2011$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H1_2012.prime <- compute.overlap(D = county3$area, G = county1_2012$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H2_2007.prime <- compute.overlap(D = county3$area, G = county2_2007$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H2_2008.prime <- compute.overlap(D = county3$area, G = county2_2008$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H2_2009.prime <- compute.overlap(D = county3$area, G = county2_2009$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H2_2010.prime <- compute.overlap(D = county3$area, G = county2_2010$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H2_2011.prime <- compute.overlap(D = county3$area, G = county2_2011$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H2_2012.prime <- compute.overlap(D = county3$area, G = county2_2012$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H3_2009.prime <- compute.overlap(D = county3$area, G = county3_2009$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H3_2010.prime <- compute.overlap(D = county3$area, G = county3_2010$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H3_2011.prime <- compute.overlap(D = county3$area, G = county3_2011$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)
H3_2012.prime <- compute.overlap(D = county3$area, G = county3_2012$area, denomflag = 2, Dxyflag = 1, Gxyflag = 1, report.period = 100)

# This normalization may not be necessary...
# But if it is, make sure we produce sparse matrices
H2 <- H2.prime / Matrix(colSums(H2.prime), nrow(H2.prime), ncol(H2.prime))
H1 <- H1.prime / Matrix(colSums(H1.prime), nrow(H1.prime), ncol(H1.prime))
H1_2006 <- H1_2006.prime / Matrix(colSums(H1_2006.prime), nrow(H1_2006.prime), ncol(H1_2006.prime))
H1_2007 <- H1_2007.prime / Matrix(colSums(H1_2007.prime), nrow(H1_2007.prime), ncol(H1_2007.prime))
H1_2008 <- H1_2008.prime / Matrix(colSums(H1_2008.prime), nrow(H1_2008.prime), ncol(H1_2008.prime))
H1_2009 <- H1_2009.prime / Matrix(colSums(H1_2009.prime), nrow(H1_2009.prime), ncol(H1_2009.prime))
H1_2010 <- H1_2010.prime / Matrix(colSums(H1_2010.prime), nrow(H1_2010.prime), ncol(H1_2010.prime))
H1_2011 <- H1_2011.prime / Matrix(colSums(H1_2011.prime), nrow(H1_2011.prime), ncol(H1_2011.prime))
H1_2012 <- H1_2012.prime / Matrix(colSums(H1_2012.prime), nrow(H1_2012.prime), ncol(H1_2012.prime))
H2_2007 <- H2_2007.prime / Matrix(colSums(H2_2007.prime), nrow(H2_2007.prime), ncol(H2_2007.prime))
H2_2008 <- H2_2008.prime / Matrix(colSums(H2_2008.prime), nrow(H2_2008.prime), ncol(H2_2008.prime))
H2_2009 <- H2_2009.prime / Matrix(colSums(H2_2009.prime), nrow(H2_2009.prime), ncol(H2_2009.prime))
H2_2010 <- H2_2010.prime / Matrix(colSums(H2_2010.prime), nrow(H2_2010.prime), ncol(H2_2010.prime))
H2_2011 <- H2_2011.prime / Matrix(colSums(H2_2011.prime), nrow(H2_2011.prime), ncol(H2_2011.prime))
H2_2012 <- H2_2012.prime / Matrix(colSums(H2_2012.prime), nrow(H2_2012.prime), ncol(H2_2012.prime))
H3_2009 <- H3_2009.prime / Matrix(colSums(H3_2009.prime), nrow(H3_2009.prime), ncol(H3_2009.prime))
H3_2010 <- H3_2010.prime / Matrix(colSums(H3_2010.prime), nrow(H3_2010.prime), ncol(H3_2010.prime))
H3_2011 <- H3_2011.prime / Matrix(colSums(H3_2011.prime), nrow(H3_2011.prime), ncol(H3_2011.prime))
H3_2012 <- H3_2012.prime / Matrix(colSums(H3_2012.prime), nrow(H3_2012.prime), ncol(H3_2012.prime))

H <- t(cbind(H1, H1_2012, H1_2011, H1_2010, H1_2009, H1_2008, H1_2007, H1_2006,
	H2_2012, H2_2011, H2_2010, H2_2009, H2_2008, H2_2007, diag(1,3109), H3_2012,
	H3_2011, H3_2010, H3_2009))

# ----- Moran's I Propagator -----
# AMR notes: GBI is almost exactly the identity matrix...
# R and Matlab return different M... the decompisition isn't supposed to be unique, is it?
Adj <- gTouches(county3$area, byid=TRUE)
for (j in 1:nrow(Adj)) {
	ss <- sum(countAdj[j,])
	if (ss > 0) {
		Adj[j,] <- Adj[j,] / ss
	}
}
Q <- diag(1,nrow(Adj)) - 0.9*Adj
Qinv <- ginv(Q)
eig.Q <- eigen(Q)

B <- cbind(t(eig.Q$vectors) %*% matrix(1, 3109, 1), diag(1,3109))
GBI <- B %*% ginv(t(B) %*% B) %*% t(B)
eig.GBI <- eigen(GBI)
M <- Re(eig.GBI$vectors)

# ----- Compute and save the C.inv matrix -----
Sconnector1 <- ArealBi2(county3$area, times = 2005, level1 = level, B = 100, srad = 0.5, trad = 0.5)
Sconnector2 <- ArealBi2(county3$area, times = 2006, level1 = level, B = 100, srad = 0.5, trad = 0.5)
Sconnector3 <- ArealBi2(county3$area, times = 2007, level1 = level, B = 100, srad = 0.5, trad = 0.5)
Sconnector4 <- ArealBi2(county3$area, times = 2008, level1 = level, B = 100, srad = 0.5, trad = 0.5)
Sconnector5 <- ArealBi2(county3$area, times = 2009, level1 = level, B = 100, srad = 0.5, trad = 0.5)
Sconnector6 <- ArealBi2(county3$area, times = 2010, level1 = level, B = 100, srad = 0.5, trad = 0.5)
Sconnector7 <- ArealBi2(county3$area, times = 2011, level1 = level, B = 100, srad = 0.5, trad = 0.5)
Sconnector8 <- ArealBi2(county3$area, times = 2012, level1 = level, B = 100, srad = 0.5, trad = 0.5)
Sconnector9 <- ArealBi2(county3$area, times = 2013, level1 = level, B = 100, srad = 0.5, trad = 0.5)
Sconnector = rbind(Sconnector1, Sconnector2, Sconnector3, Sconnector4,
	Sconnector5, Sconnector6, Sconnector7, Sconnector8, Sconnector9)
Sconnectorf <- Sconnector[,idx];

Kinv <- make_full_model_sptcovar_9(Q, M, Sconnectorf, 3109)
[P,D] <- eigen(Kinv)
D <- real(diag(D))
D(D<0) <- 0
Dinv <- D
Dinv(D>0) <- (1/D(D>0))
K <- real(P) %*% diag(Dinv) %*% real(P')
Kinv <- real(P) %*% diag(D) %*% real(P')

# ----- Compute and save Z -----
Zagg <- c(
	county1$area$eZagg,
	county1_2012$area$eZagg,
	county1_2011$area$eZagg,
	county1_2010$area$eZagg,
	county1_2009$area$eZagg,
	county1_2008$area$eZagg,
	county1_2007$area$eZagg,
	county1_2006$area$eZagg,
	county2_2012$area$eZagg,
	county2_2011$area$eZagg,
	county2_2010$area$eZagg,
	county2_2009$area$eZagg,
	county2_2008$area$eZagg,
	county2_2007$area$eZagg,
	county3$area$eZagg,
	county3_2012$area$eZagg,
	county3_2011$area$eZagg,
	county3_2010$area$eZagg,
	county3_2009$area$eZagg)

# ----- Compute and save sigmavar -----
sigmavar <- c(
	county1$area$esigmavar,
	county1_2012$area$esigmavar,
	county1_2011$area$esigmavar,
	county1_2010$area$esigmavar,
	county1_2009$area$esigmavar,
	county1_2008$area$esigmavar,
	county1_2007$area$esigmavar,
	county1_2006$area$esigmavar,
	county2_2012$area$esigmavar,
	county2_2011$area$esigmavar,
	county2_2010$area$esigmavar,
	county2_2009$area$esigmavar,
	county2_2008$area$esigmavar,
	county2_2007$area$esigmavar,
	county3$area$esigmavar,
	county3_2012$area$esigmavar,
	county3_2011$area$esigmavar,
	county3_2010$area$esigmavar,
	county3_2009$area$esigmavar)

save.image("prep.Rdata")

