if (FALSE) {
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
}

if (FALSE) {
	res1 <- get.period.data("shp/period1.shp", "shp/period1.csv", "period1")
	res2 <- get.period.data("shp/period2.shp", "shp/period2.csv", "period2")
	res3 <- get.period.data("shp/period3.shp", "shp/period3.csv", "period3")
	
	source.domains <- list(res1$area, res2$area, res3$area)
	fine.domain <- res3$area
}

# This function should produce everything needed to run the model
# o H: matrix of overlaps
# o S: matrix of basis function
# o C.inv: the solution, up to a constant, of the Higham approximation
create.domain <- function(source.domains, fine.domain, num.basis)
{
	stopifnot(is.list(source.domains))
	cl <- unlist(lapply(source.domains, class))
	stopifnot(prod(cl == "SpatialPolygonsDataFrame") == 1)
	L <- length(source.domains)

	d <- unlist(lapply(source.domains, nrow))
	H.list <- list()
	S.list <- list()

	logger("Computing overlaps\n")
	for (l in 1:L) {
		logger("Computing overlaps for source domain %d of %d\n", l, L)
		H.list[[l]] <- compute.overlap(fine.domain, source.domains[[l]], denomflag = 2,
			Dxyflag = 1, Gxyflag = 1, n_mc = 100, report.period = 10)
	}
	logger("Finished computing overlaps\n")
	
	logger("Computing basis functions\n")
	for (l in 1:L) {
		logger("Computing basis functions for source domain %d of %d\n", l, L)
		S.list[[l]] <- diag(1, d[l], d[l])
	}
	logger("Finished basis functions\n")
}
