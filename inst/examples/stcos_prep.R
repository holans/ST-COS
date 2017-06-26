library(stcos)
library(sf)
library(fields)
library(coda)

set.seed(1234)
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
	# idx.keep <- setdiff(1:nrow(area), c(idx.HI, idx.AK, idx.PR))
	idx.keep <- which(states == '24')

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
# Note: Projecting everything to latlon coordinates leads to a weird error in sf_area.
#	For some reason it's programmed differently for latlon projections.
acs5.2013 <- load.domain("shp/period3.shp", "shp/period3.csv", "period3")
# acs5.2013 <- st_transform(acs5.2013, "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# ----- Set up knots for bisquare basis -----
# Knots must be compatible (same projection) with fine-level geography
u <- st_sample(acs5.2013, size = 1000)
M <- matrix(unlist(u), length(u), 2, byrow = TRUE)
out <- cover.design(M, 250)
knots.sp <- out$design

plot(acs5.2013[,1])
points(knots.sp, pch = 2)

knots.t <- c(2012.5, 2011, 2011.5, 2011, 2010, 2010.5, 2010, 2009, 2009.5, 2009,
			 2008, 2008.5, 2008, 2007, 2007.5, 2007, 2006.5, 2006, 2005.5)
knots <- merge(knots.sp, knots.t)
names(knots) <- c("x", "y", "t")
basis <- SpaceTimeBisquareBasis$new(knots[,1], knots[,2], knots[,3], w.s = 1, w.t = 1)

# Load the domains with observations.
# If necessary, each one can be loaded at a time and added to "sp".
acs1.2013 <- load.domain("shp/period1.shp", "shp/period1.csv", "period1", crs.tx = st_crs(acs5.2013)) # ACS 1-year estimates for 2013
acs3.2013 <- load.domain("shp/period2.shp", "shp/period2.csv", "period2", crs.tx = st_crs(acs5.2013)) # ACS 3-year estimates for 2013
acs1.2012 <- load.domain("shp/period1_2012.shp", "shp/period1_2012.csv", "period1_2012", crs.tx = st_crs(acs5.2013))
acs3.2012 <- load.domain("shp/period2_2012.shp", "shp/period2_2012.csv", "period2_2012", crs.tx = st_crs(acs5.2013))
acs5.2012 <- load.domain("shp/period3_2012.shp", "shp/period3_2012.csv", "period3_2012", crs.tx = st_crs(acs5.2013))
acs1.2011 <- load.domain("shp/period1_2011.shp", "shp/period1_2011.csv", "period1_2011", crs.tx = st_crs(acs5.2013))
acs3.2011 <- load.domain("shp/period2_2011.shp", "shp/period2_2011.csv", "period2_2011", crs.tx = st_crs(acs5.2013))
acs5.2011 <- load.domain("shp/period3_2011.shp", "shp/period3_2011.csv", "period3_2011", crs.tx = st_crs(acs5.2013))
acs1.2010 <- load.domain("shp/period1_2010.shp", "shp/period1_2010.csv", "period1_2010", crs.tx = st_crs(acs5.2013))
acs3.2010 <- load.domain("shp/period2_2010.shp", "shp/period2_2010.csv", "period2_2010", crs.tx = st_crs(acs5.2013))
acs5.2010 <- load.domain("shp/period3_2010.shp", "shp/period3_2010.csv", "period3_2010", crs.tx = st_crs(acs5.2013))
acs1.2009 <- load.domain("shp/period1_2009.shp", "shp/period1_2009.csv", "period1_2009", crs.tx = st_crs(acs5.2013))
acs3.2009 <- load.domain("shp/period2_2009.shp", "shp/period2_2009.csv", "period2_2009", crs.tx = st_crs(acs5.2013))
acs5.2009 <- load.domain("shp/period3_2009.shp", "shp/period3_2009.csv", "period3_2009", crs.tx = st_crs(acs5.2013))
acs1.2008 <- load.domain("shp/period1_2008.shp", "shp/period1_2008.csv", "period1_2008", crs.tx = st_crs(acs5.2013))
acs3.2008 <- load.domain("shp/period2_2008.shp", "shp/period2_2008.csv", "period2_2008", crs.tx = st_crs(acs5.2013))
acs1.2007 <- load.domain("shp/period1_2007.shp", "shp/period1_2007.csv", "period1_2007", crs.tx = st_crs(acs5.2013))
acs3.2007 <- load.domain("shp/period2_2007.shp", "shp/period2_2007.csv", "period2_2007", crs.tx = st_crs(acs5.2013))
acs1.2006 <- load.domain("shp/period1_2006.shp", "shp/period1_2006.csv", "period1_2006", crs.tx = st_crs(acs5.2013))

# Construct a STCOSPrep object, then add space-time domains with observations
sp <- STCOSPrep$new(fine_domain = acs5.2013, fine_domain_geo_name = "GEO_ID", basis = basis, basis_mc_reps = 500)
sp$add_obs(acs1.2013, time = 2013, period = 2013, estimate_name = "EST", variance_name = "VAR", geo_name = "GEO_ID")
sp$add_obs(acs1.2012, time = 2012, period = 2012, estimate_name = "EST", variance_name = "VAR", geo_name = "GEO_ID")
sp$add_obs(acs1.2011, time = 2011, period = 2011, estimate_name = "EST", variance_name = "VAR", geo_name = "GEO_ID")
sp$add_obs(acs1.2010, time = 2010, period = 2010, estimate_name = "EST", variance_name = "VAR", geo_name = "GEO_ID")
sp$add_obs(acs1.2009, time = 2009, period = 2009, estimate_name = "EST", variance_name = "VAR", geo_name = "GEO_ID")
sp$add_obs(acs1.2008, time = 2008, period = 2008, estimate_name = "EST", variance_name = "VAR", geo_name = "GEO_ID")
sp$add_obs(acs1.2007, time = 2007, period = 2007, estimate_name = "EST", variance_name = "VAR", geo_name = "GEO_ID")
sp$add_obs(acs1.2006, time = 2006, period = 2006, estimate_name = "EST", variance_name = "VAR", geo_name = "GEO_ID")
sp$add_obs(acs3.2012, time = 2012, period = 2010:2012, estimate_name = "EST", variance_name = "VAR", geo_name = "GEO_ID")
sp$add_obs(acs3.2011, time = 2011, period = 2009:2011, estimate_name = "EST", variance_name = "VAR", geo_name = "GEO_ID")
sp$add_obs(acs3.2010, time = 2010, period = 2008:2010, estimate_name = "EST", variance_name = "VAR", geo_name = "GEO_ID")
sp$add_obs(acs3.2009, time = 2009, period = 2007:2009, estimate_name = "EST", variance_name = "VAR", geo_name = "GEO_ID")
sp$add_obs(acs3.2008, time = 2008, period = 2006:2008, estimate_name = "EST", variance_name = "VAR", geo_name = "GEO_ID")
sp$add_obs(acs3.2007, time = 2007, period = 2005:2007, estimate_name = "EST", variance_name = "VAR", geo_name = "GEO_ID")
sp$add_obs(acs5.2013, time = 2013, period = 2009:2013, estimate_name = "EST", variance_name = "VAR", geo_name = "GEO_ID")
sp$add_obs(acs5.2012, time = 2012, period = 2008:2012, estimate_name = "EST", variance_name = "VAR", geo_name = "GEO_ID")
sp$add_obs(acs5.2011, time = 2011, period = 2007:2011, estimate_name = "EST", variance_name = "VAR", geo_name = "GEO_ID")
sp$add_obs(acs5.2010, time = 2010, period = 2006:2010, estimate_name = "EST", variance_name = "VAR", geo_name = "GEO_ID")
sp$add_obs(acs5.2009, time = 2009, period = 2005:2009, estimate_name = "EST", variance_name = "VAR", geo_name = "GEO_ID")

Z <- sp$get_Z()
V <- sp$get_V()
H <- sp$get_H()

# Here's where we would think about a reduction for S.
S <- sp$get_S()

# Identity transformation as reduction (no reduction).
sp$set_basis_reduction(identity)

# Rank-revealing QR decomposition to pick out LI columns.
# The same set of columns should be used every time the reduction is needed.
li.out <- licols(as.matrix(S), tol = 0.285)
f <- function(S) { S[,li.out$idx] }
sp$set_basis_reduction(f)

S.reduced <- sp$get_reduced_S()

if (TRUE) {
	# Do a spatial-only basis expansion of fine-domain, and use this as the
	# design matrix to project away from
	sp.basis <- SpatialBisquareBasis$new(knots.sp[,1], knots.sp[,2], w = 1)
	draws.out <- draw_sp_basis_mc(R = 500, domain = acs5.2013, report.period = 100)
	X <- compute_sp_basis_mc(basis = sp.basis, domain = acs5.2013,
		s1 = draws.out$s1, s2 = draws.out$s2, report.period = 100)
	C.inv <- sp$get_Cinv(2005:2013, X)
} else {
	C.inv <- sp$get_Cinv(2005:2013)
}

# ----- Apply Gibbs sampler using MLE as initial value -----
Z.scaled <- (Z - mean(Z)) / sd(Z)
V.scaled <- V / var(Z)
H.scaled <- Diagonal(n = length(Z), x = 1/sd(Z)) %*% H
S.scaled <- Diagonal(n = length(Z), x = 1/sd(Z)) %*%  S.reduced

# Stdized
mle.out <- mle.stcos(Z.scaled, S.scaled, V.scaled, H.scaled, init = list(sig2xi = 1))
init <- list(
	sig2xi = mle.out$sig2xi.hat,
	mu_B = mle.out$mu.hat,
	eta = mle.out$eta.hat
)
gibbs.out <- gibbs.stcos.raw(Z, S.reduced, V, C.inv, H, R = 2000,
	report.period = 1000, burn = 1000, thin = 1, init = init,
	# standardize = list(center = rep(mean(Z),length(Z)), scale = 1/sqrt(V)))
	standardize = list(center = rep(mean(Z),length(Z)), scale = rep(1/sd(Z),length(Z))))

# Unstdized
mle.out <- mle.stcos(Z, S.reduced, V, H, init = list(sig2xi = 100))
init <- list(
	sig2xi = mle.out$sig2xi.hat,
	mu_B = mle.out$mu.hat,
	eta = mle.out$eta.hat
)
gibbs.out <- gibbs.stcos.raw(Z, S.reduced, V, C.inv, H, R = 2000,
	report.period = 1000, burn = 1000, thin = 1, init = init)

mu_B.mcmc <- mcmc(gibbs.out$mu_B.hist)
xi.mcmc <- mcmc(gibbs.out$xi.hist)
eta.mcmc <- mcmc(gibbs.out$eta.hist)
sig2.mcmc <- mcmc(cbind(gibbs.out$sig2mu.hist, gibbs.out$sig2xi.hist, gibbs.out$sig2K.hist))
Y.mcmc <- mcmc(gibbs.out$Y.hist)

plot(mu_B.mcmc[,1:3])
plot(eta.mcmc[,1:3])
plot(xi.mcmc[,1:3])
plot(sig2.mcmc)
plot(Y.mcmc[,1:3])

