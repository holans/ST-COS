#' ---
#' title: Analysis of City of Columbia Neighborhoods
#' author: 
#' date: "Last Updated: `r format(Sys.time(), '%B %d, %Y')`"
#' output:
#'   pdf_document:
#'     number_sections: true
#' ---

# Set so that long lines in R will be wrapped
#+ echo=FALSE, message=FALSE, warning=FALSE
require(knitr)
opts_chunk$set(results = 'asis')

#' # Overview
#' In this example, we are given four neighborhoods in the City of Columbia in
#' Boone County, Missouri. We would like to produce model-based estimates of
#' median household income for these neighborhoods based on 5-year ACS estimates
#' for block-groups in Boone County, Missouri from 2012, 2013, 2014, and 2015.
#' Therefore, the four neighborhoods will be our target supports, and the 2012,
#' 2013, 2014, 2015 year block-groups will be our source supports.
set.seed(1234)

#' # Loading Fine-level Support
#+ message=FALSE
library(jsonlite)
library(sf)
library(tigris)
library(dplyr)
library(ggplot2)

#' Load the fine-level support via shapefile using the `tigris` package.
#' For this, we will use the 2015 block-groups in Boone County, MO.
#' Convert it to an `sf` object, and transform to the projection with EPSG
#' code 3857.
#+ message=FALSE
options(tigris_use_cache = TRUE)
options(tigris_refresh = FALSE)
dom.fine <- block_groups(state = '29', county = '019', year = 2015) %>%
	st_as_sf() %>%
	st_transform(crs = 3857) %>%
	mutate(geoid = GEOID) %>%
	select(-GEOID)

#' A quick plot of the fine-level domain.
ggplot(dom.fine) +
	geom_sf(colour = "black", size = 0.05) +
	ggtitle("Boone County, Missouri") +
	theme_bw()

#' # Assemble Source Supports
#' The following loop pulls ACS direct estimates and associated MOEs from the
#' Census Bureau's Data API, and merges them into a single data frame. For more
#' information about the API, see
#' <https://www.census.gov/data/developers/guidance/api-user-guide.html>.
#' Examples of URLs to pull direct estimates and associated MOEs:
#' 
#' - <https://api.census.gov/data/2015/acs5?get=NAME,B19013_001E&for=block%20group:*&in=state:29+county:019>
#' - <https://api.census.gov/data/2015/acs5?get=NAME,B19013_001M&for=block%20group:*&in=state:29+county:019>
#' 
#' Note: ACS 2016 data is not available via the API, at the time of this writing.

year.levels <- 2012:2015
dat.list <- list()
dat.missing <- list()

for (idx in 1:length(year.levels)) {
	year <- year.levels[idx]

	est_url <- paste('https://api.census.gov/data/', year,
		'/acs5?get=NAME,B19013_001E&for=block%20group:*&in=state:29+county:019',
		sep = '')
	json_data <- fromJSON(est_url)
	est_dat <- data.frame(json_data[-1,])
	colnames(est_dat) <- json_data[1,]

	moe_url <- paste('https://api.census.gov/data/', year,
		'/acs5?get=NAME,B19013_001M&for=block%20group:*&in=state:29+county:019',
		sep = '')
	json_data <- fromJSON(moe_url)
	moe_dat <- data.frame(json_data[-1,])
	colnames(moe_dat) <- json_data[1,]

	my_dat <- est_dat %>%
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
		mutate(DirectVar = (DirectMOE / qnorm(0.95))^2)

	my_shp <- block_groups(state = '29', county = '019', year = year) %>%
		st_as_sf() %>%
		st_transform(crs = 3857)

	my_sf <- my_shp %>%
		inner_join(my_dat, by = c('STATEFP' = 'state', 'COUNTYFP' = 'county',
			'TRACTCE' = 'tract', 'BLKGRPCE' = 'blockgroup')) %>%
		select(geoid = GEOID, state = STATEFP, county = COUNTYFP,
			tract = TRACTCE, blockgroup = BLKGRPCE,
			DirectEst, DirectMOE, DirectVar)

	dat.list[[idx]] <- my_sf %>%
		filter(!is.na(DirectEst))

	dat.missing[[idx]] <- my_sf %>%
		filter(is.na(DirectEst))
}

#' Our assembled source supports are now acs5.2012, ..., acs5.2015
acs5.2012 <- dat.list[[1]]
acs5.2013 <- dat.list[[2]]
acs5.2014 <- dat.list[[3]]
acs5.2015 <- dat.list[[4]]
rm(dat.list)

#' A quick plot of one of the source supports.
ggplot(acs5.2015) +
	geom_sf(colour = "black", size = 0.05, aes(fill = DirectEst)) +
	scale_fill_distiller("DirectEst", palette = "RdYlBu") +
	theme_bw()

#' # Load Target Supports
#' Make sure to transform to the same projection as the fine-level support.
neighbs <- st_read("shp/neighborhoods.shp") %>%
	st_transform(crs = st_crs(dom.fine))

#' A quick plot of the neighborhoods.
ggplot(neighbs) +
	geom_sf(colour = "black", size = 0.05) +
	theme_bw()

#' # Prepare to fit the model
#+ message=FALSE
library(fields)
library(stcos)
library(ggforce)

#' Select spatial knots via space-filling design. Start with uniformly
#' drawn points from the fine-level geography and use `cover.design`
#' to select a subset of those points for knots.
u <- st_sample(dom.fine, size = 2000)
M <- matrix(unlist(u), length(u), 2, byrow = TRUE)
out <- cover.design(M, 500)
knots.sp <- out$design

#' Select temporal knots to be evenly spaced over the years relevant
#' to the source support years.
knots.t <- seq(2008, 2015, by = 0.5)

#' Use Cartesian join of spatial and temporal knots to obtain spatio-temporal knots.
knots <- merge(knots.sp, knots.t)
basis <- SpaceTimeBisquareBasis$new(knots[,1], knots[,2], knots[,3], w.s = 1, w.t = 1)

#' Here is a plot of the spatial knots. We plot a circle around one of the
#' points to illustrate the choice of the `w.s` argument.
knots.sp.dat <- data.frame(x = knots.sp[,1], y = knots.sp[,2], r = basis$get_rl())
g <- ggplot(dom.fine) +
	geom_sf(colour = "black", size = 0.25, fill = NA) +
	geom_point(data = knots.sp.dat, aes(x, y), lwd = 1, col = "red") +
	geom_point(data = knots.sp.dat[1,], aes(x, y), lwd = 3, col = "blue") +
	geom_circle(data = knots.sp.dat[1,], aes(x0=x, y0=y, r=r), fill = NA,
		lwd = 0.5, col = "blue") +
	labs(x = "", y = "") +
	theme_bw()
print(g)

#' # Build terms for STCOS model
#' Create an STCOSPrep object and add source supports
sp <- STCOSPrep$new(fine_domain = dom.fine, fine_domain_geo_name = "geoid",
	basis = basis, basis_mc_reps = 500)
sp$add_obs(acs5.2012, period = 2008:2012, estimate_name = "DirectEst",
	variance_name = "DirectVar", geo_name = "geoid")
sp$add_obs(acs5.2013, period = 2009:2013, estimate_name = "DirectEst",
	variance_name = "DirectVar", geo_name = "geoid")
sp$add_obs(acs5.2014, period = 2010:2014, estimate_name = "DirectEst",
	variance_name = "DirectVar", geo_name = "geoid")
sp$add_obs(acs5.2015, period = 2011:2015, estimate_name = "DirectEst",
	variance_name = "DirectVar", geo_name = "geoid")

#' Read the objects needed for MCMC
z <- sp$get_z()
v <- sp$get_v()
H <- sp$get_H()
S <- sp$get_S()

#' Dimension reduction of `S` matrix via PCA.
eig <- eigen(t(S) %*% S)
rho <- eig$values

idx.S <- which(cumsum(rho) / sum(rho) < 0.65)
Tx.S <- t(eig$vectors[idx.S,])
f <- function(S) { S %*% Tx.S }
sp$set_basis_reduction(f)
S.reduced <- sp$get_reduced_S()

#' Plot the proportion of variation captured by our selection of PCA components.
eigprops <- cumsum(rho) / sum(rho)
plot(eigprops[1:200], xlab = "Dimension", ylab = "Proportion of Variation")
abline(v = max(idx.S), lty = 2)
abline(h = eigprops[max(idx.S)], lty = 2)

#' Pick a covariance structure for random coefficients of basis expansion.
if (FALSE) {
	# Random Walk
	K.inv <- sp$get_Kinv(2008:2015)
} else if (FALSE) {
	# Spatial-only
	K.inv <- sp$get_Kinv(2008:2015, autoreg = FALSE)
} else if (TRUE) {
	# Independence
	K.inv <- diag(x = 1, nrow = ncol(S.reduced))
}

#' Standardize observations before running MCMC.
z.scaled <- (z - mean(z)) / sd(z)
v.scaled <- v / var(z)

#' # Fit the Model
#+ message=FALSE
library(coda)

#' Fit MLE; this will serve as an initial value for MCMC.
K <- solve(K.inv)
mle.out <- mle.stcos(z.scaled, v.scaled, H, S.reduced, K,
	init = list(sig2K = 1, sig2xi = 1))
init <- list(
	sig2K = mle.out$sig2K.hat,
    sig2xi = mle.out$sig2xi.hat,
    mu_B = mle.out$mu.hat
)

#' Run the Gibbs sampler.
gibbs.out <- gibbs.stcos.raw(z.scaled, v.scaled, H, S.reduced, K.inv,
	R = 10000, report.period = 2000, burn = 2000, thin = 10, init = init)
print(gibbs.out)

#' Show some trace plots to assess convergence of the sampler.
plot((mu_B.mcmc <- mcmc(gibbs.out$mu_B.hist))[,1:3])
plot((xi.mcmc <- mcmc(gibbs.out$xi.hist))[,1:3])
plot((eta.mcmc <- mcmc(gibbs.out$eta.hist))[,1:3])

varcomps.mcmc <- mcmc(cbind(
	gibbs.out$sig2mu.hist,
	gibbs.out$sig2xi.hist,
	gibbs.out$sig2K.hist
))
colnames(varcomps.mcmc) <- c("sig2mu", "sig2xi", "sig2K")
plot(varcomps.mcmc)

#' #  Produce Results on target supports
# Compute H and S matrices and get summaries of posterior distribution for E(Y).
# Use 90% significance for all credible intervals and MOEs.

targ.src <- sp$domain2model(acs5.2012, period = 2008:2012, geo_name = "geoid")
E.hat.scaled <- fitted(gibbs.out, targ.src$H, targ.src$S.reduced)
E.hat <- sd(z) * E.hat.scaled + mean(z)                   # Uncenter and unscale
acs5.2012$E.mean <- colMeans(E.hat)                       # Point estimates
acs5.2012$E.sd <- apply(E.hat, 2, sd)                     # SDs
acs5.2012$E.lo <- apply(E.hat, 2, quantile, prob = 0.05)  # Credible interval lo
acs5.2012$E.hi <- apply(E.hat, 2, quantile, prob = 0.95)  # Credible interval hi
acs5.2012$E.median <- apply(E.hat, 2, median)             # Median
acs5.2012$E.moe <- apply(E.hat, 2, sd) * qnorm(0.95)      # MOE

targ.src <- sp$domain2model(acs5.2013, period = 2009:2013, geo_name = "geoid")
E.hat.scaled <- fitted(gibbs.out, targ.src$H, targ.src$S.reduced)
E.hat <- sd(z) * E.hat.scaled + mean(z)                   # Uncenter and unscale
acs5.2013$E.mean <- colMeans(E.hat)                       # Point estimates
acs5.2013$E.sd <- apply(E.hat, 2, sd)                     # SDs
acs5.2013$E.lo <- apply(E.hat, 2, quantile, prob = 0.05)  # Credible interval lo
acs5.2013$E.hi <- apply(E.hat, 2, quantile, prob = 0.95)  # Credible interval hi
acs5.2013$E.median <- apply(E.hat, 2, median)             # Median
acs5.2013$E.moe <- apply(E.hat, 2, sd) * qnorm(0.95)      # MOE

targ.src <- sp$domain2model(acs5.2014, period = 2010:2014, geo_name = "geoid")
E.hat.scaled <- fitted(gibbs.out, targ.src$H, targ.src$S.reduced)
E.hat <- sd(z) * E.hat.scaled + mean(z)                   # Uncenter and unscale
acs5.2014$E.mean <- colMeans(E.hat)                       # Point estimates
acs5.2014$E.sd <- apply(E.hat, 2, sd)                     # SDs
acs5.2014$E.lo <- apply(E.hat, 2, quantile, prob = 0.05)  # Credible interval lo
acs5.2014$E.hi <- apply(E.hat, 2, quantile, prob = 0.95)  # Credible interval hi
acs5.2014$E.median <- apply(E.hat, 2, median)             # Median
acs5.2014$E.moe <- apply(E.hat, 2, sd) * qnorm(0.95)      # MOE

targ.src <- sp$domain2model(acs5.2015, period = 2011:2015, geo_name = "geoid")
E.hat.scaled <- fitted(gibbs.out, targ.src$H, targ.src$S.reduced)
E.hat <- sd(z) * E.hat.scaled + mean(z)                   # Uncenter and unscale
acs5.2015$E.mean <- colMeans(E.hat)                       # Point estimates
acs5.2015$E.sd <- apply(E.hat, 2, sd)                     # SDs
acs5.2015$E.lo <- apply(E.hat, 2, quantile, prob = 0.05)  # Credible interval lo
acs5.2015$E.hi <- apply(E.hat, 2, quantile, prob = 0.95)  # Credible interval hi
acs5.2015$E.median <- apply(E.hat, 2, median)             # Median
acs5.2015$E.moe <- apply(E.hat, 2, sd) * qnorm(0.95)      # MOE

targ.neighbs <- sp$domain2model(neighbs, period = 2011:2015, geo_name = "Region")
E.hat.scaled <- fitted(gibbs.out, targ.neighbs$H, targ.neighbs$S.reduced)
E.hat <- sd(z) * E.hat.scaled + mean(z)              # Uncenter and unscale
neighbs$E.median <- apply(E.hat, 2, median)
neighbs$E.moe <- apply(E.hat, 2, sd) * qnorm(0.95)
neighbs$E.mean <- apply(E.hat, 2, mean)
neighbs$E.sd <- apply(E.hat, 2, sd)
neighbs$E.lo <- apply(E.hat, 2, quantile, prob = 0.05)
neighbs$E.hi <- apply(E.hat, 2, quantile, prob = 0.95)

#' The objective of our analysis - predictions on the four target neighborhoods.
print(neighbs)

#' # Plot Results
#+ message=FALSE
library(gridExtra)
library(ggrepel)

#' Maps of direct and model-based 2015 5-year estimates.
lim.est <- range(acs5.2015$DirectEst, acs5.2015$E.mean)
g <- ggplot(acs5.2015) +
	geom_sf(colour = "black", size = 0.05, aes(fill = DirectEst)) +
	ggtitle("Median Household Income\nfor Boone County",
		subtitle = "2015 5yr ACS Direct Estimates") +
	scale_fill_distiller("DirectEst", palette = "RdYlBu", limits = lim.est) +
	theme_bw()
h <- ggplot(acs5.2015) +
	geom_sf(colour = "black", size = 0.05, aes(fill = E.mean)) +
	ggtitle("Median Household Income\nfor Boone County",
		subtitle = "2015 5yr Model Estimates") +
	scale_fill_distiller("E.mean", palette = "RdYlBu", limits = lim.est) +
	theme_bw()
k <- grid.arrange(g,h, ncol = 2)

#' Scatter plots comparing direct and model-based 5-year estimates for
#' 2012, ..., 2015.
g2012 <- ggplot(acs5.2012, aes(x=DirectEst, y=E.mean)) +
	geom_point(size = 2) +
	geom_abline(intercept = 0, slope = 1, color="red",
		linetype="dashed", size=1.2) +
	ggtitle("2012 5yr ACS Direct Estimates") +
	labs(x = "Direct Estimate", y = "Model-Based Estimate") +
	theme_bw()
g2013 <- ggplot(acs5.2013, aes(x=DirectEst, y=E.mean)) +
	geom_point(size = 2) +
	geom_abline(intercept = 0, slope = 1, color="red",
		linetype="dashed", size=1.2) +
	ggtitle("2013 5yr ACS Direct Estimates") +
	labs(x = "Direct Estimate", y = "Model-Based Estimate") +
	theme_bw()
g2014 <- ggplot(acs5.2014, aes(x=DirectEst, y=E.mean)) +
	geom_point(size = 2) +
	geom_abline(intercept = 0, slope = 1, color="red",
		linetype="dashed", size=1.2) +
	ggtitle("2014 5yr ACS Direct Estimates") +
	labs(x = "Direct Estimate", y = "Model-Based Estimate") +
	theme_bw()
g2015 <- ggplot(acs5.2015, aes(x=DirectEst, y=E.mean)) +
	geom_point(size = 2) +
	geom_abline(intercept = 0, slope = 1, color="red",
		linetype="dashed", size=1.2) +
	ggtitle("2015 5yr ACS Direct Estimates") +
	labs(x = "Direct Estimate", y = "Model-Based Estimate") +
	theme_bw()
k <- grid.arrange(g2012, g2013, g2014, g2015, nrow = 2, ncol = 2)

#' Plot neighborhood areas (target supports) among ACS 5-year direct estimates;
#' this gives a sense of whether the model-based esimtates are reasonable.
#' This map takes a bit of preparation.

Central <- neighbs[1,]
East <- neighbs[2,]
North <- neighbs[3,]
Paris <- neighbs[4,]
Outlier <- acs5.2015[idx,]
Missing1 <- dat.missing[[4]][1,]
Missing2 <- dat.missing[[4]][2,]
Missing3 <- dat.missing[[4]][3,]
Missing4 <- dat.missing[[4]][4,]

# Prevent `sf` package warnings like "st_centroid assumes attributes are
# constant over geometries of x"
st_agr(Central) <- "constant"
st_agr(East) <- "constant"
st_agr(North) <- "constant"
st_agr(Paris) <- "constant"
st_agr(Outlier) <- "constant"
st_agr(Missing1) <- "constant"
st_agr(Missing2) <- "constant"
st_agr(Missing3) <- "constant"
st_agr(Missing4) <- "constant"

Central.coord <- st_coordinates(st_centroid(Central))
East.coord <- st_coordinates(st_centroid(East))
North.coord <- st_coordinates(st_centroid(North))
Paris.coord <- st_coordinates(st_centroid(Paris))
Outlier.coord <- st_coordinates(st_centroid(Outlier))
Missing1.coord <- st_coordinates(st_centroid(Missing1))
Missing2.coord <- st_coordinates(st_centroid(Missing2))
Missing3.coord <- st_coordinates(st_centroid(Missing3))
Missing4.coord <- st_coordinates(st_centroid(Missing4))

g <- ggplot(acs5.2015) +
	geom_sf(colour = "black", size = 0.05, aes(fill = DirectEst)) +
	ggtitle("Median Household Income for Boone County",
		subtitle = "ACS 2015 5yr Direct Estimates") +
	scale_fill_distiller("DirectEst", palette = "RdYlBu", limits = lim.est) +
	geom_sf(data = neighbs, fill = "black") +
	geom_label_repel(data = st_centroid(East), nudge_x = 20000, nudge_y = 130000,
		aes(x=East.coord[1], y=East.coord[2], label="East")) +
	geom_label_repel(data = st_centroid(Central), nudge_x = -10000, nudge_y = 130000,
		aes(x=Central.coord[1], y=Central.coord[2], label="Central")) +
	geom_label_repel(data = st_centroid(North), nudge_x = 0, nudge_y = 100000,
		aes(x=North.coord[1], y=North.coord[2], label="North")) +
	geom_label_repel(data = st_centroid(Paris), nudge_x = 10000, nudge_y = 100000,
		aes(x=Paris.coord[1], y=Paris.coord[2], label="Paris")) +
	geom_label_repel(data = st_centroid(Outlier), nudge_x = 0, nudge_y = -50000,
		aes(x=Outlier.coord[1], y=Outlier.coord[2], label="Outlier")) +
	geom_label_repel(data = st_centroid(Missing1), nudge_x = 100000, nudge_y = -50000,
		aes(x=Missing1.coord[1], y=Missing1.coord[2], label="Missing1")) +
	geom_label_repel(data = st_centroid(Missing2), nudge_x = 100000, nudge_y = -36000,
		aes(x=Missing2.coord[1], y=Missing2.coord[2], label="Missing2")) +
	geom_label_repel(data = st_centroid(Missing3), nudge_x = -100000, nudge_y = -50000,
		aes(x=Missing3.coord[1], y=Missing3.coord[2], label="Missing3")) +
	geom_label_repel(data = st_centroid(Missing4), nudge_x = -100000, nudge_y = -38000,
		aes(x=Missing4.coord[1], y=Missing4.coord[2], label="Missing4")) +
	xlab(NULL) +
	ylab(NULL) +
	theme_bw()
print(g)
