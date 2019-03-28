#' ---
#' title: Analysis of City of Columbia Neighborhoods
#' author: 
#' date: "Last Updated: `r format(Sys.time(), '%B %d, %Y')`"
#' output:
#'   html_document:
#'     number_sections: true
#' ---

#+ echo = FALSE
options(width = 80)

#' # Overview
#' In this example, we are given four neighborhoods in the City of Columbia in
#' Boone County, Missouri. We would like to produce model-based estimates of
#' median household income for these neighborhoods based on 5-year ACS estimates
#' for block-groups in Boone County, Missouri from years 2013, 2014, ..., 2017.
#' Therefore, the four neighborhoods will be our target supports, and the 2013 -
#' 2017 year block-groups will be our source supports.
set.seed(1234)

#' Load the fine-level support via shapefile using the `tigris` package.
#' For this, we will use the 2017 block-groups in Boone County, MO.
#' Convert it to an `sf` object, and transform to the projection with EPSG
#' code 3857.
#+ message=FALSE
library(tigris)
library(ggplot2)
library(dplyr)
options(tigris_use_cache = TRUE)
options(tigris_refresh = FALSE)
dom.fine <- block_groups(state = '29', county = '019', year = 2017) %>%
	st_as_sf() %>%
	st_transform(crs = 3857) %>%
	mutate(geoid = GEOID) %>%
	select(-GEOID)

#' A quick plot of the fine-level domain.
ggplot(dom.fine) +
	geom_sf(colour = "black", size = 0.05) +
	ggtitle("Boone County, Missouri") +
	theme_bw()

#' Load previously constructed source supports.
data(acs_sf)

#' A quick plot of one of the source supports.
ggplot(acs5_2017) +
	geom_sf(colour = "black", size = 0.05, aes(fill = DirectEst)) +
	scale_fill_distiller("DirectEst", palette = "RdYlBu") +
	theme_bw()

#' # Load Target Supports
#' Make sure to transform to the same projection as the fine-level support.
data(columbia_neighbs)
neighbs <- columbia_neighbs %>%
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
knots.t <- seq(2009, 2017, by = 0.5)

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
sp$add_obs(acs5_2013, period = 2009:2013, estimate_name = "DirectEst",
	variance_name = "DirectVar", geo_name = "geoid")
sp$add_obs(acs5_2014, period = 2010:2014, estimate_name = "DirectEst",
	variance_name = "DirectVar", geo_name = "geoid")
sp$add_obs(acs5_2015, period = 2011:2015, estimate_name = "DirectEst",
	variance_name = "DirectVar", geo_name = "geoid")
sp$add_obs(acs5_2016, period = 2012:2016, estimate_name = "DirectEst",
	variance_name = "DirectVar", geo_name = "geoid")
sp$add_obs(acs5_2017, period = 2013:2017, estimate_name = "DirectEst",
	variance_name = "DirectVar", geo_name = "geoid")

#' Read the objects needed for MCMC
z <- sp$get_z()
v <- sp$get_v()
H <- sp$get_H()
S <- sp$get_S()
idx.missing <- which(is.na(z))
idx.nonmissing <- which(!is.na(z))

#' Dimension reduction of `S` matrix via PCA.
eig <- eigen(t(S) %*% S)
rho <- eig$values

idx.S <- which(cumsum(rho) / sum(rho) < 0.65)
Tx <- eig$vectors[,idx.S]
f <- function(S) { S %*% Tx }
sp$set_basis_reduction(f)
S.reduced <- sp$get_reduced_S()

#' Plot the proportion of variation captured by our selection of PCA components.
eigprops <- cumsum(rho) / sum(rho)
plot(eigprops[1:200], xlab = "Dimension", ylab = "Proportion of Variation")
abline(v = max(idx.S), lty = 2)
abline(h = eigprops[max(idx.S)], lty = 2)

#' Pick a covariance structure for random coefficients of basis expansion.
if (FALSE) {
	# Moran’s I Basis
	K.inv <- sp$get_Kinv(2009:2017, method = "moran")
} else if (FALSE) {
	# Random Walk
	K.inv <- sp$get_Kinv(2009:2017, method = "randomwalk")
} else if (FALSE) {
	# Spatial-only (CAR)
	K.inv <- sp$get_Kinv(2009:2017, method = "car")
} else if (TRUE) {
	# Independence
	K.inv <- sp$get_Kinv(2009:2017, method = "independence")
}

#' Standardize observations before running MCMC.
z.mean <- mean(z, na.rm = TRUE)
z.sd <- sd(z, na.rm = TRUE)
z.scaled <- (z - z.mean) / z.sd
v.scaled <- v / z.sd^2

#' # Fit the Model
#+ message=FALSE
library(coda)

#' Fit MLE; this will serve as an initial value for MCMC.
K <- solve(K.inv)
mle.out <- mle.stcos(
	z = z.scaled[idx.nonmissing],
	v = v.scaled[idx.nonmissing],
	H = H[idx.nonmissing,],
	S = S.reduced[idx.nonmissing,],
	K = K,
	init = list(sig2K = 1, sig2xi = 1)
)
init <- list(
	sig2K = mle.out$sig2K.hat,
    sig2xi = mle.out$sig2xi.hat,
    mu_B = mle.out$mu.hat
)

#' Run the Gibbs sampler.
gibbs.out <- gibbs.stcos.raw(
	z = z.scaled[idx.nonmissing],
	v = v.scaled[idx.nonmissing],
	H = H[idx.nonmissing,],
	S = S.reduced[idx.nonmissing,],
	K.inv = K.inv,
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
#' Compute `H` and `S` matrices and get summaries of posterior distribution for E(Y).
#' Use 90% significance for all credible intervals and MOEs.
append_results <- function(dat_sf, period, geo_name, alpha = 0.10) {
	out <- sp$domain2model(dat_sf, period = period, geo_name = geo_name)
	E.hat.scaled <- fitted(gibbs.out, out$H, out$S.reduced)
	E.hat <- z.sd * E.hat.scaled + z.mean                      # Uncenter and unscale
	dat_sf$E.mean <- colMeans(E.hat)                           # Point estimates
	dat_sf$E.sd <- apply(E.hat, 2, sd)                         # SDs
	dat_sf$E.lo <- apply(E.hat, 2, quantile, prob = alpha/2)   # Credible interval lo
	dat_sf$E.hi <- apply(E.hat, 2, quantile, prob = 1-alpha/2) # Credible interval hi
	dat_sf$E.median <- apply(E.hat, 2, median)                 # Median
	dat_sf$E.moe <- apply(E.hat, 2, sd) * qnorm(1-alpha/2)     # MOE
	return(dat_sf)
}

acs5_2013 <- append_results(acs5_2013, period = 2009:2013, geo_name = "geoid")
acs5_2014 <- append_results(acs5_2014, period = 2010:2014, geo_name = "geoid")
acs5_2015 <- append_results(acs5_2015, period = 2011:2015, geo_name = "geoid")
acs5_2016 <- append_results(acs5_2016, period = 2012:2016, geo_name = "geoid")
acs5_2017 <- append_results(acs5_2017, period = 2013:2017, geo_name = "geoid")
neighbs <- append_results(neighbs, 2013:2017, geo_name = "Region")

#' The objective of our analysis - predictions on the four target neighborhoods.
print(neighbs)

#' # Plot Results
#+ message=FALSE
library(gridExtra)
library(ggrepel)

#' Maps of direct and model-based 2017 5-year estimates.
lim.est <- range(acs5_2017$DirectEst, acs5_2017$E.mean)
g <- ggplot(acs5_2017) +
	geom_sf(colour = "black", size = 0.05, aes(fill = DirectEst)) +
	ggtitle("Median Household Income\nfor Boone County",
		subtitle = "2017 5yr ACS Direct Estimates") +
	scale_fill_distiller("DirectEst", palette = "RdYlBu", limits = lim.est) +
	theme_bw()
h <- ggplot(acs5_2017) +
	geom_sf(colour = "black", size = 0.05, aes(fill = E.mean)) +
	ggtitle("Median Household Income\nfor Boone County",
		subtitle = "2017 5yr Model Estimates") +
	scale_fill_distiller("E.mean", palette = "RdYlBu", limits = lim.est) +
	theme_bw()
k <- grid.arrange(g,h, ncol = 2)

#' Scatter plots comparing direct and model-based 5-year estimates for
#' 2013, ..., 2017.
scatter.list <- list()
years <- 2013:2017
for (idx in 1:length(years)) {
	year <- years[idx]
	obj <- get(sprintf("acs5_%d", year))
	g <- ggplot(obj, aes(x=DirectEst, y=E.mean)) +
		geom_point(size = 2) +
		geom_abline(intercept = 0, slope = 1, color="red",
			linetype="dashed", size=1.2) +
		ggtitle(sprintf("%d 5yr ACS Direct Estimates", year)) +
		labs(x = "Direct Estimate", y = "Model-Based Estimate") +
		theme_bw()
	scatter.list[[idx]] <- g
}
marrangeGrob(scatter.list, nrow = 3, ncol = 2)

idx.missing2017 <- which(is.na(acs5_2017$DirectEst))

#' Plot neighborhood areas (target supports) among ACS 5-year direct estimates;
#' this gives a sense of whether the model-based esimtates are reasonable.
#' This map takes a bit of preparation.

Central <- neighbs[1,]
East <- neighbs[2,]
North <- neighbs[3,]
Paris <- neighbs[4,]
Missing1 <- acs5_2017[idx.missing2017[1],]
Missing2 <- acs5_2017[idx.missing2017[2],]
Missing3 <- acs5_2017[idx.missing2017[3],]
Missing4 <- acs5_2017[idx.missing2017[4],]

# Prevent `sf` package warnings like "st_centroid assumes attributes are
# constant over geometries of x"
st_agr(Central) <- "constant"
st_agr(East) <- "constant"
st_agr(North) <- "constant"
st_agr(Paris) <- "constant"
st_agr(Missing1) <- "constant"
st_agr(Missing2) <- "constant"
st_agr(Missing3) <- "constant"
st_agr(Missing4) <- "constant"

Central.coord <- st_coordinates(st_centroid(Central))
East.coord <- st_coordinates(st_centroid(East))
North.coord <- st_coordinates(st_centroid(North))
Paris.coord <- st_coordinates(st_centroid(Paris))
Missing1.coord <- st_coordinates(st_centroid(Missing1))
Missing2.coord <- st_coordinates(st_centroid(Missing2))
Missing3.coord <- st_coordinates(st_centroid(Missing3))
Missing4.coord <- st_coordinates(st_centroid(Missing4))

g <- ggplot(acs5_2017) +
	geom_sf(colour = "black", size = 0.05, aes(fill = DirectEst)) +
	ggtitle("Median Household Income for Boone County",
		subtitle = "ACS 2017 5yr Direct Estimates") +
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
