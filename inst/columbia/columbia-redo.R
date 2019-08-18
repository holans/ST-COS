#' ---
#' title: Analysis of City of Columbia Neighborhoods
#' author: 
#' date: "Last Updated: `r format(Sys.time(), '%B %d, %Y')`"
#' output:
#'   html_document:
#'	 number_sections: true
#' ---

#' Let's try to redo the Columbia example without the STCOSPrep class. This may
#' require changes to the package. I'm thinking it will result in a more modular
#' package that can be reused for other types of analysis than fitting the model
#' in the 2015 Stat paper...

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
dom_fine = block_groups(state = '29', county = '019', year = 2017) %>%
	st_as_sf() %>%
	st_transform(crs = 3857) %>%
	mutate(geoid = GEOID) %>%
	select(-GEOID)
n = nrow(dom_fine)

#' A quick plot of the fine-level domain.
ggplot(dom_fine) +
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
neighbs = columbia_neighbs %>%
	st_transform(crs = st_crs(dom_fine))

#' A quick plot of the neighborhoods.
ggplot(neighbs) +
	geom_sf(colour = "black", size = 0.05) +
	theme_bw()

#' Gather the data for source supports
#' The elements in `period_list` correspond to `source_list`, so that
#' `period_list[[l]]` describes the lookback period for `source_list[[l]]`.
source_list = list(
	acs5_2013 %>% filter(!is.na(DirectEst)) %>% filter(!is.na(DirectVar)),
	acs5_2014 %>% filter(!is.na(DirectEst)) %>% filter(!is.na(DirectVar)),
	acs5_2015 %>% filter(!is.na(DirectEst)) %>% filter(!is.na(DirectVar)),
	acs5_2016 %>% filter(!is.na(DirectEst)) %>% filter(!is.na(DirectVar)),
	acs5_2017 %>% filter(!is.na(DirectEst)) %>% filter(!is.na(DirectVar))
)
period_list = list(2009:2013, 2010:2014, 2011:2015, 2012:2016, 2013:2017)
times_seq = 2009:2017
L = length(source_list)
T = length(times_seq)

# TBD: No non-NA observation intersects with some particular areas
# I think we need to drop them from the analysis
# Does this help with the MLE issue?
U = rbind(
	overlap_matrix(source_list[[1]], dom_fine, proportion = FALSE),
	overlap_matrix(source_list[[2]], dom_fine, proportion = FALSE),
	overlap_matrix(source_list[[3]], dom_fine, proportion = FALSE),
	overlap_matrix(source_list[[4]], dom_fine, proportion = FALSE),
	overlap_matrix(source_list[[5]], dom_fine, proportion = FALSE)
)
idx = which(colSums(U) < 10)
dom_fine = dom_fine[-idx,]
n = nrow(dom_fine)

#' # Prepare to fit the model
#+ message=FALSE
library(fields)
library(stcos)
library(ggforce)

#' Select spatial knots via space-filling design. Start with uniformly
#' drawn points from the fine-level geography and use `cover.design`
#' to select a subset of those points for knots.
u = st_sample(dom_fine, size = 2000)
M = matrix(unlist(u), length(u), 2, byrow = TRUE)
# out = cover.design(M, 500)
out = cover.design(M, 200)
knots_sp = out$design

#' Select temporal knots to be evenly spaced over the years relevant
#' to the source support years.
knots_t = seq(2009, 2017, by = 0.5)

#' Use Cartesian join of spatial and temporal knots to obtain spatio-temporal knots.
knots = merge(knots_sp, knots_t)

#' Create an ArealSpaceTimeBisquareBasis object with our knot points
basis1 = ArealSpaceTimeBisquareBasis$new(knots[,1], knots[,2], knots[,3],
	w_s = 1, w_t = 1, mc_reps = 200)

#' Here is a plot of the spatial knots. We plot a circle around one of the
#' points to illustrate the choice of the `w.s` argument.
rad = basis1$get_basis_spt()$get_rl()
knots_sp_dat = data.frame(x = knots_sp[,1], y = knots_sp[,2], r = rad)
g = ggplot(dom_fine) +
	geom_sf(colour = "black", size = 0.25, fill = NA) +
	geom_point(data = knots_sp_dat, aes(x, y), lwd = 1, col = "red") +
	geom_point(data = knots_sp_dat[1,], aes(x, y), lwd = 3, col = "blue") +
	geom_circle(data = knots_sp_dat[1,], aes(x0=x, y0=y, r=r), fill = NA,
		lwd = 0.5, col = "blue") +
	labs(x = "", y = "") +
	theme_bw()
print(g)

#' Build components for STCOS model

# Compute overlap matrix H
H = Matrix(0, 0, n)
for (l in 1:L) {
	logger("Computing overlap matrix for domain %d of %d\n", l, L)
	H_l = overlap_matrix(source_list[[l]], dom_fine)
	H = rbind(H, H_l)
}
N = nrow(H)

# Compute basis function matrix S
S_full = Matrix(0, 0, nrow(knots))
for (l in 1:L) {
	logger("Computing basis matrix for domain %d of %d\n", l, L)
	S_l = basis1$compute(source_list[[l]], period_list[[l]])
	S_full = rbind(S_full, S_l)
}

# Extract the direct estimates and variance estimates
z = unlist(Map(function(x) { x$DirectEst }, source_list))
v = unlist(Map(function(x) { x$DirectVar }, source_list))

#' Do a PCA reduction on S to reduce its dimension
# TBD: Do we want to use package to return limited number of eigencomponents?
eig = eigen(t(S_full) %*% S_full)
idx.S = which(cumsum(eig$values) / sum(eig$values) < 0.65)
Ts = eig$vectors[,idx.S]
S = S_full %*% Ts
r = ncol(S)

#' Plot the proportion of variation captured by our selection of PCA components.
eigprops = cumsum(eig$values) / sum(eig$values)
plot(eigprops[1:200], xlab = "Dimension", ylab = "Proportion of Variation")
abline(v = max(idx.S), lty = 2)
abline(h = eigprops[max(idx.S)], lty = 2)

#' Compute basis function on fine-level domain, which is required for some
#' structures of K.
S_fine_full = Matrix(0, 0, nrow(knots))
for (l in 1:length(times_seq)) {
	logger("Computing basis matrix for fine-level domain, time %d\n", times_seq[l])
	S_l = basis1$compute(dom_fine, times_seq[l])
	S_fine_full = rbind(S_fine_full, S_l)
}
S_fine = S_fine_full %*% Ts

# Compute adjacency matrix for fine-level support
A = adjacency_matrix(dom_fine)
W = 1/rowSums(A) * A
tau = 0.9
Q = Diagonal(n,1) - tau*W
Qinv = solve(Q)

#' Pick a covariance structure for random coefficients of basis expansion.
method = "independence"

if (method == "moran") {
	# Assume covariance structure with M computed via Moran's I basis
	# This method requires an X matrix. Create an X matrix using spatial-only basis.
	# We need to reduce its dimension, just as we did with S
	basis2 = ArealSpatialBisquareBasis$new(knots_sp[,1], knots_sp[,2], w = 1, mc_reps = 200)
	X_full = basis2$compute(dom_fine)
	eig2 = eigen(t(X_full) %*% X_full)
	cumsum(eig2$values) / sum(eig2$values)
	X = X_full %*% eig2$vectors[,1:10]
	K = cov_approx_moran(Qinv, X, S_fine, lag_max = T)
} else if (method == "randomwalk") {
	# Random Walk
	K = cov_approx_randwalk(Qinv, S_fine, lag_max = T)
} else if (method == "car") {
	# Spatial-only (CAR)
	# Assume covariance structure without dependence over time
	K = cov_approx_blockdiag(Qinv, S_fine, lag_max = T)
} else if (method == "independence") {
	# Independence
	K = Diagonal(n = r)
}
K_inv = solve(K)

#' Standardize observations before running MCMC.
z_mean = mean(z)
z_sd = sd(z)
z_scaled = (z - z_mean) / z_sd
v_scaled = v / z_sd^2

#' # Fit the Model

#' Fit the model with Stan
library(rstan)
stan_dat = list(
	N = N, n = n, r = r, z = z_scaled, v = v_scaled, H = as.matrix(H),
	S = as.matrix(S), K = as.matrix(K),
	alpha_K = 1, beta_K = 2, alpha_xi = 1, beta_xi = 2, alpha_mu = 1, beta_mu = 2
)
fit = stan(file = "stcos.stan", data = stan_dat, iter = 2000, chains = 1,
	verbose = TRUE, sample_file = "stan_draws.csv")
print(fit, par = c("eta", "sig2mu", "sig2xi", "sig2K"))

traceplot(fit, par = c("sig2mu", "sig2xi", "sig2K"), nrow = 3, ncol = 1)
traceplot(fit, par = sprintf("mu[%d]", 1:6))
traceplot(fit, par = c("eta"))
traceplot(fit, par = c("lp__"))

stan_hist(fit, par = c("sig2mu", "sig2xi", "sig2K"), bins = 30, nrow = 3, ncol = 1)
stan_hist(fit, par = sprintf("mu[%d]", 1:6), bins = 30)
stan_hist(fit, par = c("eta"), bins = 30)
stan_hist(fit, par = c("lp__"), bins = 30)

#+ message=FALSE
library(coda)

#' Fit MLE; this will serve as an initial value for MCMC.
mle_out = mle_stcos(z = z_scaled, v = v_scaled, H = H, S = S, K = K,
	init = list(sig2K = 1, sig2xi = 1))
init = list(
	sig2K = mle_out$sig2K_hat,
	sig2xi = mle_out$sig2xi_hat,
	muB = mle_out$mu_hat
)
hyper = list(a_sig2K = 1, b_sig2K = 2, a_sig2xi = 1, b_sig2xi = 2,
	a_sig2mu = 1, b_sig2mu = 2)

#' Run the Gibbs sampler.
gibbs_out = gibbs_stcos_raw(z = z_scaled, v = v_scaled, H = H, S = S,
	K_inv = K_inv, R = 10000, report_period = 2000, burn = 2000,
	thin = 10, init = init, hyper = hyper)
print(gibbs_out)

#' Show some trace plots to assess convergence of the sampler.
plot((muB_mcmc = mcmc(gibbs_out$muB_hist))[,1:3])
plot((xi_mcmc = mcmc(gibbs_out$xi_hist))[,1:3])
plot((eta_mcmc = mcmc(gibbs_out$eta_hist))[,1:3])

varcomps_mcmc = mcmc(cbind(
	gibbs_out$sig2mu_hist,
	gibbs_out$sig2xi_hist,
	gibbs_out$sig2K_hist
))
colnames(varcomps_mcmc) = c("sig2mu", "sig2xi", "sig2K")
plot(varcomps_mcmc)

#' #  Produce Results on target supports
#' Compute `H` and `S` matrices and get summaries of posterior distribution for E(Y).
#' Use 90% significance for all credible intervals and MOEs.
append_results = function(dat_sf, period, alpha = 0.10) {
	H_new = overlap_matrix(dat_sf, dom_fine)
	S_new_full = basis1$compute(dat_sf, period)
	S_new = S_new_full %*% Ts

	E_hat_scaled = fitted(gibbs_out, H_new, S_new)
	E_hat = z_sd * E_hat_scaled + z_mean					  # Uncenter and unscale
	dat_sf$E_mean = colMeans(E_hat)						      # Point estimates
	dat_sf$E_sd = apply(E_hat, 2, sd)						  # SDs
	dat_sf$E_lo = apply(E_hat, 2, quantile, prob = alpha/2)   # Credible interval lo
	dat_sf$E_hi = apply(E_hat, 2, quantile, prob = 1-alpha/2) # Credible interval hi
	dat_sf$E_median = apply(E_hat, 2, median)				  # Median
	dat_sf$E_moe = apply(E_hat, 2, sd) * qnorm(1-alpha/2)	  # MOE
	return(dat_sf)
}

acs5_2013 = append_results(acs5_2013, period = 2009:2013)
acs5_2014 = append_results(acs5_2014, period = 2010:2014)
acs5_2015 = append_results(acs5_2015, period = 2011:2015)
acs5_2016 = append_results(acs5_2016, period = 2012:2016)
acs5_2017 = append_results(acs5_2017, period = 2013:2017)
neighbs = append_results(neighbs, period = 2013:2017)

#' The objective of our analysis - predictions on the four target neighborhoods.
print(neighbs)

#' # Plot Results
#+ message=FALSE
library(gridExtra)
library(ggrepel)

#' Maps of direct and model-based 2017 5-year estimates.
lim_est = range(acs5_2017$DirectEst, acs5_2017$E_mean)
g = ggplot(acs5_2017) +
	geom_sf(colour = "black", size = 0.05, aes(fill = DirectEst)) +
	ggtitle("Median Household Income\nfor Boone County",
		subtitle = "2017 5yr ACS Direct Estimates") +
	scale_fill_distiller("DirectEst", palette = "RdYlBu", limits = lim_est) +
	theme_bw()
h = ggplot(acs5_2017) +
	geom_sf(colour = "black", size = 0.05, aes(fill = E_mean)) +
	ggtitle("Median Household Income\nfor Boone County",
		subtitle = "2017 5yr Model Estimates") +
	scale_fill_distiller("E_mean", palette = "RdYlBu", limits = lim_est) +
	theme_bw()
k = grid.arrange(g, h, ncol = 2)

#' Scatter plots comparing direct and model-based 5-year estimates for
#' 2013, ..., 2017.
scatter_list = list()
years = 2013:2017
for (idx in 1:length(years)) {
	year = years[idx]
	obj = get(sprintf("acs5_%d", year))
	g = ggplot(obj, aes(x=DirectEst, y=E_mean)) +
		geom_point(size = 2) +
		geom_abline(intercept = 0, slope = 1, color="red",
			linetype="dashed", size=1.2) +
		ggtitle(sprintf("%d 5yr ACS Direct Estimates", year)) +
		labs(x = "Direct Estimate", y = "Model-Based Estimate") +
		theme_bw()
	scatter_list[[idx]] = g
}
marrangeGrob(scatter_list, nrow = 3, ncol = 2)

idx_missing2017 = which(is.na(acs5_2017$DirectEst))

#' Plot neighborhood areas (target supports) among ACS 5-year direct estimates;
#' this gives a sense of whether the model-based esimtates are reasonable.
#' This map takes a bit of preparation.

Central = neighbs[1,]
East = neighbs[2,]
North = neighbs[3,]
Paris = neighbs[4,]
Missing1 = acs5_2017[idx_missing2017[1],]
Missing2 = acs5_2017[idx_missing2017[2],]
Missing3 = acs5_2017[idx_missing2017[3],]
Missing4 = acs5_2017[idx_missing2017[4],]

# Prevent `sf` package warnings like "st_centroid assumes attributes are
# constant over geometries of x"
st_agr(Central) = "constant"
st_agr(East) = "constant"
st_agr(North) = "constant"
st_agr(Paris) = "constant"
st_agr(Missing1) = "constant"
st_agr(Missing2) = "constant"
st_agr(Missing3) = "constant"
st_agr(Missing4) = "constant"

Central_coord = st_coordinates(st_centroid(Central))
East_coord = st_coordinates(st_centroid(East))
North_coord = st_coordinates(st_centroid(North))
Paris_coord = st_coordinates(st_centroid(Paris))
Missing1_coord = st_coordinates(st_centroid(Missing1))
Missing2_coord = st_coordinates(st_centroid(Missing2))
Missing3_coord = st_coordinates(st_centroid(Missing3))
Missing4_coord = st_coordinates(st_centroid(Missing4))

g = ggplot(acs5_2017) +
	geom_sf(colour = "black", size = 0.05, aes(fill = DirectEst)) +
	ggtitle("Median Household Income for Boone County",
		subtitle = "ACS 2017 5yr Direct Estimates") +
	scale_fill_distiller("DirectEst", palette = "RdYlBu", limits = lim_est) +
	geom_sf(data = neighbs, fill = "black") +
	geom_label_repel(data = st_centroid(East), nudge_x = 20000, nudge_y = 130000,
		aes(x=East_coord[1], y=East_coord[2], label="East")) +
	geom_label_repel(data = st_centroid(Central), nudge_x = -10000, nudge_y = 130000,
		aes(x=Central_coord[1], y=Central_coord[2], label="Central")) +
	geom_label_repel(data = st_centroid(North), nudge_x = 0, nudge_y = 100000,
		aes(x=North_coord[1], y=North_coord[2], label="North")) +
	geom_label_repel(data = st_centroid(Paris), nudge_x = 10000, nudge_y = 100000,
		aes(x=Paris_coord[1], y=Paris_coord[2], label="Paris")) +
	geom_label_repel(data = st_centroid(Missing1), nudge_x = 100000, nudge_y = -50000,
		aes(x=Missing1_coord[1], y=Missing1_coord[2], label="Missing1")) +
	geom_label_repel(data = st_centroid(Missing2), nudge_x = 100000, nudge_y = -36000,
		aes(x=Missing2_coord[1], y=Missing2_coord[2], label="Missing2")) +
	geom_label_repel(data = st_centroid(Missing3), nudge_x = -100000, nudge_y = -50000,
		aes(x=Missing3_coord[1], y=Missing3_coord[2], label="Missing3")) +
	geom_label_repel(data = st_centroid(Missing4), nudge_x = -100000, nudge_y = -38000,
		aes(x=Missing4_coord[1], y=Missing4_coord[2], label="Missing4")) +
	xlab(NULL) +
	ylab(NULL) +
	theme_bw()
print(g)
