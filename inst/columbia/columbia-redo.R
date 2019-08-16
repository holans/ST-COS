#' ---
#' title: Analysis of City of Columbia Neighborhoods
#' author: 
#' date: "Last Updated: `r format(Sys.time(), '%B %d, %Y')`"
#' output:
#'   html_document:
#'     number_sections: true
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
dom.fine = block_groups(state = '29', county = '019', year = 2017) %>%
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
neighbs = columbia_neighbs %>%
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
u = st_sample(dom.fine, size = 2000)
M = matrix(unlist(u), length(u), 2, byrow = TRUE)
out = cover.design(M, 500)
knots.sp = out$design

#' Select temporal knots to be evenly spaced over the years relevant
#' to the source support years.
knots.t = seq(2009, 2017, by = 0.5)

#' Use Cartesian join of spatial and temporal knots to obtain spatio-temporal knots.
knots = merge(knots.sp, knots.t)
basis = SpaceTimeBisquareBasis$new(knots[,1], knots[,2], knots[,3], w.s = 1, w.t = 1)

#' Here is a plot of the spatial knots. We plot a circle around one of the
#' points to illustrate the choice of the `w.s` argument.
knots.sp.dat = data.frame(x = knots.sp[,1], y = knots.sp[,2], r = basis$get_rl())
g = ggplot(dom.fine) +
	geom_sf(colour = "black", size = 0.25, fill = NA) +
	geom_point(data = knots.sp.dat, aes(x, y), lwd = 1, col = "red") +
	geom_point(data = knots.sp.dat[1,], aes(x, y), lwd = 3, col = "blue") +
	geom_circle(data = knots.sp.dat[1,], aes(x0=x, y0=y, r=r), fill = NA,
		lwd = 0.5, col = "blue") +
	labs(x = "", y = "") +
	theme_bw()
print(g)

#' Build components for STCOS model
area_list = list(acs5_2013, acs5_2014, acs5_2015, acs5_2016, acs5_2017)
period_list = list(2009:2013, 2010:2014, 2011:2015, 2012:2016, 2013:2017)
times_all = 2009:2017
L = length(area_list)

# Compute overlap matrix H
# TBD: Should we normalize in compute.overlap?
# TBD: Should we rename function to overlap or overlap_matrix or cos_matrix?
# TBD: maybe we should set this up so we don't need the awkward transpose?
# TBD: Can we avoid even having this function at all?
# set_agr("constant") suppresses unhelpful warnings when calling st_intersection
H = Matrix(0, 0, n)
for (l in 1:L) {
	logger("Computing overlap matrix for domain %d\n", l)
	H_l = intersection_matrix(area_list[[l]], dom.fine)
	
	# TBD: Make sure this is normalized and transposed correctly...
	H = rbind(H, apply(H_l, 1, normalize))
}

# Compute basis function matrix S
# TBD: Maybe areal basis functions should be classes? Does it make sense?
# Or over-engineered?
S = Matrix(0, 0, nrow(knots))
for (l in 1:L) {
	logger("Computing basis matrix for domain %d\n", l)
	S_l = compute_spt_basis_mc(basis = basis, domain = area_list[[l]],
		R = 50, period = period_list[[l]], report.period = 100)
	S = rbind(S, S_l)
}

#' Do a PCA reduction on S to reduce its dimension
# TBD: Do we want to use package to return limited number of eigencomponents?
eig = eigen(t(S) %*% S)
idx.S = which(cumsum(eig$values) / sum(eig$values) < 0.65)
Tx = eig$vectors[,idx.S]
S_reduced = S %*% Tx

#' Plot the proportion of variation captured by our selection of PCA components.
eigprops = cumsum(eig$values) / sum(eig$values)
plot(eigprops[1:200], xlab = "Dimension", ylab = "Proportion of Variation")
abline(v = max(idx.S), lty = 2)
abline(h = eigprops[max(idx.S)], lty = 2)

# Record the dimensions for model components
N = nrow(H)
n = nrow(dom.fine)
r = ncol(S_reduced)
T = length(times_all)

# Extract the direct estimates and variance estimates
z = numeric(0)
v = numeric(0)
for (l in 1:L) {
	z = c(z, area_list[[l]]$DirectEst)
	v = c(v, area_list[[l]]$DirectVar)
}

idx_missing = which(is.na(z))
idx_nonmissing = which(!is.na(z))

S_fine = Matrix(0, 0, nrow(knots))
for (l in 1:length(times_all)) {
	logger("Computing basis matrix for fine-level domain, time %d\n", times_all[l])
	S_l = compute_spt_basis_mc(basis = basis, domain = dom.fine,
		R = 50, period = times_all[l], report.period = 100)
	S_fine = rbind(S_fine, S_l)
}
S_fine_reduced = S_fine %*% Tx

# Compute adjacency matrix for fine-level support
out = st_touches(dom.fine, dom.fine)
A = adjList2Matrix(out)
countAdj = Matrix(0, nrow(A), ncol(A))
s = rowSums(A)
for (j in 1:n) {
	if (s[j] > 0) {
		countAdj[j,] = A[j,] / s[j]
	}
}
Q = Diagonal(n,1) - 0.9*countAdj
Qinv = solve(Q)

#' Pick a covariance structure for random coefficients of basis expansion.
method = "moran"

if (method == "moran") {
	# Assume covariance structure with M computed via Moran's I basis
	K_inv = sp$get_Kinv(2009:2017, method = "moran")

	# Create an X matrix using spatial-only basis
	basis_sp = SpatialBisquareBasis$new(knots.sp[,1], knots.sp[,2], w = 1)
	X_full = compute_sp_basis_mc(basis = basis_sp, domain = dom.fine,
		R = 500, report.period = 100)
	eig = eigen(t(X_full) %*% X_full)
	cumsum(eig$values) / sum(eig$values)
	X = X_full %*% eig$vectors[,1:10]

	P_perp = Diagonal(nrow(X),1) - X %*% solve(t(X) %*% X, t(X))
	eig = eigen(P_perp, symmetric = TRUE)
	M = Re(eig$vectors)
	M = (M + t(M)) / 2
	Sigma = sptcovar.vectautoreg(Qinv, M, S_fine_reduced, lag_max = T)

	# TBD: I don't think this is right... the sptcovar functions are already returning the minimizer K ...
	# The "covariance_approximant" function I have now is doing some kind of pseudo-inverse ...
	K_inv = covariance_approximant(Sigma, S_fine_reduced)
} else if (method == "randomwalk") {
	# Random Walk
	# Assume covariance structure with M as identity matrix
	M = Diagonal(n,1)
	Sigma = sptcovar.randwalk(Qinv, M, S_fine_reduced, lag_max = T)

	# TBD: Is this right??
	K_inv = covariance_approximant(Sigma, S_fine_reduced)
} else if (method == "car") {
	# Spatial-only (CAR)
	# Assume covariance structure without dependence over time
	Sigma = sptcovar.indep(Qinv, S_fine_reduced, lag_max = T)
	
	# TBD: Is this right??
	K_inv = covariance_approximant(Sigma, S_fine_reduced)
} else if (method == "independence") {
	# Independence
	K_inv = Diagonal(n = N)
}

#' Standardize observations before running MCMC.
z_mean = mean(z, na.rm = TRUE)
z_sd = sd(z, na.rm = TRUE)
z_scaled = (z - z_mean) / z_sd
v_scaled = v / z_sd^2

#' # Fit the Model
#+ message=FALSE
library(coda)

#' Fit MLE; this will serve as an initial value for MCMC.
K = solve(K_inv)
mle.out = mle.stcos(
	z = z_scaled[idx_nonmissing],
	v = v_scaled[idx_nonmissing],
	H = H[idx_nonmissing,],
	S = S_reduced[idx_nonmissing,],
	K = K,
	init = list(sig2K = 1, sig2xi = 1)
)
init = list(
	sig2K = mle.out$sig2K.hat,
    sig2xi = mle.out$sig2xi.hat,
    mu_B = mle.out$mu.hat
)

#' Run the Gibbs sampler.
gibbs.out = gibbs.stcos.raw(
	z = z.scaled[idx_nonmissing],
	v = v.scaled[idx_nonmissing],
	H = H[idx_nonmissing,],
	S = S_reduced[idx_nonmissing,],
	K.inv = K_inv,
	R = 10000, report.period = 2000, burn = 2000, thin = 10, init = init)
print(gibbs.out)

#' Show some trace plots to assess convergence of the sampler.
plot((mu_B.mcmc = mcmc(gibbs.out$mu_B.hist))[,1:3])
plot((xi.mcmc = mcmc(gibbs.out$xi.hist))[,1:3])
plot((eta.mcmc = mcmc(gibbs.out$eta.hist))[,1:3])

varcomps.mcmc = mcmc(cbind(
	gibbs.out$sig2mu.hist,
	gibbs.out$sig2xi.hist,
	gibbs.out$sig2K.hist
))
colnames(varcomps.mcmc) = c("sig2mu", "sig2xi", "sig2K")
plot(varcomps.mcmc)

#' #  Produce Results on target supports
#' Compute `H` and `S` matrices and get summaries of posterior distribution for E(Y).
#' Use 90% significance for all credible intervals and MOEs.
append_results = function(dat_sf, period, alpha = 0.10) {
	H_raw = intersection_matrix(dat_sf, dom.fine)
	H_new = Matrix(apply(H_raw, 1, normalize))

	S_new = compute_spt_basis_mc(basis = basis, domain = dat_sf,
		R = 50, period = period, report.period = 100)
	S_new_reduced = S_new %*% Tx

	E.hat.scaled = fitted(gibbs.out, H_new, S_new_reduced)
	E.hat = z.sd * E.hat.scaled + z.mean                      # Uncenter and unscale
	dat_sf$E.mean = colMeans(E.hat)                           # Point estimates
	dat_sf$E.sd = apply(E.hat, 2, sd)                         # SDs
	dat_sf$E.lo = apply(E.hat, 2, quantile, prob = alpha/2)   # Credible interval lo
	dat_sf$E.hi = apply(E.hat, 2, quantile, prob = 1-alpha/2) # Credible interval hi
	dat_sf$E.median = apply(E.hat, 2, median)                 # Median
	dat_sf$E.moe = apply(E.hat, 2, sd) * qnorm(1-alpha/2)     # MOE
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
lim.est = range(acs5_2017$DirectEst, acs5_2017$E.mean)
g = ggplot(acs5_2017) +
	geom_sf(colour = "black", size = 0.05, aes(fill = DirectEst)) +
	ggtitle("Median Household Income\nfor Boone County",
		subtitle = "2017 5yr ACS Direct Estimates") +
	scale_fill_distiller("DirectEst", palette = "RdYlBu", limits = lim.est) +
	theme_bw()
h = ggplot(acs5_2017) +
	geom_sf(colour = "black", size = 0.05, aes(fill = E.mean)) +
	ggtitle("Median Household Income\nfor Boone County",
		subtitle = "2017 5yr Model Estimates") +
	scale_fill_distiller("E.mean", palette = "RdYlBu", limits = lim.est) +
	theme_bw()
k = grid.arrange(g,h, ncol = 2)

#' Scatter plots comparing direct and model-based 5-year estimates for
#' 2013, ..., 2017.
scatter.list = list()
years = 2013:2017
for (idx in 1:length(years)) {
	year = years[idx]
	obj = get(sprintf("acs5_%d", year))
	g = ggplot(obj, aes(x=DirectEst, y=E.mean)) +
		geom_point(size = 2) +
		geom_abline(intercept = 0, slope = 1, color="red",
			linetype="dashed", size=1.2) +
		ggtitle(sprintf("%d 5yr ACS Direct Estimates", year)) +
		labs(x = "Direct Estimate", y = "Model-Based Estimate") +
		theme_bw()
	scatter.list[[idx]] = g
}
marrangeGrob(scatter.list, nrow = 3, ncol = 2)

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

Central.coord = st_coordinates(st_centroid(Central))
East.coord = st_coordinates(st_centroid(East))
North.coord = st_coordinates(st_centroid(North))
Paris.coord = st_coordinates(st_centroid(Paris))
Missing1.coord = st_coordinates(st_centroid(Missing1))
Missing2.coord = st_coordinates(st_centroid(Missing2))
Missing3.coord = st_coordinates(st_centroid(Missing3))
Missing4.coord = st_coordinates(st_centroid(Missing4))

g = ggplot(acs5_2017) +
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
