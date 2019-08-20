library(tigris)
library(ggplot2)
library(dplyr)
library(fields)
library(rstan)
library(stcos)

set.seed(1234)

# ----- Load data for analysis -----
# Load the fine-level support via shapefile using the `tigris` package.
# We use the 2017 block-groups in Boone County, MO. Convert to an `sf`
# object, and transform to the projection with EPSG code 3857.
options(tigris_use_cache = TRUE)
options(tigris_refresh = FALSE)
dom_fine = block_groups(state = '29', county = '019', year = 2017) %>%
	st_as_sf() %>%
	st_transform(crs = 3857) %>%
	mutate(geoid = GEOID) %>%
	select(-GEOID)
n = nrow(dom_fine)

# Load previously constructed source supports.
data(acs_sf)

# Load Target Supports
# Make sure to transform to the same projection as the fine-level support.
data(columbia_neighbs)
neighbs = columbia_neighbs %>% st_transform(crs = st_crs(dom_fine))

# Gather the data for source supports
# Remove observations with NA for the direct estimate or variance.
# The elements in `period_list` correspond to `source_list`, so that
# `period_list[[l]]` describes the lookback period for `source_list[[l]]`.
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
N = sum(unlist(Map(nrow, source_list)))

# ----- Remove unused areas in fine-level domain -----
# For some areas in the fine-level domain, there is very little overlap
# area with source supports (at least the ones having non-NA values
# which will be used to fit the model). We will drop these areas
# from the fine-level domain to avoid rank-deficiency issues with the
# H matrix later.
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

# ----- Assemble design and variance matrices needed to fit model -----
# Select spatial knots via space-filling design. Start with uniformly
# drawn points from the fine-level geography and use `cover.design`
# to select a subset of those points for knots.
u = st_sample(dom_fine, size = 2000)
M = matrix(unlist(u), length(u), 2, byrow = TRUE)
out = cover.design(M, 200)
knots_sp = out$design

# Select temporal knots to be evenly spaced over the years relevant
# to the source support years.
knots_t = seq(2009, 2017, by = 0.5)

# Use Cartesian join of spatial and temporal knots to obtain spatio-temporal knots.
knots = merge(knots_sp, knots_t)
r_full = nrow(knots)

# Create an ArealSpaceTimeBisquareBasis object with our knot points
bs_spt = ArealSpaceTimeBisquareBasis$new(knots[,1], knots[,2], knots[,3],
	w_s = 1, w_t = 1, mc_reps = 200)

# Compute overlap matrix H
H = Matrix(0, 0, n)
for (l in 1:L) {
	logger("Computing overlap matrix for domain %d of %d\n", l, L)
	H_l = overlap_matrix(source_list[[l]], dom_fine)
	H = rbind(H, H_l)
}

# Compute basis function matrix S
S_full = Matrix(0, 0, r_full)
for (l in 1:L) {
	logger("Computing basis matrix for domain %d of %d\n", l, L)
	S_l = bs_spt$compute(source_list[[l]], period_list[[l]])
	S_full = rbind(S_full, S_l)
}

# Extract the direct estimates and variance estimates
z = unlist(Map(function(x) { x$DirectEst }, source_list))
v = unlist(Map(function(x) { x$DirectVar }, source_list))

# Reduce the dimension of `S_full` using PCA
eig = eigen(t(S_full) %*% S_full)
idx_S = which(cumsum(eig$values) / sum(eig$values) < 0.65)
Ts = eig$vectors[,idx_S]
S = S_full %*% Ts
r = ncol(S)

# Compute basis function on fine-level domain, which is required for some
# structures of K.
S_fine_full = Matrix(0, 0, r_full)
for (l in 1:length(times_seq)) {
	logger("Computing basis matrix for fine-level domain, time %d\n", times_seq[l])
	S_l = bs_spt$compute(dom_fine, times_seq[l])
	S_fine_full = rbind(S_fine_full, S_l)
}
S_fine = S_fine_full %*% Ts

# Compute adjacency matrix for fine-level support
A = adjacency_matrix(dom_fine)
aa = rowSums(A) + (rowSums(A) == 0)
W = 1/aa * A
tau = 0.9
Q = Diagonal(n,1) - tau*W
Qinv = solve(Q)

# Construct K with a "Random Walk" covariance structure.
K = cov_approx_randwalk(Qinv, S_fine, lag_max = T)
K_inv = solve(K)

#' Standardize observations before MCMC.
z_mean = mean(z)
z_sd = sd(z)
z_scaled = (z - z_mean) / z_sd
v_scaled = v / z_sd^2

# ----- Fit the model with Stan and produce outputs -----
stan_dat = list(
	N = N, n = n, r = r, z = z_scaled, v = v_scaled, H = as.matrix(H),
	S = as.matrix(S), K = as.matrix(K),
	alpha_K = 1, beta_K = 2, alpha_xi = 1, beta_xi = 2, alpha_mu = 1, beta_mu = 2
)
fit = stan(file = "stcos.stan", data = stan_dat, iter = 2000, chains = 1,
	verbose = TRUE)
stan_out = extract(fit, pars = c("mu", "eta"))

append_results = function(dat_sf, period, alpha = 0.10) {
	H_new = overlap_matrix(dat_sf, dom_fine)
	S_new_full = bs_spt$compute(dat_sf, period)
	S_new = S_new_full %*% Ts

	# Get draws of the mean E(Y), then uncenter and unscale
	EY_scaled = stan_out$mu %*% t(H_new) + stan_out$eta %*% t(S_new)
	A = z_sd * EY_scaled + z_mean

	dat_sf$E_mean = colMeans(A)                           # Point estimates
	dat_sf$E_sd = apply(A, 2, sd)                         # SDs
	dat_sf$E_lo = apply(A, 2, quantile, prob = alpha/2)   # Credible interval lo
	dat_sf$E_hi = apply(A, 2, quantile, prob = 1-alpha/2) # Credible interval hi
	dat_sf$E_median = apply(A, 2, median)                 # Median
	dat_sf$E_moe = apply(A, 2, sd) * qnorm(1-alpha/2)     # MOE
	return(dat_sf)
}

acs5_2013 = append_results(acs5_2013, period = 2009:2013)
acs5_2014 = append_results(acs5_2014, period = 2010:2014)
acs5_2015 = append_results(acs5_2015, period = 2011:2015)
acs5_2016 = append_results(acs5_2016, period = 2012:2016)
acs5_2017 = append_results(acs5_2017, period = 2013:2017)
neighbs = append_results(neighbs, period = 2013:2017)

# The objective of our analysis - predictions on the four target neighborhoods.
print(neighbs)

# ----- Run the Gibbs sampler and produce outputs -----
hyper = list(a_sig2K = 1, b_sig2K = 2, a_sig2xi = 1, b_sig2xi = 2,
	a_sig2mu = 1, b_sig2mu = 2)
gibbs_out = gibbs_stcos(z = z_scaled, v = v_scaled, H = H, S = S,
	K_inv = K_inv, R = 10000, report_period = 2000, burn = 2000,
	thin = 10, hyper = hyper)
print(gibbs_out)

append_results = function(dat_sf, period, alpha = 0.10) {
	H_new = overlap_matrix(dat_sf, dom_fine)
	S_new_full = bs_spt$compute(dat_sf, period)
	S_new = S_new_full %*% Ts

	# Get draws of the mean E(Y), then uncenter and unscale
	EY_scaled = fitted(gibbs_out, H_new, S_new)
	A = z_sd * EY_scaled + z_mean

	dat_sf$E_mean = colMeans(A)                           # Point estimates
	dat_sf$E_sd = apply(A, 2, sd)                         # SDs
	dat_sf$E_lo = apply(A, 2, quantile, prob = alpha/2)   # Credible interval lo
	dat_sf$E_hi = apply(A, 2, quantile, prob = 1-alpha/2) # Credible interval hi
	dat_sf$E_median = apply(A, 2, median)                 # Median
	dat_sf$E_moe = apply(A, 2, sd) * qnorm(1-alpha/2)     # MOE
	return(dat_sf)
}

acs5_2013 = append_results(acs5_2013, period = 2009:2013)
acs5_2014 = append_results(acs5_2014, period = 2010:2014)
acs5_2015 = append_results(acs5_2015, period = 2011:2015)
acs5_2016 = append_results(acs5_2016, period = 2012:2016)
acs5_2017 = append_results(acs5_2017, period = 2013:2017)
neighbs = append_results(neighbs, period = 2013:2017)

# The objective of our analysis - predictions on the four target neighborhoods.
print(neighbs)
