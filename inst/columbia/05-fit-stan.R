# Some of the coda function names conflict with rstan function names, so
# let's unload coda.
unloadNamespace("coda")
library(rstan)

# Fit the model with Stan
stan_dat = list(
	N = N, n = n, r = r, z = z_scaled, v = v_scaled, H = as.matrix(H),
	S = as.matrix(S), K = as.matrix(K),
	alpha_K = 1, beta_K = 2, alpha_xi = 1, beta_xi = 2, alpha_mu = 1, beta_mu = 2
)
stan_out = stan(file = "stcos.stan", data = stan_dat, iter = 2000, chains = 2)
print(stan_out, par = c("eta", "sig2mu", "sig2xi", "sig2K"))

# Check traceplots and histograms of draws
traceplot(stan_out, par = c("sig2mu", "sig2xi", "sig2K"), nrow = 3, ncol = 1)
traceplot(stan_out, par = sprintf("mu[%d]", 1:6))
traceplot(stan_out, par = c("eta"))
traceplot(stan_out, par = c("lp__"))

stan_hist(stan_out, par = c("sig2mu", "sig2xi", "sig2K"), bins = 30, nrow = 3, ncol = 1)
stan_hist(stan_out, par = sprintf("mu[%d]", 1:6), bins = 30)
stan_hist(stan_out, par = c("eta"), bins = 30)
stan_hist(stan_out, par = c("lp__"), bins = 30)

# Produce results on target supports.
# Compute H and S matrices and get summaries of posterior distribution for E(Y).
# Use 90% significance for all credible intervals and MOEs, following the Census
# Bureau standard.
stan_draws = extract(stan_out, pars = c("mu", "eta"), permuted = TRUE)
append_results = function(dat_sf, period, alpha = 0.10) {
	H_new = overlap_matrix(dat_sf, dom_fine)              # New overlap
	S_new_full = bs_spt$compute(dat_sf, period)           # New basis fn
	S_new = S_new_full %*% Tx_S                           # Reduce dimension

	EY_scaled = stan_draws$mu %*% t(H_new) +
		stan_draws$eta %*% t(S_new)                       # Get draws of E(Y)
	A = sd(z) * EY_scaled + mean(z)                       # Uncenter and unscale

	dat_sf$E_mean = colMeans(A)                           # Point estimates
	dat_sf$E_sd = apply(A, 2, sd)                         # SDs
	dat_sf$E_lo = apply(A, 2, quantile, prob = alpha/2)   # Interval lo
	dat_sf$E_hi = apply(A, 2, quantile, prob = 1-alpha/2) # Interval hi
	dat_sf$E_median = apply(A, 2, median)                 # Median
	dat_sf$E_moe = apply(A, 2, sd) * qnorm(1-alpha/2)     # MOE
	return(dat_sf)
}

acs5_2013_out = append_results(acs5_2013, period = 2009:2013)
acs5_2014_out = append_results(acs5_2014, period = 2010:2014)
acs5_2015_out = append_results(acs5_2015, period = 2011:2015)
acs5_2016_out = append_results(acs5_2016, period = 2012:2016)
acs5_2017_out = append_results(acs5_2017, period = 2013:2017)
nb_out = append_results(neighbs, period = 2013:2017)
