library(coda)

# Run the Gibbs sampler.
hyper = list(a_sig2K = 1, b_sig2K = 2, a_sig2xi = 1, b_sig2xi = 2,
	a_sig2mu = 1, b_sig2mu = 2)
gibbs_out = gibbs_stcos(z = z_scaled, v = v_scaled, H = H, S = S,
	Kinv = Kinv, R = 10000, report_period = 2000, burn = 2000,
	thin = 10, hyper = hyper)
print(gibbs_out)

# Show some trace plots to assess convergence of the sampler.
plot((muB_mcmc = mcmc(gibbs_out$muB_hist))[,1:3])
plot((xi_mcmc = mcmc(gibbs_out$xi_hist))[,1:3])
plot((eta_mcmc = mcmc(gibbs_out$eta_hist))[,1:3])

varcomps_mcmc = mcmc(data.frame(
	sig2mu = gibbs_out$sig2mu_hist,
	sig2xi = gibbs_out$sig2xi_hist,
	sig2K = gibbs_out$sig2K_hist
))
plot(varcomps_mcmc)

# Same plots as above, but save individual plots to files.
for (j in 1:ncol(varcomps_mcmc)) {
	pdf(sprintf("trace-%s.pdf", colnames(varcomps_mcmc)[j]), width = 5, height = 4)
	plot(varcomps_mcmc[,j], trace = TRUE, density = FALSE,
		main = sprintf("Trace of %s", colnames(varcomps_mcmc)[j]))
	dev.off()

	pdf(sprintf("density-%s.pdf", colnames(varcomps_mcmc)[j]), width = 5, height = 4)
	plot(varcomps_mcmc[,j], trace = FALSE, density = TRUE,
		main = sprintf("Density of %s", colnames(varcomps_mcmc)[j]))
	dev.off()
}

# Produce results on target supports.
# Compute H and S matrices and get summaries of posterior distribution for E(Y).
# Use 90% significance for all credible intervals and MOEs, following the Census
# Bureau standard.
append_results = function(dat_sf, period, alpha = 0.10) {
	H_new = overlap_matrix(dat_sf, dom_fine)              # New overlap
	S_new_full = bs_spt$compute(dat_sf, period)           # New basis fn
	S_new = S_new_full %*% Tx_S                           # Reduce dimension

	EY_scaled = fitted(gibbs_out, H_new, S_new)           # Get draws of E(Y)
	A = sd(z) * EY_scaled + mean(z)                       # Uncenter and unscale

	dat_sf$E_mean = colMeans(A)						      # Point estimates
	dat_sf$E_sd = apply(A, 2, sd)						  # SDs
	dat_sf$E_lo = apply(A, 2, quantile, prob = alpha/2)   # Credible interval lo
	dat_sf$E_hi = apply(A, 2, quantile, prob = 1-alpha/2) # Credible interval hi
	dat_sf$E_median = apply(A, 2, median)				  # Median
	dat_sf$E_moe = apply(A, 2, sd) * qnorm(1-alpha/2)	  # MOE
	return(dat_sf)
}

acs5_2013_out = append_results(acs5_2013, period = 2009:2013)
acs5_2014_out = append_results(acs5_2014, period = 2010:2014)
acs5_2015_out = append_results(acs5_2015, period = 2011:2015)
acs5_2016_out = append_results(acs5_2016, period = 2012:2016)
acs5_2017_out = append_results(acs5_2017, period = 2013:2017)
nb_out = append_results(neighbs, period = 2013:2017)
