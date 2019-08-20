library(tigris)
library(ggplot2)
library(dplyr)
library(fields)
library(stcos)
library(ggforce)

set.seed(1234)

# ----- Load data for analysis -----
# Load the fine-level support via shapefile using the `tigris` package.
# We use the 2017 block-groups in Boone County, MO. Convert to an `sf`
# object, and transform to the projection with EPSG code 3857.
options(tigris_use_cache = TRUE)
options(tigris_refresh = FALSE)
dom.fine = block_groups(state = '29', county = '019', year = 2017) %>%
	st_as_sf() %>%
	st_transform(crs = 3857) %>%
	mutate(geoid = GEOID) %>%
	select(-GEOID)

# Load previously constructed source supports.
data(acs_sf)

# Load Target Supports
# Make sure to transform to the same projection as the fine-level support.
data(columbia_neighbs)
neighbs = columbia_neighbs %>% st_transform(crs = st_crs(dom.fine))

# ----- Prepare basis for model -----
# Select spatial knots via space-filling design. Start with uniformly
# drawn points from the fine-level geography and use `cover.design`
# to select a subset of those points for knots.
u = st_sample(dom.fine, size = 2000)
M = matrix(unlist(u), length(u), 2, byrow = TRUE)
out = cover.design(M, 500)
knots.sp = out$design

# Select temporal knots to be evenly spaced over the years relevant
# to the source support years.
knots.t = seq(2009, 2017, by = 0.5)

# Use Cartesian join of spatial and temporal knots to obtain spatio-temporal knots.
knots = merge(knots.sp, knots.t)
r_full = nrow(knots)

# Create an SpaceTimeBisquareBasis object with our knot points.
# The STCOSPrep object makes use of this basis when pass it a support.
basis = SpaceTimeBisquareBasis$new(knots[,1], knots[,2], knots[,3], w.s = 1, w.t = 1)

# ----- Build components for STCOS model -----
# Create an STCOSPrep object and add source supports
sp = STCOSPrep$new(fine_domain = dom.fine, fine_domain_geo_name = "geoid",
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

# Read the objects needed for MCMC
z = sp$get_z()
v = sp$get_v()
H = sp$get_H()
S = sp$get_S()

# Dimension reduction of `S` matrix via PCA.
eig = eigen(t(S) %*% S)
rho = eig$values

idx.S = which(cumsum(rho) / sum(rho) < 0.65)
Tx = eig$vectors[,idx.S]
f = function(S) { S %*% Tx }
sp$set_basis_reduction(f)
S.reduced = sp$get_reduced_S()

# Pick a covariance structure for random coefficients of basis expansion.
# Construct K with a "Random Walk" covariance structure.
K.inv = sp$get_Kinv(2009:2017, method = "randomwalk")
K = solve(K.inv)

# Standardize observations before MCMC.
z.mean = mean(z, na.rm = TRUE)
z.sd = sd(z, na.rm = TRUE)
z.scaled = (z - z.mean) / z.sd
v.scaled = v / z.sd^2

# ----- Run the Gibbs sampler and produce outputs -----
gibbs.out = gibbs.stcos.raw(z = z.scaled, v = v.scaled, H = H,
	S = S.reduced, K.inv = K.inv, R = 10000, report.period = 2000,
	burn = 2000, thin = 10)
print(gibbs.out)

append_results = function(dat_sf, period, geo_name, alpha = 0.10) {
	out = sp$domain2model(dat_sf, period = period, geo_name = geo_name)
	E.hat.scaled = fitted(gibbs.out, out$H, out$S.reduced)
	E.hat = z.sd * E.hat.scaled + z.mean                      # Uncenter and unscale
	dat_sf$E.mean = colMeans(E.hat)                           # Point estimates
	dat_sf$E.sd = apply(E.hat, 2, sd)                         # SDs
	dat_sf$E.lo = apply(E.hat, 2, quantile, prob = alpha/2)   # Credible interval lo
	dat_sf$E.hi = apply(E.hat, 2, quantile, prob = 1-alpha/2) # Credible interval hi
	dat_sf$E.median = apply(E.hat, 2, median)                 # Median
	dat_sf$E.moe = apply(E.hat, 2, sd) * qnorm(1-alpha/2)     # MOE
	return(dat_sf)
}

acs5_2013 = append_results(acs5_2013, period = 2009:2013, geo_name = "geoid")
acs5_2014 = append_results(acs5_2014, period = 2010:2014, geo_name = "geoid")
acs5_2015 = append_results(acs5_2015, period = 2011:2015, geo_name = "geoid")
acs5_2016 = append_results(acs5_2016, period = 2012:2016, geo_name = "geoid")
acs5_2017 = append_results(acs5_2017, period = 2013:2017, geo_name = "geoid")
neighbs = append_results(neighbs, 2013:2017, geo_name = "Region")

# The objective of our analysis - predictions on the four target neighborhoods.
print(neighbs)
