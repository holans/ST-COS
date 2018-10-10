# Analysis of City of Columbia data
set.seed(1234)

# ----- Load Fine-level Support -----
library(jsonlite)
library(sf)
library(tigris)
library(dplyr)

dom.fine <- block_groups(state = '29', county = '019', year = 2015) %>% 
	st_as_sf() %>%
	st_transform(crs = 3857) %>%
	mutate(geoid = GEOID) %>%
	dplyr::select(-GEOID)
plot(dom.fine[,1])

# ----- Assemble Source Supports -----

# https://www.census.gov/data/developers/guidance/api-user-guide.html
# https://api.census.gov/data/2015/acs5?get=NAME,B19013_001E&for=block%20group:*&in=state:29+county:019
# https://api.census.gov/data/2015/acs5?get=NAME,B19013_001M&for=block%20group:*&in=state:29+county:019
# Note: ACS 2016 data is not yet available via Census API

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
		dplyr::select(state, county, tract, blockgroup = `block group`,
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
		dplyr::select(geoid = GEOID, state = STATEFP, county = COUNTYFP,
			tract = TRACTCE, blockgroup = BLKGRPCE,
			DirectEst, DirectMOE, DirectVar)

	dat.list[[idx]] <- my_sf %>%
		filter(!is.na(DirectEst))

	dat.missing[[idx]] <- my_sf %>%
		filter(is.na(DirectEst))
}

acs5.2012 <- dat.list[[1]]
acs5.2013 <- dat.list[[2]]
acs5.2014 <- dat.list[[3]]
acs5.2015 <- dat.list[[4]]
rm(dat.list)

plot(acs5.2015[,'DirectEst'])

# ----- Load Target Supports -----
neighbs <- st_read("shp/neighborhoods.shp") %>%
	st_transform(crs = st_crs(dom.fine))
plot(neighbs)

# ----- Prepare to fit the model -----
library(fields)
library(stcos)
library(ggplot2)
library(ggforce)

# Spatial knots selected via space-filling design
u <- st_sample(dom.fine, size = 2000)
M <- matrix(unlist(u), length(u), 2, byrow = TRUE)
out <- cover.design(M, 500)
knots.sp <- out$design

# Temporal knots are selected to be evenly spaced
knots.t <- seq(2008, 2015, by = 0.5)

# Combined spatio-temporal knots
knots <- merge(knots.sp, knots.t)
basis <- SpaceTimeBisquareBasis$new(knots[,1], knots[,2], knots[,3], w.s = 1, w.t = 1)

knots.sp.dat <- data.frame(x = knots.sp[,1], y = knots.sp[,2], r = basis$get_rl())
g <- ggplot(dom.fine) +
	geom_sf(colour = "black", size = 0.25, fill = NA) +
	geom_point(data = knots.sp.dat, aes(x, y), lwd = 1, col = "red") +
	geom_point(data = knots.sp.dat[1,], aes(x, y), lwd = 3, col = "blue") +
	geom_circle(data = knots.sp.dat[1,], aes(x0=x, y0=y, r=r), fill = NA, lwd = 0.5, col = "blue") +
	labs(x = "", y = "") +
	theme_bw()
print(g)
ggsave(g, filename = "spatial-knots.pdf", width = 4, height = 7)

# Create an STCOSPrep object and add source supports
sp <- STCOSPrep$new(fine_domain = dom.fine, fine_domain_geo_name = "geoid", basis = basis, basis_mc_reps = 500)
sp$add_obs(acs5.2012, period = 2008:2012, estimate_name = "DirectEst", variance_name = "DirectVar", geo_name = "geoid")
sp$add_obs(acs5.2013, period = 2009:2013, estimate_name = "DirectEst", variance_name = "DirectVar", geo_name = "geoid")
sp$add_obs(acs5.2014, period = 2010:2014, estimate_name = "DirectEst", variance_name = "DirectVar", geo_name = "geoid")
sp$add_obs(acs5.2015, period = 2011:2015, estimate_name = "DirectEst", variance_name = "DirectVar", geo_name = "geoid")

# Read the objects needed for MCMC
z <- sp$get_z()
v <- sp$get_v()
H <- sp$get_H()
S <- sp$get_S()

# Some dimension reduction
eig <- eigen(t(S) %*% S)
rho <- eig$values

idx.S <- which(cumsum(rho) / sum(rho) < 0.65)
Tx.S <- t(eig$vectors[idx.S,])
f <- function(S) { S %*% Tx.S }
sp$set_basis_reduction(f)
S.reduced <- sp$get_reduced_S()

pdf("pca-reduction.pdf", width = 7, height = 5)
eigprops <- cumsum(rho) / sum(rho)
plot(eigprops[1:200], xlab = "Dimension", ylab = "Proportion of Variation")
abline(v = max(idx.S), lty = 2)
abline(h = eigprops[max(idx.S)], lty = 2)
dev.off()

# Pick a covariance structure for random coefficients of basis expansion
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

# Std'ize before MCMC
z.scaled <- (z - mean(z)) / sd(z)
v.scaled <- v / var(z)

# ----- Fit the Model -----
library(coda)

# Use MLE as initial value for MCMC
# This seems to fail sometimes when we don't have a lot of observations...
K <- solve(K.inv)
mle.out <- mle.stcos(z.scaled, v.scaled, H, S.reduced, K, init = list(sig2K = 1, sig2xi = 1))
init <- list(
	sig2K = mle.out$sig2K.hat,
    sig2xi = mle.out$sig2xi.hat,
    mu_B = mle.out$mu.hat
)

# Gibbs Sampler
gibbs.out <- gibbs.stcos.raw(z.scaled, v.scaled, H, S.reduced, K.inv, R = 10000,
	report.period = 2000, burn = 2000, thin = 10, init = init)
print(gibbs.out)

pdf("trace.pdf")
plot((mu_B.mcmc <- mcmc(gibbs.out$mu_B.hist))[,1:3])
plot((xi.mcmc <- mcmc(gibbs.out$xi.hist))[,1:3])
plot((eta.mcmc <- mcmc(gibbs.out$eta.hist))[,1:3])
dev.off()

pdf("trace-varcomps.pdf", width = 7, height = 5)
varcomps.mcmc <- mcmc(cbind(
	gibbs.out$sig2mu.hist,
	gibbs.out$sig2xi.hist,
	gibbs.out$sig2K.hist
))
colnames(varcomps.mcmc) <- c("sig2mu", "sig2xi", "sig2K")
plot(varcomps.mcmc)
dev.off()

# ----- Produce Results on target supports -----
# Compute H and S matrices and get summaries of posterior distribution for E(Y)
# Use 90% significance for all credible intervals and MOEs

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
print(neighbs)

# Save progress here before attempting to make plots (that can always be done later)
save.image("results.Rdata")

# ----- Make some plots -----
library(gridExtra)
library(ggrepel)

lim.est <- range(acs5.2015$DirectEst, acs5.2015$E.mean)
g <- ggplot(acs5.2015) +
	geom_sf(colour = "black", size = 0.05, aes(fill = DirectEst)) +
	ggtitle("Median Household Income\nfor Boone County", subtitle = "2015 5yr ACS Direct Estimates") +
	scale_fill_distiller("DirectEst", palette = "RdYlBu", limits = lim.est) +
	theme_bw()
h <- ggplot(acs5.2015) +
	geom_sf(colour = "black", size = 0.05, aes(fill = E.mean)) +
	ggtitle("Median Household Income\nfor Boone County", subtitle = "2015 5yr Model Estimates") +
	scale_fill_distiller("E.mean", palette = "RdYlBu", limits = lim.est) +
	theme_bw()
k <- grid.arrange(g,h, ncol = 2)
ggsave(k, filename = "compare2015.pdf")
ggsave(g, filename = "compare2015-direct.pdf", width = 5)
ggsave(h, filename = "compare2015-model.pdf", width = 5)

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
ggsave(g2015, filename = "compare-scatter.pdf", width = 5, height = 5)


idx <- which.max(abs(acs5.2015$DirectEst - acs5.2015$E.mean))
acs5.2015$geoid_outlier <- ""
acs5.2015$geoid_outlier[idx] <- acs5.2015$geoid[idx]
g <- ggplot(acs5.2015, aes(x=DirectEst, y=E.mean)) +
	geom_point(size = 2) +
	geom_abline(intercept = 0, slope = 1, color="red", 
		linetype="dashed", size=1.2) +
	ggtitle("2015 5-year Estimates") +
	labs(x = "Direct Estimate", y = "Model-Based Estimate") +
	theme_bw() +
	geom_text(aes(label=geoid_outlier), hjust=0.5, vjust=-1)
print(g)
ggsave(g, filename = "compare2015-scatter.pdf", width = 5, height = 5)

Central <- neighbs[1,]
East <- neighbs[2,]
North <- neighbs[3,]
Paris <- neighbs[4,]
Outlier <- acs5.2015[idx,]
Missing1 <- dat.missing[[4]][1,]
Missing2 <- dat.missing[[4]][2,]
Missing3 <- dat.missing[[4]][3,]
Missing4 <- dat.missing[[4]][4,]

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
	ggtitle("Median Household Income for Boone County", subtitle = "ACS 2015 5yr Direct Estimates") +
	scale_fill_distiller("DirectEst", palette = "RdYlBu", limits = lim.est) +
	geom_sf(data = neighbs, fill = "black") +
	geom_label_repel(data = st_centroid(East), nudge_x = 20000, nudge_y = 130000, aes(x=East.coord[1], y=East.coord[2], label="East")) +
	geom_label_repel(data = st_centroid(Central), nudge_x = -10000, nudge_y = 130000, aes(x=Central.coord[1], y=Central.coord[2], label="Central")) +
	geom_label_repel(data = st_centroid(North), nudge_x = 0, nudge_y = 100000, aes(x=North.coord[1], y=North.coord[2], label="North")) +
	geom_label_repel(data = st_centroid(Paris), nudge_x = 10000, nudge_y = 100000, aes(x=Paris.coord[1], y=Paris.coord[2], label="Paris")) +
	geom_label_repel(data = st_centroid(Outlier), nudge_x = 0, nudge_y = -50000, aes(x=Outlier.coord[1], y=Outlier.coord[2], label="Outlier")) +
	geom_label_repel(data = st_centroid(Missing1), nudge_x = 100000, nudge_y = -50000, aes(x=Missing1.coord[1], y=Missing1.coord[2], label="Missing1")) +
	geom_label_repel(data = st_centroid(Missing2), nudge_x = 100000, nudge_y = -36000, aes(x=Missing2.coord[1], y=Missing2.coord[2], label="Missing2")) +
	geom_label_repel(data = st_centroid(Missing3), nudge_x = -100000, nudge_y = -50000, aes(x=Missing3.coord[1], y=Missing3.coord[2], label="Missing3")) +
	geom_label_repel(data = st_centroid(Missing4), nudge_x = -100000, nudge_y = -38000, aes(x=Missing4.coord[1], y=Missing4.coord[2], label="Missing4")) +
	xlab(NULL) +
	ylab(NULL) +
	theme_bw()
print(g)
ggsave(g, filename = "areas-of-interest-map.pdf")


