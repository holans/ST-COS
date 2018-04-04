library(fields)
library(ggplot2)
library(ggforce)
library(dplyr)
library(readr)
library(sf)
library(stcos)

num.sp.knots <- 500
radius <- 1
pr.rho <- 0.75
gridn <- 16
sample.prop <- 0.005
reps <- 1

state <- 42
extract <- function(dat) {
	dat[dat$STATE == state,]
}

# ----- Load space-time domains with direct survey estimates/variances -----
# Fine-level support
county_acs_5yr2015 <- st_read("../acs_county_S1901/merged/county_acs_5yr2015.shp")
dom.fine <- extract(county_acs_5yr2015)

# ----- Set up knots for bisquare basis -----
# Knots must be compatible (same projection) with fine-level geography
u <- st_sample(dom.fine, size = 5000)
M <- matrix(unlist(u), length(u), 2, byrow = TRUE)
out <- cover.design(M, num.sp.knots)
knots.sp <- out$design

knots.t <- c(2005, 2005.5, 2006, 2006.5, 2007, 2007.5, 2008, 2008.5,
	2009, 2009.5, 2010, 2010.5, 2011, 2011.5, 2012, 2012.5, 2013, 2013.5,
	2014, 2014.5, 2015)
knots <- merge(knots.sp, knots.t)
names(knots) <- c("x", "y", "t")
basis <- SpaceTimeBisquareBasis$new(knots[,1], knots[,2], knots[,3], w.s = radius, w.t = 1)

pdf("design.pdf")
plot(dom.fine[,1], col = "white")
points(knots.sp)
dev.off()

#knots.sp.dat <- data.frame(x = knots.sp[,1], y = knots.sp[,2], r = basis$get_rl())
#g <- ggplot(dom.fine) +
#	geom_sf(colour = "black", size = 0.05) +
#	geom_point(data = knots.sp.dat, aes(x, y), lwd = 1, col = "blue") +
#	geom_circle(data = knots.sp.dat, aes(x0=x, y0=y, r=r), fill = NA, lwd = 0.1, col = "blue") +
#	labs(x = "", y = "") +
#	theme_bw()
#print(g)
#ggsave("design.pdf", width = 7, height = 5)
#dev.off()

plan.db <- read_csv("../planning_2016//pdb2016bgv8_us.csv")
head(plan.db)
# plan.db <- read.csv("planning/pdb2016bgv8_us.csv")
# head(plan.db)

cb_2016_42_bg <- st_read("cb_2016_42_bg_500k/cb_2016_42_bg_500k.shp")
cb_2016_42_bg <- st_transform(cb_2016_42_bg, crs = st_crs(dom.fine))

plan.db <- plan.db %>%
	mutate(State = as.character(State)) %>%
	mutate(County = as.character(County)) %>%
	mutate(Tract = as.character(Tract)) %>%
	mutate(Block_group = as.character(Block_group)) %>%
	mutate(Med_HHD = sub("\\$", "", Med_HHD_Inc_BG_ACS_10_14)) %>%
	mutate(Med_HHD = sub(",", "", Med_HHD)) %>%
	mutate(Med_HHD = as.numeric(Med_HHD)) %>%
	mutate(Med_HHD_MOE = sub("\\$", "", Med_HHD_Inc_BG_ACSMOE_10_14)) %>%
	mutate(Med_HHD_MOE = sub(",", "", Med_HHD_MOE)) %>%
	mutate(Med_HHD_MOE = as.numeric(Med_HHD_MOE)) %>%
	dplyr::select(State, County, Tract, Block_group, Med_HHD, Med_HHD_MOE, Tot_Housing_Units_CEN_2010)

cb_2016_42_bg <- cb_2016_42_bg %>%
	mutate(STATEFP = as.character(STATEFP)) %>%
	mutate(COUNTYFP = as.character(COUNTYFP)) %>%
	mutate(TRACTCE = as.character(TRACTCE)) %>%
	mutate(BLKGRPCE = as.character(BLKGRPCE)) %>%
	left_join(plan.db, by = c("STATEFP" = "State", "COUNTYFP" = "County", "TRACTCE" = "Tract", "BLKGRPCE" = "Block_group"))

# It takes quite a while to draw from areal units, so we'll just draw the
# locations once at the beginning of the simulation. We'll redraw the survey'd
# variables within each rep

N.total <- sum(cb_2016_42_bg$Tot_Housing_Units_CEN_2010[
	cb_2016_42_bg$Tot_Housing_Units_CEN_2010 > 0 &
	!is.na(cb_2016_42_bg$Med_HHD) &
	!is.na(cb_2016_42_bg$Med_HHD_MOE)])
na <- rep(NA, N.total)
sim.popn <- data.frame(xcoord = na, ycoord = na, state = NA,
	county = NA, tract = NA, block_group = NA, geo_id = NA,
	mean_y = NA, sd_y = NA, y = NA)
geo2units <- list()
last.idx <- 0
for (i in 1:nrow(cb_2016_42_bg)) {
	if (i %% 100 == 0) {
		logger("Generating data for block group %d of %d\n", i, nrow(cb_2016_42_bg))
	}

	N.blg <- cb_2016_42_bg$Tot_Housing_Units_CEN_2010[i]
	if (N.blg == 0 || is.na(cb_2016_42_bg$Med_HHD[i]) || is.na(cb_2016_42_bg$Med_HHD_MOE[i])) {
		logger("Skipping block group %d of %d with GEOID %s\n", i, nrow(cb_2016_42_bg), cb_2016_42_bg$GEOID[i])
	} else {
		pp <- rArea(N.blg, cb_2016_42_bg[i,], blocksize = 2*N.blg)

		idx <- 1:N.blg + last.idx
		sim.popn$xcoord[idx] <- pp[,1]
		sim.popn$ycoord[idx] <- pp[,2]
		sim.popn$state[idx] <- cb_2016_42_bg$STATEFP[i]
		sim.popn$county[idx] <- cb_2016_42_bg$COUNTYFP[i]
		sim.popn$tract[idx] <- cb_2016_42_bg$TRACTCE[i]
		sim.popn$block_group[idx] <- cb_2016_42_bg$BLKGRPCE[i]
		sim.popn$geo_id[idx] <- cb_2016_42_bg$GEOID[i]
		
		sim.popn$mean_y[idx] <- cb_2016_42_bg$Med_HHD[i]
		sim.popn$sd_y[idx] <- cb_2016_42_bg$Med_HHD_MOE[i] / 1.645

		geo2units[[cb_2016_42_bg$GEOID[i]]] <- idx
		last.idx <- last.idx + N.blg
	}
}
save.image("popn.Rdata")

# ----- Construct a STCOSPrep object and add space-time domains -----
# Do this without having generated DirectEst and DirectVar first
year.levels <- 2005:2015
sp <- STCOSPrep$new(fine_domain = dom.fine, fine_domain_geo_name = "GEO_ID", basis = basis, basis_mc_reps = 500)

source.supps.1yr <- list()
for (j in 1:length(year.levels)) {
	year <- year.levels[j]
	supp <- dom.fine
	supp$DirectEst <- NA
	supp$DirectVar <- NA
	lookback <- year
	sp$add_obs(supp, period = lookback, estimate_name = "DirectEst",
		variance_name = "DirectVar", geo_name = "GEO_ID")
	source.supps.1yr[[as.character(year)]] <- supp
}

source.supps.3yr <- list()
for (j in 1:length(year.levels)) {
	if (j < 3) next
	year <- year.levels[j]
	supp <- dom.fine
	supp$DirectEst <- NA
	supp$DirectVar <- NA
	lookback <- seq(year-2, year)
	sp$add_obs(supp, period = lookback, estimate_name = "DirectEst",
		variance_name = "DirectVar", geo_name = "GEO_ID")
	source.supps.3yr[[as.character(year)]] <- supp
}

source.supps.5yr <- list()
for (j in 1:length(year.levels)) {
	if (j < 5) next
	year <- year.levels[j]
	supp <- dom.fine
	supp$DirectEst <- NA
	supp$DirectVar <- NA
	lookback <- seq(year-4, year)
	sp$add_obs(supp, period = lookback, estimate_name = "DirectEst",
		variance_name = "DirectVar", geo_name = "GEO_ID")
	source.supps.5yr[[as.character(year)]] <- supp
}

# ----- Dimension reduction -----
logger("Runing dimension reduction on S")
S <- sp$get_S()
eig <- eigen(t(S) %*% S)
rho <- eig$values

pdf("cols.pdf", width = 5, height = 5)
plot(cumsum(rho) / sum(rho), ylab = "Proportion of Eigenvalues", xlab = "Number of Eigenvalues")

idx.S <- which(cumsum(rho) / sum(rho) < 0.60)
abline(v = max(idx.S), lty = 2, col = "red", lwd = 2)
axis(side = 3, at = max(idx.S), lwd = 2, padj = 0.5, hadj = 0.75)

idx.S <- which(cumsum(rho) / sum(rho) < 0.75)
abline(v = max(idx.S), lty = 2, col = "red", lwd = 2)
axis(side = 3, at = max(idx.S), lwd = 2, padj = -1)

idx.S <- which(cumsum(rho) / sum(rho) < 0.90)
abline(v = max(idx.S), lty = 2, col = "red", lwd = 2)
axis(side = 3, at = max(idx.S), lwd = 2, padj = 0.5)
dev.off()

idx.S <- which(cumsum(rho) / sum(rho) < pr.rho)
Tx.S <- t(eig$vectors[idx.S,])
f <- function(S) { S %*% Tx.S }
sp$set_basis_reduction(f)

# ----- Compute Kinv -----
S.reduced <- sp$get_reduced_S()
K.inv <- diag(x = 1, nrow = ncol(S.reduced))

# ----- Set up target supports -----
logger("Setting up target geography\n")
cb_2016_us_cd115 <- st_read("cb_2016_us_cd115_500k/cb_2016_us_cd115_500k.shp") %>%
	filter(STATEFP == '42') %>%
	st_transform(crs = st_crs(dom.fine))
# cb_2016_us_cd115 <- st_transform(cb_2016_us_cd115, crs = st_crs(dom.fine))
plot(dom.fine[,1])
plot(cb_2016_us_cd115[,1], col = "NA", add = TRUE, lwd = 2, lty = 2)
plot(cb_2016_us_cd115[,4])

# See https://gis.stackexchange.com/questions/225157/generate-rectangular-fishnet-or-vector-grid-cells-shapefile-in-r/243585
# Remove squares that intersect with the outside of the state
# Also remove Erie County (049) because much of it is water
state.grid <- st_make_grid(extract(dom.fine), n = c(gridn, gridn), what = 'polygons') %>%
	st_sf('geometry' = ., data.frame('GEO_ID' = 1:length(.))) %>%
	st_transform(st_crs(dom.fine))
d1 <- state.grid %>%
	st_join(county_acs_5yr2015) %>%
	filter(STATE != '42' | COUNTY == '049' | is.na(COUNTY)) %>%
	group_by(GEO_ID.x) %>%
	summarize(count = length(GEO_ID.x)) %>%
	mutate(GEO_ID = GEO_ID.x) %>%
	select(GEO_ID, count)
idx.drop <- which(state.grid$GEO_ID %in% d1$GEO_ID)
state.grid <- state.grid[-idx.drop,]
plot(dom.fine[,1])
plot(state.grid[,1], col = NA, add = TRUE)

# See https://gis.stackexchange.com/questions/225157/generate-rectangular-fishnet-or-vector-grid-cells-shapefile-in-r/243585
# Remove squares that intersect with the outside of the state
# Also remove Erie County (049) because much of it is water
state.grid.finer <- st_make_grid(extract(dom.fine), n = c(gridn*2, gridn*2), what = 'polygons') %>%
	st_sf('geometry' = ., data.frame('GEO_ID' = 1:length(.))) %>%
	st_transform(st_crs(dom.fine))
d1 <- state.grid.finer %>%
	st_join(county_acs_5yr2015) %>%
	filter(STATE != '42' | COUNTY == '049' | is.na(COUNTY)) %>%
	group_by(GEO_ID.x) %>%
	summarize(count = length(GEO_ID.x)) %>%
	mutate(GEO_ID = GEO_ID.x) %>%
	select(GEO_ID, count)
idx.drop <- which(state.grid.finer$GEO_ID %in% d1$GEO_ID)
state.grid.finer <- state.grid.finer[-idx.drop,]
plot(dom.fine[,1])
plot(state.grid.finer[,1], col = NA, add = TRUE)

# See https://gis.stackexchange.com/questions/225157/generate-rectangular-fishnet-or-vector-grid-cells-shapefile-in-r/243585
# Remove squares that intersect with the outside of the state
# Also remove Erie County (049) because much of it is water
state.grid.coarser <- st_make_grid(extract(dom.fine), n = c(gridn/2, gridn/2), what = 'polygons') %>%
	st_sf('geometry' = ., data.frame('GEO_ID' = 1:length(.))) %>%
	st_transform(st_crs(dom.fine))
d1 <- state.grid.coarser %>%
	st_join(county_acs_5yr2015) %>%
	filter(STATE != '42' | COUNTY == '049' | is.na(COUNTY)) %>%
	group_by(GEO_ID.x) %>%
	summarize(count = length(GEO_ID.x)) %>%
	mutate(GEO_ID = GEO_ID.x) %>%
	select(GEO_ID, count)
idx.drop <- which(state.grid.coarser$GEO_ID %in% d1$GEO_ID)
state.grid.coarser <- state.grid.coarser[-idx.drop,]
plot(dom.fine[,1])
plot(state.grid.coarser[,1], col = NA, add = TRUE)

target.out <- sp$domain2model(cb_2016_us_cd115, period = 2016, geo_name = 'GEOID')
H.targ <- target.out$H
S.reduced.targ <- target.out$S.reduced

target2.out <- sp$domain2model(dom.fine, period = 2016, geo_name = 'GEO_ID')
H.targ2 <- target2.out$H
S.reduced.targ2 <- target2.out$S.reduced

target3.out <- sp$domain2model(state.grid, period = 2016, geo_name = 'GEO_ID')
H.targ3 <- target3.out$H
S.reduced.targ3 <- target3.out$S.reduced

target4.out <- sp$domain2model(state.grid.finer, period = 2016, geo_name = 'GEO_ID')
H.targ4 <- target4.out$H
S.reduced.targ4 <- target4.out$S.reduced

target5.out <- sp$domain2model(state.grid.coarser, period = 2016, geo_name = 'GEO_ID')
H.targ5 <- target5.out$H
S.reduced.targ5 <- target5.out$S.reduced

RES1 <- matrix(NA, nrow(H.targ), reps)
RES2 <- matrix(NA, nrow(H.targ2), reps)
RES3 <- matrix(NA, nrow(H.targ3), reps)
RES4 <- matrix(NA, nrow(H.targ4), reps)
RES5 <- matrix(NA, nrow(H.targ5), reps)

for (rep in 1:reps)
{
	logger("Starting simulation rep %d\n", rep)

	logger("Generating unit level data\n")
	# We'll generate data at the block-group level, but model at the county level.
	# Need population sizes for each block-group
	# Need average and stddev income for block-groups
	# Need shape file for block-groups

	# Assume mean_y and sd_y are parameters for lognormal rvs.
	# Solve for mean and sd of normal rvs.
	# Then generate normal rvs representing the data on a log scale
	var.normal <- log(1 + sim.popn$sd_y^2 / sim.popn$mean_y^2)
	mean.normal <- log(sim.popn$mean_y) - 1/2 * var.normal
	sim.popn$y <- rnorm(nrow(sim.popn), mean.normal, sqrt(var.normal))
	# hist(exp(sim.popn$y[sim.popn$y < log(200000)]))

	# Aggregate simulated y to each target support
	totals1 <- sim.popn %>%
		st_as_sf(coords = c("xcoord", "ycoord"), crs = st_crs(dom.fine)) %>%
		st_join(cb_2016_us_cd115) %>%
		group_by(state, CD115FP) %>%
		summarize(y_bar = mean(y), N = length(y))
	st_geometry(totals1) <- NULL
	order(cb_2016_us_cd115$CD115FP)
	totals1 <- totals1[order(cb_2016_us_cd115$CD115FP),]

	totals2 <- sim.popn %>%
		group_by(state, county) %>%
		summarize(y_bar = mean(y), N = length(y))
	totals2 <- totals2[order(dom.fine$COUNTY),]

	totals3 <- sim.popn %>%
		st_as_sf(coords = c("xcoord", "ycoord"), crs = st_crs(dom.fine)) %>%
		st_join(state.grid) %>%
		filter(!is.na(GEO_ID)) %>%
		group_by(GEO_ID) %>%
		summarize(y_bar = mean(y), N = length(y))
	st_geometry(totals3) <- NULL
	totals3 <- totals3[order(state.grid$GEO_ID),]

	totals4 <- sim.popn %>%
		st_as_sf(coords = c("xcoord", "ycoord"), crs = st_crs(dom.fine)) %>%
		st_join(state.grid.finer) %>%
		filter(!is.na(GEO_ID)) %>%
		group_by(GEO_ID) %>%
		summarize(y_bar = mean(y), N = length(y))
	st_geometry(totals4) <- NULL
	totals4 <- totals4[order(state.grid.finer$GEO_ID),]

	totals5 <- sim.popn %>%
		st_as_sf(coords = c("xcoord", "ycoord"), crs = st_crs(dom.fine)) %>%
		st_join(state.grid.coarser) %>%
		filter(!is.na(GEO_ID)) %>%
		group_by(GEO_ID) %>%
		summarize(y_bar = mean(y), N = length(y))
	st_geometry(totals5) <- NULL
	totals5 <- totals5[order(state.grid.coarser$GEO_ID),]


	# ----- Take a sample to get estimates and variances -----
	# First let's assume that all estimates can be disclosed.
	logger("Sample to get estimates\n")

	# TBD: Use appropriate HT estimates for the mean and variance
	logger("Taking samples for 1-year estimates\n")
	sizes <- sim.popn %>%
		group_by(county) %>%
		summarize(popn = length(y)) %>%
		mutate(sample = ceiling(sample.prop * popn))
	na <- rep(NA, sum(sizes$sample) * length(year.levels))
	all.samples <- data.frame(year = na, county = na, unit = na) %>%
		mutate(county = as.character(county))
	last.idx <- 0
	for (j in 1:length(year.levels)) {
		year <- year.levels[j]
		for (l in 1:nrow(supp)) {
			county <- supp$COUNTY[l]
			idx.county <- which(sim.popn$county == county)
			idx.sample <- sample(idx.county, size = sizes$sample[l], replace = TRUE)
			source.supps.1yr[[as.character(year)]]$DirectEst[l] <- mean(sim.popn$y[idx.sample])
			source.supps.1yr[[as.character(year)]]$DirectVar[l] <- var(sim.popn$y[idx.sample])

			ind <- 1:sizes$sample[l] + last.idx
			all.samples$year[ind] <- year
			all.samples$county[ind] <- county
			all.samples$unit[ind] <- idx.sample
			last.idx <- last.idx + sizes$sample[l]
		}
	}

	logger("Taking samples for 3-year estimates\n")
	for (j in 1:length(year.levels)) {
		if (j < 3) { next }
		year <- year.levels[j]
		lookback <- year.levels[seq(j-2, j)]
		sample.3yr <- all.samples %>%
			filter(year %in% lookback)
		dat.estimates <- sim.popn[sample.3yr$unit,] %>%
			group_by(county) %>%
			summarize(direct_est = mean(y), direct_var = var(y))
		supp <- dom.fine %>%
			mutate(COUNTY = as.character(COUNTY)) %>%
			inner_join(dat.estimates, by = c("COUNTY" = "county")) %>%
			mutate(DirectEst = direct_est) %>%
			mutate(DirectVar = direct_var) %>%
			dplyr::select(-c(direct_est, direct_var))
		source.supps.3yr[[as.character(year)]]$DirectEst <- supp$DirectEst
		source.supps.3yr[[as.character(year)]]$DirectVar <- supp$DirectVar
	}

	logger("Taking samples for 5-year estimates\n")
	for (j in 1:length(year.levels)) {
		if (j < 5) { next }
		year <- year.levels[j]
		lookback <- year.levels[seq(j-4, j)]
		sample.5yr <- all.samples %>%
			filter(year %in% lookback)
		dat.estimates <- sim.popn[sample.5yr$unit,] %>%
			group_by(county) %>%
			summarize(direct_est = mean(y), direct_var = var(y))
		supp <- dom.fine %>%
			mutate(COUNTY = as.character(COUNTY)) %>%
			inner_join(dat.estimates, by = c("COUNTY" = "county")) %>%
			mutate(DirectEst = direct_est) %>%
			mutate(DirectVar = direct_var) %>%
			dplyr::select(-c(direct_est, direct_var))
		source.supps.5yr[[as.character(year)]]$DirectEst <- supp$DirectEst
		source.supps.5yr[[as.character(year)]]$DirectVar <- supp$DirectVar
	}

	# save.image("temp.Rdata")

	Z <- c(
		unlist(Map(f = function(x) { x$DirectEst }, source.supps.1yr)),
		unlist(Map(f = function(x) { x$DirectEst }, source.supps.3yr)),
		unlist(Map(f = function(x) { x$DirectEst }, source.supps.5yr))
	)
	V <- c(
		unlist(Map(f = function(x) { x$DirectVar }, source.supps.1yr)),
		unlist(Map(f = function(x) { x$DirectVar }, source.supps.3yr)),
		unlist(Map(f = function(x) { x$DirectVar }, source.supps.5yr))
	)
	H <- sp$get_H()

	# Stdize before MCMC
	D <- Diagonal(n = length(Z), x = 1/sd(Z))
	Z.scaled <- (Z - mean(Z)) / sd(Z)
	V.scaled <- V / var(Z)

	# ----- Run MCMC -----
	# Compute MLE as initial value
	mle.out <- mle.stcos(Z.scaled, S.reduced, V.scaled, H, init = list(sig2xi = 1))
	init <- list(
		sig2xi = mle.out$sig2xi.hat,
		mu_B = mle.out$mu.hat,
		eta = mle.out$eta.hat
	)
	gibbs.out <- gibbs.stcos.raw(Z.scaled, S.reduced, V.scaled, K.inv, H,
		R = 10000, report.period = 1000, burn = 1000, thin = 10,
		init = init)
	if (FALSE) {
		print(gibbs.out)
		plot(sig2mu.mcmc <- mcmc(gibbs.out$sig2mu.hist))
		plot(sig2xi.mcmc <- mcmc(gibbs.out$sig2xi.hist))
		plot(sig2K.mcmc <- mcmc(gibbs.out$sig2K.hist))
	}

	# Save results for this rep. At the end we can compare posterior expectation
	# of Y on each target support of interest with point level data. For example,
	# we can compute the true mean at the target support level, or recompute
	# estimates of the mean using sampled data.
	alpha <- 0.05

	# Posterior distribution for E(Y), applied to target support 1 (CDs)
	E.hat.scaled <- fitted(gibbs.out, H.targ, S.reduced.targ)
	E.hat <- sd(Z) * E.hat.scaled + mean(Z)
	RES1[,rep] <- colMeans(E.hat)
	if (FALSE) {
		cb_2016_us_cd115$E.mean <- colMeans(E.hat)
		cb_2016_us_cd115$E.sd <- apply(E.hat, 2, sd)
		cb_2016_us_cd115$E.lo <- apply(E.hat, 2, quantile, prob = alpha/2)
		cb_2016_us_cd115$E.hi <- apply(E.hat, 2, quantile, prob = 1 - alpha/2)
		plot(exp(cb_2016_us_cd115$y_bar), exp(cb_2016_us_cd115$E.mean))
		plot(cb_2016_us_cd115$y_bar, cb_2016_us_cd115$E.mean)
		plot(cb_2016_us_cd115[,c("y_bar","E.mean")])
	}
	
	# Posterior distribution for E(Y), applied to target support 2 (counties)
	E.hat.scaled <- fitted(gibbs.out, H.targ2, S.reduced.targ2)
	E.hat <- sd(Z) * E.hat.scaled + mean(Z)
	RES2[,rep] <- colMeans(E.hat)
	if (FALSE) {
		dom.fine$E.mean <- colMeans(E.hat)
		dom.fine$E.sd <- apply(E.hat, 2, sd)
		dom.fine$E.lo <- apply(E.hat, 2, quantile, prob = alpha/2)
		dom.fine$E.hi <- apply(E.hat, 2, quantile, prob = 1 - alpha/2)
		plot(exp(dom.fine$y_bar), exp(dom.fine$E.mean))
		plot(dom.fine$y_bar, dom.fine$E.mean)
		plot(dom.fine[,c("y_bar","E.mean")])
	}

	# Posterior distribution for E(Y), applied to target support 3 (grid spaces)
	E.hat.scaled <- fitted(gibbs.out, H.targ3, S.reduced.targ3)
	E.hat <- sd(Z) * E.hat.scaled + mean(Z)
	RES3[,rep] <- colMeans(E.hat)
	if (FALSE) {
		state.grid$E.mean <- colMeans(E.hat)
		state.grid$E.sd <- apply(E.hat, 2, sd)
		state.grid$E.lo <- apply(E.hat, 2, quantile, prob = alpha/2)
		state.grid$E.hi <- apply(E.hat, 2, quantile, prob = 1 - alpha/2)
		plot(exp(state.grid$y_bar), exp(state.grid$E.mean))
		plot(state.grid$y_bar, state.grid$E.mean)
		plot(state.grid[,c("y_bar","E.mean")])
	}

	# Posterior distribution for E(Y), applied to target support 4 (finer grid spaces)
	E.hat.scaled <- fitted(gibbs.out, H.targ4, S.reduced.targ4)
	E.hat <- sd(Z) * E.hat.scaled + mean(Z)
	RES4[,rep] <- colMeans(E.hat)
	if (FALSE) {
		state.grid.finer$E.mean <- colMeans(E.hat)
		state.grid.finer$E.sd <- apply(E.hat, 2, sd)
		state.grid.finer$E.lo <- apply(E.hat, 2, quantile, prob = alpha/2)
		state.grid.finer$E.hi <- apply(E.hat, 2, quantile, prob = 1 - alpha/2)
		plot(exp(state.grid.finer$y_bar), exp(state.grid.finer$E.mean))
		plot(state.grid.finer$y_bar, state.grid.finer$E.mean)
		plot(state.grid.finer[,c("y_bar","E.mean")])
	}

	# Posterior distribution for E(Y), applied to target support 5 (coarser grid spaces)
	E.hat.scaled <- fitted(gibbs.out, H.targ5, S.reduced.targ5)
	E.hat <- sd(Z) * E.hat.scaled + mean(Z)
	RES5[,rep] <- colMeans(E.hat)
	if (FALSE) {
		state.grid.coarser$E.mean <- colMeans(E.hat)
		state.grid.coarser$E.sd <- apply(E.hat, 2, sd)
		state.grid.coarser$E.lo <- apply(E.hat, 2, quantile, prob = alpha/2)
		state.grid.coarser$E.hi <- apply(E.hat, 2, quantile, prob = 1 - alpha/2)
		plot(exp(state.grid.coarser$y_bar), exp(state.grid.coarser$E.mean))
		plot(state.grid.coarser$y_bar, state.grid.coarser$E.mean)
		plot(state.grid.coarser[,c("y_bar","E.mean")])
	}

	# save.image("temp.Rdata")
}

totals1$E_hat <- RES1
totals2$E_hat <- RES2
totals3$E_hat <- RES3
totals4$E_hat <- RES4
totals5$E_hat <- RES5
totals1$err <- RES1 - totals1$y_bar
totals2$err <- RES2 - totals2$y_bar
totals3$err <- RES3 - totals3$y_bar
totals4$err <- RES4 - totals4$y_bar
totals5$err <- RES5 - totals5$y_bar

totals1_sf <- cb_2016_us_cd115 %>%
	inner_join(totals1, by = c("STATEFP" = "state", "CD115FP" = "CD115FP"))
totals2_sf <- dom.fine %>%
	inner_join(totals2, by = c("STATE" = "state", "COUNTY" = "county"))
totals3_sf <- state.grid %>%
	inner_join(totals3, by = c("GEO_ID" = "GEO_ID"))
totals4_sf <- state.grid.finer %>%
	inner_join(totals4, by = c("GEO_ID" = "GEO_ID"))
totals5_sf <- state.grid.coarser %>%
	inner_join(totals5, by = c("GEO_ID" = "GEO_ID"))

plot(totals1_sf[,c("y_bar", "E_hat")])
plot(totals2_sf[,c("y_bar", "E_hat")])
plot(totals3_sf[,c("y_bar", "E_hat")])
plot(totals4_sf[,c("y_bar", "E_hat")])
plot(totals5_sf[,c("y_bar", "E_hat")])

plot(totals1_sf[,c("err")])
plot(totals2_sf[,c("err")])
plot(totals3_sf[,c("err")])
plot(totals4_sf[,c("err")])
plot(totals5_sf[,c("err")])

plot(totals1_sf$y_bar, totals1_sf$E_hat); abline(c(0,1))
plot(totals2_sf$y_bar, totals2_sf$E_hat); abline(c(0,1))
plot(totals3_sf$y_bar, totals3_sf$E_hat); abline(c(0,1))
plot(totals4_sf$y_bar, totals4_sf$E_hat); abline(c(0,1))
plot(totals5_sf$y_bar, totals5_sf$E_hat); abline(c(0,1))

save.image("results.Rdata")
