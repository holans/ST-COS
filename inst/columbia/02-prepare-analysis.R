# -------------------------------------------
# ggforce is used to plot circles in ggplots
# fields is used for cover.design function
# sf, dplyr, and ggplot2 are used throughout
# -------------------------------------------
library(sf)
library(ggplot2)
library(ggforce)
library(dplyr)
library(fields)
library(stcos)

set.seed(1234)

# Load previously constructed source supports.
data(acs_sf)

# There are some NA values of DirectEst and DirectVar in the data.
# Let's create source supports from these objects with the NAs filtered out.
source_2013 = acs5_2013 %>% filter(!is.na(DirectEst) & !is.na(DirectVar))
source_2014 = acs5_2014 %>% filter(!is.na(DirectEst) & !is.na(DirectVar))
source_2015 = acs5_2015 %>% filter(!is.na(DirectEst) & !is.na(DirectVar))
source_2016 = acs5_2016 %>% filter(!is.na(DirectEst) & !is.na(DirectVar))
source_2017 = acs5_2017 %>% filter(!is.na(DirectEst) & !is.na(DirectVar))

# We plan to use the acs5_2017 geography as the fine-level support.
# For some areas in our planned fine-level domain, there is very little
# overlap area with anything in the NA-filtered source supports. We will
# drop these areas from the fine-level domain to avoid rank-deficiency
# issues with the H matrix.
U = rbind(
	overlap_matrix(source_2013, acs5_2017, proportion = FALSE),
	overlap_matrix(source_2014, acs5_2017, proportion = FALSE),
	overlap_matrix(source_2015, acs5_2017, proportion = FALSE),
	overlap_matrix(source_2016, acs5_2017, proportion = FALSE),
	overlap_matrix(source_2017, acs5_2017, proportion = FALSE)
)
dom_fine = acs5_2017 %>%
	mutate(keep = (colSums(U) >= 10)) %>%
	filter(keep == TRUE) %>%
	select(-c("DirectEst", "DirectMOE", "DirectVar", "keep"))
n = nrow(dom_fine)

# A quick plot of the fine-level domain.
ggplot(dom_fine) +
	geom_sf(colour = "black", size = 0.05) +
	ggtitle("Boone County, Missouri") +
	theme_bw()

# A quick plot of one of the source supports.
ggplot(acs5_2017) +
	geom_sf(colour = "black", size = 0.05, aes(fill = DirectEst)) +
	scale_fill_distiller("DirectEst", palette = "RdYlBu") +
	theme_bw()

# Load previously constructed target supports.
# Make sure to transform to the same projection as the fine-level support.
data(columbia_neighbs)
neighbs = columbia_neighbs %>% st_transform(crs = st_crs(dom_fine))

# A quick plot of the neighborhoods.
ggplot(neighbs) +
	geom_sf(colour = "black", size = 0.05) +
	theme_bw()

# Compute overlap matrix H.
H = rbind(
	overlap_matrix(source_2013, dom_fine),
	overlap_matrix(source_2014, dom_fine),
	overlap_matrix(source_2015, dom_fine),
	overlap_matrix(source_2016, dom_fine),
	overlap_matrix(source_2017, dom_fine)
)
N = nrow(H)

# Select spatial knots via space-filling design. Start with uniformly
# drawn points from the fine-level geography and use `cover.design`
# to select a subset of points to use as knots.
u = st_sample(dom_fine, size = 2000)
P = matrix(unlist(u), length(u), 2, byrow = TRUE)
out = fields::cover.design(P, 200)
knots_sp = out$design

# Evenly spaced points can also be achieved with hexagonal sampling in the sf
# package.
u = st_sample(dom_fine, 200, type = "hexagonal")
knots_sp_alt = matrix(unlist(u), length(u), 2, byrow = TRUE)

# Select temporal knots to be evenly spaced over the years relevant
# to the source support years.
knots_t = seq(2009, 2017, by = 0.5)

# Use Cartesian join of spatial and temporal knots to obtain spatio-temporal knots.
knots = merge(knots_sp, knots_t)

# Choose a spatial radius for the basis. The following choice puts
# on a scale appropriate to the projection used in the fine-level
# domain. Adjust w_s to scale it up or down.
ws_tilde = 1
D = dist(knots_sp)
w_s = ws_tilde * quantile(D[D > 0], prob = 0.05, type = 1)

# Create an ArealSpaceTimeBisquareBasis object with our knot points.
bs_spt = ArealSpaceTimeBisquareBasis$new(knots[,1], knots[,2], knots[,3],
	w_s = as.numeric(w_s), w_t = 1, mc_reps = 500)

# Here is a plot of the spatial knots. We plot a circle around one of the
# points to illustrate the choice of the `w_s` argument.
rad = bs_spt$get_basis_spt()$get_ws()
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
ggsave(file = "spatial-knots-cover.pdf", g, height = 5, width = 4)

# Compare the above to the points from hexogonal sampling
ws_tilde = 1
D = dist(knots_sp_alt)
w_s_alt = ws_tilde * quantile(D[D > 0], prob = 0.05, type = 1)
knots_sp_dat = data.frame(x = knots_sp_alt[,1], y = knots_sp_alt[,2], r = as.numeric(w_s_alt))
g = ggplot(dom_fine) +
	geom_sf(colour = "black", size = 0.25, fill = NA) +
	geom_point(data = knots_sp_dat, aes(x, y), lwd = 1, col = "red") +
	geom_point(data = knots_sp_dat[1,], aes(x, y), lwd = 3, col = "blue") +
	geom_circle(data = knots_sp_dat[1,], aes(x0=x, y0=y, r=r), fill = NA,
		lwd = 0.5, col = "blue") +
	labs(x = "", y = "") +
	theme_bw()
print(g)
ggsave(file = "spatial-knots-hex.pdf", g, height = 5, width = 4)

# Compute basis function matrix S for source supports.
# This can become time-consuming when there are a large number of areal units,
# or many mc_reps are requested.
S_full = rbind(
	bs_spt$compute(source_2013, 2009:2013),
	bs_spt$compute(source_2014, 2010:2014),
	bs_spt$compute(source_2015, 2011:2015),
	bs_spt$compute(source_2016, 2012:2016),
	bs_spt$compute(source_2017, 2013:2017)
)

# Compute basis function on fine-level domain. This may be required for
# computing K. As with the previous step, this can be time-consuming in
# larger problems.
S_fine_full = rbind(
	bs_spt$compute(dom_fine, 2009),
	bs_spt$compute(dom_fine, 2010),
	bs_spt$compute(dom_fine, 2011),
	bs_spt$compute(dom_fine, 2012),
	bs_spt$compute(dom_fine, 2013),
	bs_spt$compute(dom_fine, 2014),
	bs_spt$compute(dom_fine, 2015),
	bs_spt$compute(dom_fine, 2016),
	bs_spt$compute(dom_fine, 2017)
)

# Extract the direct estimates and variance estimates
z = c(source_2013$DirectEst, source_2014$DirectEst, source_2015$DirectEst,
	source_2016$DirectEst, source_2017$DirectEst)
v = c(source_2013$DirectVar, source_2014$DirectVar, source_2015$DirectVar,
	source_2016$DirectVar, source_2017$DirectVar)

# Do a PCA reduction on S to reduce its dimension
eig = eigen(t(S_full) %*% S_full)
idx_S = which(cumsum(eig$values) / sum(eig$values) < 0.65)
Tx_S = eig$vectors[,idx_S]

# Plot the proportion of variation captured by our selection of PCA components.
pdf("pca-reduction.pdf", width = 6, height = 5)
eigprops = cumsum(eig$values) / sum(eig$values)
plot(eigprops[1:100], xlab = "Dimension", ylab = "Proportion of Variation")
abline(v = max(idx_S), lty = 2)
abline(h = eigprops[max(idx_S)], lty = 2)
dev.off()

# Apply the same reduction to S_full and S_fine_full.
S = S_full %*% Tx_S
S_fine = S_fine_full %*% Tx_S
r = ncol(S)

# Get adjacency matrix for fine-level support and use it to compute
# a covariance structure (the nonsingular type) of a CAR process.
# This may be required for computing K.
A = adjacency_matrix(dom_fine)
aa = rowSums(A) + (rowSums(A) == 0)
W = 1/aa * A
tau = 0.9
Q = Diagonal(n,1) - tau*W
Qinv = solve(Q)

# Try all of the covariance structures for K, then pick one for the analysis.
K_list = list()

# Compute K by Random Walk method, which assumes CAR structure for obs
# within time points and random walk across time points.
K_list[["randwalk"]] = cov_approx_randwalk(Qinv, S_fine)

# Compute K by spatial-only method, which assumes CAR structure for obs
# within time points and independence across time points.
K_list[["blockdiag"]] = cov_approx_blockdiag(Qinv, S_fine)

# Assume K is the identity matrix.
K_list[["identity"]] = Diagonal(n = r)

K = K_list[["randwalk"]]
Kinv = solve(K)

# Standardize observations before running MCMC.
z_scaled = (z - mean(z)) / sd(z)
v_scaled = v / var(z)
