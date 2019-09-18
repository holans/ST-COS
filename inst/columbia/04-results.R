library(gridExtra)
library(ggrepel)

# The objective of our analysis - predictions on the four target neighborhoods.
print(nb_out)

# Plot maps of direct- and model-based 2017 5-year estimates.
lim_est = range(acs5_2017_out$DirectEst, acs5_2017_out$E_mean)
g = ggplot(acs5_2017_out) +
	geom_sf(colour = "black", size = 0.05, aes(fill = DirectEst)) +
	ggtitle("Median Household Income\nfor Boone County",
		subtitle = "2017 5yr ACS Direct Estimates") +
	scale_fill_distiller("DirectEst", palette = "RdYlBu", limits = lim_est) +
	theme_bw()
h = ggplot(acs5_2017_out) +
	geom_sf(colour = "black", size = 0.05, aes(fill = E_mean)) +
	ggtitle("Median Household Income\nfor Boone County",
		subtitle = "2017 5yr Model Estimates") +
	scale_fill_distiller("E_mean", palette = "RdYlBu", limits = lim_est) +
	theme_bw()
k = grid.arrange(g, h, ncol = 2)
ggsave("compare2017-direct.pdf", g, width = 7, height = 7)
ggsave("compare2017-model.pdf", h, width = 7, height = 7)

# Scatter plots comparing direct and model-based 5-year estimates for
# 2013, ..., 2017.
scatter_list = list()
years = 2013:2017
for (idx in 1:length(years)) {
	year = years[idx]
	obj = get(sprintf("acs5_%d_out", year))
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

ggsave("compare2014-scatter.pdf", scatter_list[[2]], width = 5, height = 5)
ggsave("compare2017-scatter.pdf", scatter_list[[5]], width = 5, height = 5)

# Plot neighborhood areas (target supports) among ACS 5-year direct estimates.
# This gives a sense of whether the model-based esimtates are reasonable.
# The remainder of the code in this file is working toward this plot.

idx_missing2017 = which(is.na(acs5_2017_out$DirectEst))
Central = nb_out[1,]
East = nb_out[2,]
North = nb_out[3,]
Paris = nb_out[4,]
Missing1 = acs5_2017_out[idx_missing2017[1],]
Missing2 = acs5_2017_out[idx_missing2017[2],]
Missing3 = acs5_2017_out[idx_missing2017[3],]
Missing4 = acs5_2017_out[idx_missing2017[4],]

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

g = ggplot(acs5_2017_out) +
	geom_sf(colour = "black", size = 0.05, aes(fill = DirectEst)) +
	ggtitle("Median Household Income for Boone County",
		subtitle = "ACS 2017 5yr Direct Estimates") +
	scale_fill_distiller("DirectEst", palette = "RdYlBu", limits = lim_est) +
	geom_sf(data = nb_out, fill = "black") +
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
ggsave("areas-of-interest-map.pdf", g)
