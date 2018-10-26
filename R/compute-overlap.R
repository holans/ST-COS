# Compute the overlap between two sets of areas (sf objects): D and G.
# geo.name.D and geo.name.G should be strings that identify the names of
# columns that uniquely identify the geography, for D and G respectively.
compute.overlap <- function(D, G, geo.name.D, geo.name.G)
{
	nd <- nrow(D)
	ng <- nrow(G)

	# Suppress unhelpful warnings when calling st_intersection
	st_agr(D) <- "constant"
	st_agr(G) <- "constant"

	# Compute the overlaps
	INT <- st_intersection(D, G)
	
	# If geo.name.D == geo.name.G, the column of INT corresponding to geo.name.G
	# will have a ".1" appended to it.
	if (geo.name.D == geo.name.G) {
		geo.name.G.INT <- sprintf("%s.1", geo.name.G)
		stopifnot(geo.name.G.INT %in% colnames(INT))
	} else {
		geo.name.G.INT <- geo.name.G
	}
	DF <- data.frame(
		D_GEO_ID = INT[[geo.name.D]],
		G_GEO_ID = INT[[geo.name.G.INT]],
		area.overlap = as.numeric(st_area(INT))
	)

	# Translate from GEO_IDs to row numbers of D and G
	ID.D <- data.frame(row_D = 1:nd, D_GEO_ID = D[[geo.name.D]])
	ID.G <- data.frame(row_G = 1:ng, G_GEO_ID = G[[geo.name.G]])

	# Create a sparse matrix with the overlaps, using the indices
	DF2 <- merge(ID.D, DF, by = "D_GEO_ID")
	DF3 <- merge(ID.G, DF2, by = "G_GEO_ID")
	sparseMatrix(i = DF3$row_D, j = DF3$row_G, x = DF3$area.overlap, dims = c(nd, ng))
}

