compute.overlap.v3 <- function(D, G)
{
	nd <- nrow(D)
	ng <- nrow(G)

	# Compute the overlaps
	INT <- st_intersection(D, G)
	DF <- data.frame(
		D_GEO_ID = INT$GEO_ID,
		G_GEO_ID = INT$GEO_ID.1,
		area.overlap = as.numeric(st_area(INT))
	)

	# Translate from GEO_IDs to row numbers of D and G
	ID.D <- data.frame(row_D = 1:nd, D_GEO_ID = D$GEO_ID)
	ID.G <- data.frame(row_G = 1:ng, G_GEO_ID = G$GEO_ID)

	# Create a sparse matrix with the overlaps, using the indices
	DF2 <- merge(ID.D, DF, by = "D_GEO_ID")
	DF3 <- merge(ID.G, DF2, by = "G_GEO_ID")
	sparseMatrix(i = DF3$row_D, j = DF3$row_G, x = DF3$area.overlap, dims = c(nd, ng))
}

