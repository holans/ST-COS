# HMATRIX_OVERLAP.M 
# calculate area of overlap of data polygons and grid polygons
# to form a change of support matrix, H
#
#   D - structural array of length nd that has .X and .Y locations for 
#        data polygons
#   
#   G - structural array of length ng that has .X and .Y locations for 
#         grid polygons [NOTE: for school districts, use .LON and .LAT 
#         instead of .X and .Y
#  
#   denomflag: 1 - denom = 1
#              2 - D area in the denom
#              3 - G area in the denom
#
#   Dxyflag: 1 - use .X and .Y, else use .LON and .LAT
#   Gxyflag: 1 - use .X and .Y, else use .LON and .LAT
#
#   n_mc: number of Monte Carlo samples
#
#  modification of 7/25/2012 version for structural
#    arrays; 01/09/2013 ckw
# -------------------------------------------------
compute.overlap.v2 <- function(D, G, denomflag = 1, report.period = nrow(D) + 1)
{
	nd <- nrow(D)
	ng <- nrow(G)
	H <- Matrix(0, nd, ng)

	D.xmin <- numeric(nd)
	D.xmax <- numeric(nd)
	D.ymin <- numeric(nd)
	D.ymax <- numeric(nd)
	D.area <- numeric(nd)
	for (k in 1:nd) {
		dd <- D[k,]
		D.area[k] <- st_area(dd)
		D.xmin[k] <- st_bbox(dd)[1]
		D.ymin[k] <- st_bbox(dd)[2]
		D.xmax[k] <- st_bbox(dd)[3]
		D.ymax[k] <- st_bbox(dd)[4]
	}

	G.xmin <- numeric(ng)
	G.xmax <- numeric(ng)
	G.ymin <- numeric(ng)
	G.ymax <- numeric(ng)
	G.area <- numeric(ng)
	for (k in 1:ng) {
		gg <- G[k,]
		G.area[k] <- st_area(gg)
		G.xmin[k] <- st_bbox(gg)[1]
		G.ymin[k] <- st_bbox(gg)[2]
		G.xmax[k] <- st_bbox(gg)[3]
		G.ymax[k] <- st_bbox(gg)[4]
	}

	for (k in 1:nd) {
		if (k %% report.period == 0) {
			logger("Location %d of %d\n", k, nd)
		}
		dd <- D[k,]

		# Only consider the areas in G whose bounding box intersects with area dd
		idx1 <- which(D.xmin[k] <= G.xmin & G.xmin <= D.xmax[k])
		idx2 <- which(D.xmin[k] <= G.xmax & G.xmax <= D.xmax[k])
		idx3 <- which(D.ymin[k] <= G.ymin & G.ymin <= D.ymax[k])
		idx4 <- which(D.ymin[k] <= G.ymax & G.ymax <= D.ymax[k])
		idx <- sort(unique(c(idx1, idx2, idx3, idx4)))

		for (j in idx) {
			gg <- G[j,]
			pi <- st_intersection(dd, gg)
			DGoverlap <- max(as.numeric(st_area(pi)), 0)

			if (denomflag == 1) {
				H[k,j] <- DGoverlap
			} else if (denomflag == 2) {
				H[k,j] <- DGoverlap / D.area[k]
			} else if (denomflag == 3) {
				H[k,j] <- DGoverlap / G.area[j]
			} else {
				stop("denomflag must be 1, 2, or 3")
			}
		}
	}

	return(H)
}

