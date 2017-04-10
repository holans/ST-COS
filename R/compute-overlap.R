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
compute.overlap <- function(D, G, denomflag = 1, report.period = nrow(D@data) + 1)
{
	nd <- nrow(D@data)
	ng <- nrow(G@data)
	H <- Matrix(0, nd, ng)

	D.xmin <- numeric(nd)
	D.xmax <- numeric(nd)
	D.ymin <- numeric(nd)
	D.ymax <- numeric(nd)
	D.area <- numeric(nd)
	for (k in 1:nd) {
		dd <- D[k,]
		D.area[k] <- gArea(dd)
		D.xmin[k] <- dd@bbox[1,1]
		D.xmax[k] <- dd@bbox[1,2]
		D.ymin[k] <- dd@bbox[2,1]
		D.ymax[k] <- dd@bbox[2,2]
	}

	G.xmin <- numeric(ng)
	G.xmax <- numeric(ng)
	G.ymin <- numeric(ng)
	G.ymax <- numeric(ng)
	G.area <- numeric(ng)
	for (k in 1:ng) {
		gg <- G[k,]
		G.area[k] <- gArea(gg)
		G.xmin[k] <- gg@bbox[1,1]
		G.xmax[k] <- gg@bbox[1,2]
		G.ymin[k] <- gg@bbox[2,1]
		G.ymax[k] <- gg@bbox[2,2]		
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

			int <- gIntersection(dd, gg)
			DGoverlap <- 0
			if (length(int) > 0) {
				DGoverlap <- gArea(int)
			}

			# We have to be vigilant about cleaning up after gIntersection,
			# or memory use will accumulate quickly.
			# rm(int); gc()

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

