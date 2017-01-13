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

compute.overlap <- function(D, G, denomflag, Dxyflag, Gxyflag, n_mc, report.period = 10)
{
	nd <- nrow(D@data)
	ng <- nrow(G@data)
	H <- Matrix(0, nd, ng)

	D.areas <- numeric(nd)
	D.xmin <- numeric(nd)
	D.xmax <- numeric(nd)
	D.ymin <- numeric(nd)
	D.ymax <- numeric(nd)
	for (k in 1:nd) {
		dd <- D[k,]
		D.areas[k] <- gArea(dd)
		D.xmin[k] <- dd@bbox[1,1]
		D.xmax[k] <- dd@bbox[1,2]
		D.ymin[k] <- dd@bbox[2,1]
		D.ymax[k] <- dd@bbox[2,2]
	}

	G.areas <- numeric(ng)
	G.xmin <- numeric(ng)
	G.xmax <- numeric(ng)
	G.ymin <- numeric(ng)
	G.ymax <- numeric(ng)
	for (k in 1:ng) {
		gg <- G[k,]
		G.areas[k] <- gArea(gg)
		G.xmin[k] <- gg@bbox[1,1]
		G.xmax[k] <- gg@bbox[1,2]
		G.ymin[k] <- gg@bbox[2,1]
		G.ymax[k] <- gg@bbox[2,2]		
	}

	for (k in 1:nd) {
		if (k %% report.period == 0) {
			logger("Location %d of %d\n", k, nd)
		}

		# if (Dxyflag == 1) {
		# 	dx <- D[k]$X
		# 	dy <- D[k]$Y
		# } else {
		# 	dx <- D[k]$LON
		# 	dy <- D[k]$LAT
		# }
		
		dd <- D[k,]
		Darea <- gArea(dd)

		#dd.xmin <- dd@bbox[1,1]
		#dd.xmax <- dd@bbox[1,2]
		#dd.ymin <- dd@bbox[2,1]
		#dd.ymax <- dd@bbox[2,2]

		# Only consider the areas in G whose bounding box intersects with area dd
		b1 <- which(D.xmin[k] <= G.xmin & G.xmin <= D.xmax[k])
		b2 <- which(D.xmin[k] <= G.xmax & G.xmax <= D.xmax[k])
		b3 <- which(D.ymin[k] <= G.ymin & G.ymin <= D.ymax[k])
		b4 <- which(D.ymin[k] <= G.ymax & G.ymax <= D.ymax[k])
		idx <- sort(unique(c(b1, b2, b3, b4)))
		
		for (j in idx) {
			# if (Gxyflag == 1) {
			# 	gx <- G[j]$X
			# 	gy <- G[j]$Y
			# } else {
			# 	gx <- G[j]$LON
			# 	gy <- G[j]$LAT
			# }
			
			gg <- G[j,]
			#Garea <- gArea(gg)

			#gg.xmin <- gg@bbox[1,1]
			#gg.xmax <- gg@bbox[1,2]
			#gg.ymin <- gg@bbox[2,1]
			#gg.ymax <- gg@bbox[2,2]

			int <- gIntersection(dd, gg)
			if (length(int) == 0) {
				DGoverlap <- 0
			} else {
				DGoverlap <- gArea(int)
			}				

			# res <- MCAREA(dx, dy, gx, gy, n_mc)
			# Darea <- res$Darea
			# Garea <- res$Garea
			# DGoverlap <- res$DGoverlap

			if (denomflag == 1) {
				# 1 in denom
				H[k,j] <- DGoverlap
			} else if (denomflag == 2) {
				# D area in denom
				H[k,j] <- DGoverlap / Darea
			} else if (denomflag == 3) {
				# G area in denom
				H[k,j] <- DGoverlap / Garea
			} else {
				stop("denomflag must be 1, 2, or 3")
			}
		}
	}

	return(H)
}

