# Keep drawing until we get n draws
# Also, might want to consider just picking points on a grid
sample_Ds <- function(n, A)
{
	stopifnot(class(A) == "SpatialPolygonsDataFrame")

	lon.low <- A@bbox[1,1]
	lon.hi <- A@bbox[1,2]
	lat.low <- A@bbox[2,1]
	lat.hi <- A@bbox[2,2]

	uu <- matrix(NA, 0, 2)
	n.cur <- 0

	while (n.cur < n) {
		u.lon <- runif(n, lon.low, lon.hi)
		u.lat <- runif(n, lat.low, lat.hi)
		u <- data.frame(Longitude = u.lon, Latitude = u.lat)
		coordinates(u) <- ~ Longitude + Latitude
		proj4string(u) <- proj4string(A)

		# Only accept the ones inside the target area
		idx <- which(!is.na(over(u, A)[,1]))
		uu <- rbind(uu, u@coords[idx,])
		n.cur <- nrow(uu)
	}

	# plot(A)
	# points(uu, col = "red")
	return(uu[1:n,])
}

ArealBi2 <- function(D, times, level1, level2 = matrix(NA, 0, 0),
	level3 = matrix(NA, 0, 0), B, srad, trad, report.period = 1)
{
	r <- nrow(level1) + nrow(level2) + nrow(level3)
	n <- length(D)
	S <- Matrix(0, n, r)
	T <- length(times)

	# s.lon and s.lat are n x B matrices
	# Note that we don't need to generate different draws over time for the same space
	s.lon <- matrix(NA, n, B)
	s.lat <- matrix(NA, n, B)

	for (i in 1:n) {
		p <- sample_Ds(B, D[i])
		s.lon[i,] <- p[,1]
		s.lat[i,] <- p[,2]
	}

	for (b in 1:B) {
		if (b %% report.period == 0) { logger("Iteration %d\n", b) }
		# TBD: I switched the order of the coordinates here. Was that a mistake?
		P <- cbind(s.lon[,b], s.lat[,b])
		for (t in 1:T) {
			Sq <- Create_S(P, times[t], level1, level2, level3, srad, trad)
			S <- S + Sq
		}
	}

	S / (B * T)
}

Create_S <- function(crd, times, level1, level2, level3, w, tw)
{
	if (nrow(level2) > 0 || nrow(level3) > 0) {
		stop("Levels 2 and 3 not supported in this implementation")
	}

	if (nrow(level1) == 0) {
		return(matrix(0, nrow(crd), 0))
	}

	# AMR: Not sure why Jon sorted this way, but we need to do it to get
	# similar results
	idx <- order(level1[,2])
	level1 <- level1[idx,]

	# TBD: This is kind of wasteful to be recomputing over and over
	G <- as.matrix(dist(level1))
	G[lower.tri(G)] <- NA
	rl <- w * quantile(G[G > 0], 0.05, na.rm = TRUE)
	S <- matrix(0, nrow(crd), nrow(level1))

	for (j in 1:nrow(crd)) {
		h <- sqrt((level1[,1] - crd[j,1])^2 + (level1[,2] - crd[j,2])^2)
		ht <- sqrt((level1[,3] - times)^2)
		s <- (1 - (h / rl)^2 - (ht / tw)^2)^2
		# s[h > rl | ht > tw] <- 0		# TBD: Shouldn't we be doing this??
		s[h > rl] <- 0
		S[j,] <- s
	}

	return(S)
}
