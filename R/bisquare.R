# This is a naive version. We can be more efficient by gridding.
# uu is a q x d matrix, where each row is a spatial point
# cc is an r x d matrix of spatial cutpoints
# gg is an m.t vector of temporal cutpoints
# Return an q x r x d array of bases
local.bisquare.basis.points <- function(uu, times, w.s, w.t, cc, gg)
{
	if (!is.matrix(uu)) { uu <- matrix(uu, nrow = 1) }
	q <- nrow(uu)
	d <- ncol(uu)
	r <- nrow(cc)
	stopifnot(d == ncol(cc))
	stopifnot(length(times) == 1)

	m.t <- length(gg)
	Psi <- array(0, dim = c(q, r, m.t))

	for (i in 1:q) {
		u <- uu[i,]
		dist.u <- sqrt(rowSums((matrix(u, r, d, byrow = TRUE) - cc)^2))
		dist.t <- sqrt((times - gg)^2)
		idx.u <- which(dist.u <= w.s)
		idx.t <- which(dist.t <= w.t)
		if (length(idx.u) > 0) {
			for (idx in idx.t) {
				Psi[i, idx.u, idx] <- (1 - dist.u[idx.u]^2/w.s^2 - dist.t[idx]^2/w.t^2)^2
			}
		}
	}

	return(Psi)
}

# For a single spatial area, sample points uniformly within the area,
# call local.bisquare.basis.points, and take the sum.
# How do R authors commonly represent a spatial area? (E.g. as a polygon?)
local.bisquare.basis.area <- function(A, times, w.s, w.t, cc, gg, n = 500)
{
	stopifnot(length(times) == 1)
	stopifnot(class(A) == "SpatialPolygonsDataFrame")

	lon.low <- A@bbox[1,1]
	lon.hi <- A@bbox[1,2]
	lat.low <- A@bbox[2,1]
	lat.hi <- A@bbox[2,2]

	u.lon <- runif(n, lon.low, lon.hi)
	u.lat <- runif(n, lat.low, lat.hi)
	u <- data.frame(Longitude = u.lon, Latitude = u.lat)
	coordinates(u) <- ~ Longitude + Latitude
	proj4string(u) <- proj4string(A)

	# Now reject the ones outside of the area
	idx <- which(!is.na(over(u, A)$GEO_ID))
	uu <- u@coords[idx,]

	# plot(A)
	# points(u[-idx], col = "black")
	# points(u[idx], col = "red")

	Psi <- local.bisquare.basis.points(uu, times, w.s, w.t, cc, gg)
	apply(Psi, c(2,3), mean)
}

# A domain D is a geography D and a time times
local.bisquare.basis.domain <- function(D, times, w.s, w.t, cc, gg, n = 500)
{
	stopifnot(length(times) == 1)
	n <- length(D)
	m.t <- length(gg)
	m.c <- nrow(cc)
	S <- Matrix(0, n, m.c)

	# TBD: Check this against Jon's code. I just made this up so far.

	for (k in 1:length(times)) {
		for (i in 1:n) {
			logger("Computing bisquare basis for time %f, area %d in domain\n", times[k], i)
			ret <- local.bisquare.basis.area(D[i,], times[k], w.s, w.t, cc, gg, n)
			S[i,] <- S[i,] + as.numeric(ret)
		}
	}

	return(S / m.t)
}
