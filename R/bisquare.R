# This is a naive version. We can be more efficient by gridding.
# uu is a q x d matrix, where each row is a spatial point
# cc is an r x d matrix of spatial cutpoints
# gg is an m.t vector of temporal cutpoints
# Return an q x r x d array of bases
local.bisquare.basis.points <- function(uu, tt, w.s, w.t, cc, gg)
{
	if (!is.matrix(uu)) { uu <- matrix(uu, nrow = 1) }
	q <- nrow(uu)
	d <- ncol(uu)

	r <- nrow(cc)
	stopifnot(d == ncol(cc))

	m.t <- length(gg)
	Psi <- array(0, dim = c(q, r, m.t))

	for (i in 1:q) {
		u <- uu[i,]
		dist.u <- rowSums((matrix(u, r, d, byrow = TRUE) - cc)^2)
		dist.t <- (tt - gg)^2
		idx.u <- which(dist.u <= w.s)
		idx.t <- which(dist.t <= w.t)
		for (idx in idx.t) {
			Psi[i, idx.u, idx] <- (1 - dist.u[idx.u]/w.s^2 - dist.t[idx]/w.t^2)^2
		}
	}

	return(Psi)
}

# Now D should somehow be a list of n.D areas
# For each area, we should sample points uniformly within the area,
# call local.bisquare.basis.points, and take the sum.
# How do R authors commonly represent a spatial area? (E.g. as a polygon?)
local.bisquare.basis.areas <- function(D, tt, w.s, w.t, cc, gg) {
	stop("TBD!")
}
