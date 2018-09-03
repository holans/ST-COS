compute_sp_basis_mc <- function(basis, domain, R, report.period = 100)
{
	n <- nrow(domain)
	r <- basis$get_dim()
	S <- Matrix(0, n, r)

	for (j in 1:n) {
		if (j %% report.period == 0) {
			logger("Computing basis for area %d of %d\n", j, n)
		}

		# Request a few more samples than we'll need, to prevent the loop in rArea.
		P <- rArea(R, domain[j,], blocksize = ceiling(1.2*R))
		S[j,] <- S[j,] + colSums(basis$compute(P[,1], P[,2]))
	}

	return(S / R)
}

compute_spt_basis_mc <- function(basis, domain, R, period, report.period = 100)
{
	n <- nrow(domain)
	r <- basis$get_dim()
	S <- Matrix(0, n, r)
	T <- length(period)

	for (j in 1:n) {
		if (j %% report.period == 0) {
			logger("Computing basis for area %d of %d\n", j, n)
		}

		# Drawing samples from an area is more time consuming than computing
		# basis function. Let's reuse samples over multiple lookbacks.
		# Request a few more samples than we'll need, to prevent the loop in rArea.
		P <- rArea(R, domain[j,], blocksize = ceiling(1.2*R))

		for (t in 1:T) {
			S[j,] <- S[j,] + colSums(basis$compute(P[,1], P[,2], period[t]))
		}
	}

	return( S / (R*T) )
}
