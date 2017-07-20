compute_sp_basis_mc <- function(basis, domain, R, report.period = 100)
{
	n <- nrow(domain)
	r <- basis$get_dim()
	S <- Matrix(0, n, r)
	areas <- as.numeric(st_area(domain))

	for (j in 1:n) {
		if (j %% report.period == 0) {
			logger("Computing basis for area %d of %d\n", j, n)
		}

		for (t in 1:T) {
			P <- rArea(R, domain[j,])
			S[j,] <- S[j,] + colSums(basis$compute(P[,1], P[,2]))
		}

		S[j,] <- S[j,] / areas[j]
	}

	return(S / R)
}

compute_spt_basis_mc <- function(basis, domain, R, period, report.period = 100)
{
	n <- nrow(domain)
	r <- basis$get_dim()
	S <- Matrix(0, n, r)
	T <- length(period)
	areas <- as.numeric(st_area(domain))

	for (j in 1:n) {
		if (j %% report.period == 0) {
			logger("Computing basis for area %d of %d\n", j, n)
		}

		for (t in 1:T) {
			P <- rArea(R, domain[j,])
			S[j,] <- S[j,] + colSums(basis$compute(P[,1], P[,2], period[t]))
		}

		S[j,] <- S[j,] / areas[j]
	}

	return( S / (R*T) )
}
