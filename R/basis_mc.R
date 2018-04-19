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

		# This block is just for debugging, and should be taken out soon
		if (FALSE) {
			browser()
			plot(dom.fine[,1], col = NA)
			points(P[,], pch = ".")
			points(t(P[1,]), pch = 18, col = "red", cex = 3)
			points(knots, pch = 16, col = "blue", cex = 3)

			t <- 1
			temp <- basis$compute(P[,1], P[,2], period[t])
			idx <- which(temp[1,] > 0)

			idx.counties <- which(dom.fine$NAME %in% c("Adams", "York", "Cumberland", "Franklin"))

			t2 <- cbind(knots[idx,], round(temp[1,idx], 6))
			plot(dom.fine[idx.counties,1], col = NA)
			points(t(P[1,]), pch = 18, col = "red", cex = 3)
			points(knots, pch = 16, col = "blue", cex = 3)
			points(t2[5,1:2], pch = 16, col = t2[,3], cex = 2)

			norm1 <- (P[1,1] - knots[idx[35],1])^2
			norm2 <- (P[1,2] - knots[idx[35],2])^2
			norm3 <- (period[t] - knots[idx[35],3])^2
			(norm1 + norm2) < basis$get_rl()^2
			norm3 <= basis$get_wt()^2
			(2 - (norm1 + norm2) / basis$get_rl()^2 - norm3 / basis$get_wt()^2)^2

			norm1 <- (P[1,1] - knots[idx[8],1])^2
			norm2 <- (P[1,2] - knots[idx[8],2])^2
			norm3 <- (period[t] - knots[idx[8],3])^2
			(norm1 + norm2) < basis$get_rl()^2
			norm3 <= basis$get_wt()^2
			(2 - (norm1 + norm2) / basis$get_rl()^2 - norm3 / basis$get_wt()^2)^2

			basis.new <- SpaceTimeBisquareBasis$new(knots[,1], knots[,2], knots[,3], w.s = radius, w.t = 1)
			round(basis.new$compute(P[1,1], P[1,2], period[t])[idx], 6)
		}

		for (t in 1:T) {
			S[j,] <- S[j,] + colSums(basis$compute(P[,1], P[,2], period[t]))
		}
	}

	return( S / (R*T) )
}
