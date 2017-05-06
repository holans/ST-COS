# Computes 
#   K <- SpSinvSp %*% FullCovar %*% t(SpSinvSp)
# where
#   SpSinvSp <- solve(t(S) %*% S) %*% t(S)
#   FullCovar is the autocovariance of VAR(1), the first lag_max blocks
# We don't need to construct FullCovar in its entirety to do this operation.
#
# K minimizes || FullCovar - S X S^T ||_Frob over symmetric psd matrices X
sptcovar <- function(Q, M, S, n, lag_max)
{
	# d <- nrow(S)
	# FullCovar <- Matrix(0, d, d)
	# for (i in 1:9) {
	# 	for (j in 1:9) {
	# 		logger("Working on block (%d, %d) of FullCovar\n", i, j)
	# 		postmult <- diag(1,n)
	# 		if (i == j) {
	# 			idx <- seq(n*(i-1)+1, n*i)
	# 			FullCovar[idx, idx] <- 0.5 * Q
	# 		} else if (j > i) {
	# 			for (k in (i+1):j) {
	# 				postmult <- M %*% postmult
	# 			}
	# 			idx.i <- seq(n*(i-1)+1, n*i)
	# 			idx.j <- seq(n*(j-1)+1, n*j)
	# 			FullCovar[idx.i, idx.j] <- postmult %*% (Q %*% t(postmult))
	# 	    }
	# 	}
	# }
	# 
	# FullCovar <- t(FullCovar) + FullCovar
	# SpS <- t(S) %*% S
	# SpSinv <- solve(SpS)
	# SpSinvSp <- SpSinv %*% t(S)
	# K <- (SpSinvSp %*% FullCovar) %*% t(SpSinvSp)

	d <- nrow(S)
	r <- ncol(S)
	SpS <- t(S) %*% S
	SpSinv <- solve(SpS)

	K <- matrix(0, r, r)
	C <- as.matrix(Q)
	for (i in 1:lag_max) {
		logger("Working on block (%d, %d) with lag %d\n", i, i, 0)
		idx <- seq(n*(i-1)+1, n*i)
		logger("\tidx = [%d,%d]\n", min(idx), max(idx))
		A.i <- as.matrix(SpSinv %*% t(S[idx,]))
		K <- K + A.i %*% (C %*% t(A.i))
	}

	M.h <- diag(1,n)
	for (h in 1:(lag_max-1)) {
		M.h <- M %*% M.h
		C <- as.matrix(M.h %*% (Q %*% t(M.h)))

		for (i in seq(1, lag_max-h)) {
			j <- i + h
			logger("Working on block (%d, %d) with lag %d\n", i, j, h)
			idx.i <- seq(n*(i-1)+1, n*i)
			idx.j <- seq(n*(j-1)+1, n*j)
			logger("\tidx.i = [%d,%d], idx.j = [%d,%d]\n", min(idx.i), max(idx.i), min(idx.j), max(idx.j))
			A.i <- as.matrix(SpSinv %*% t(S[idx.i,]))
			A.j <- as.matrix(SpSinv %*% t(S[idx.j,]))
			B <- A.i %*% (C %*% t(A.j))
			K <- K + 2*B
		}
	}

	return(K)
}
