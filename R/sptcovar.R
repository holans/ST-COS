# I'm not sure what's going on here; just ported Jon's code
sptcovar9 <- function(Q, M, S, n)
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

	# We don't need to construct the FullCovar in its entirety to do this operation:
	# SpSinvSp <- SpSinv %*% t(S)
	# K <- SpSinvSp %*% FullCovar %*% t(SpSinvSp)
	K <- matrix(0, r, r)
	C <- as.matrix(Q)
	for (i in 1:9) {
		logger("Working on block (%d, %d)\n", i, i)
		idx <- seq(n*(i-1)+1, n*i)
		A.i <- as.matrix(SpSinv %*% t(S[idx,]))
		K <- K + A.i %*% (C %*% t(A.i))
	}

	# for (i in 1:9) {
	# 	for (j in setdiff(1:9, 1:i)) {
	# 		logger("Working on block (%d, %d)\n", i, j)
	# 		postmult <- diag(1,n)
	# 		for (k in (i+1):j) {
	# 			postmult <- M %*% postmult
	# 		}
	# 		C <- as.matrix(postmult %*% (Q %*% t(postmult)))
	# 
	# 		idx.i <- seq(n*(i-1)+1, n*i)
	# 		idx.j <- seq(n*(j-1)+1, n*j)
	# 		A.i <- as.matrix(SpSinv %*% t(S[idx.i,]))
	# 		A.j <- as.matrix(SpSinv %*% t(S[idx.j,]))
	# 
	# 		# K <- K + A.i %*% (C %*% t(A.j))
	# 		# K <- K + A.j %*% (t(C) %*% t(A.i))
	# 		K <- K + 2 * A.i %*% (C %*% t(A.j))
	# 	}
	# }

	M.h <- diag(1,n)
	for (h in 1:8) {
		M.h <- M %*% M.h
		C <- as.matrix(M.h %*% (Q %*% t(M.h)))

		for (i in seq(1, 9-h)) {
			j <- i + h
			logger("Working on block (%d, %d)\n", i, j)

			idx.i <- seq(n*(i-1)+1, n*i)
			idx.j <- seq(n*(j-1)+1, n*j)
			A.i <- as.matrix(SpSinv %*% t(S[idx.i,]))
			A.j <- as.matrix(SpSinv %*% t(S[idx.j,]))
			K <- K + 2 * A.i %*% (C %*% t(A.j))
		}
	}

	return(K)
}
