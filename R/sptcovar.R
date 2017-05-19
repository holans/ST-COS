# Computes 
#   K <- SpSinvSp %*% FullCovar %*% t(SpSinvSp)
# where
#   SpSinvSp <- solve(t(S) %*% S) %*% t(S)
#   FullCovar is the autocovariance of VAR(1), the first lag_max blocks
# We don't need to construct FullCovar in its entirety to do this operation.
#
# K minimizes || FullCovar - S X S^T ||_Frob over symmetric psd matrices X
sptcovar <- function(Qinv, M, S, lag_max)
{
	n <- nrow(Qinv)
	d <- nrow(S)
	r <- ncol(S)
	SpS <- t(S) %*% S

	if (rankMatrix(SpS) < ncol(SpS)) {
		warning("The matrix (S' S) is rank-deficient. Consider reducing the dimension of S")
	}

	SpSinv <- ginv(SpS)

	# Get the autocovariance at lag 0
	# We can compute the other lags as we need them (without storing all of them)
	Gamma <- covVAR1(M, Qinv, lag_max = 0)
	C <- Gamma[,,1]
	K <- matrix(0, r, r)
	for (i in 1:lag_max) {
		idx <- seq(n*(i-1)+1, n*i)
		A.i <- as.matrix(SpSinv %*% t(S[idx,]))
		K <- K + A.i %*% (C %*% t(A.i))
	}

	M.h <- diag(1,n)
	for (h in 1:(lag_max-1)) {
		M.h <- M %*% M.h
		C <- as.matrix(M.h %*% (Gamma[,,1] %*% t(M.h)))

		for (i in seq(1, lag_max-h)) {
			j <- i + h
			idx.i <- seq(n*(i-1)+1, n*i)
			idx.j <- seq(n*(j-1)+1, n*j)
			A.i <- as.matrix(SpSinv %*% t(S[idx.i,]))
			A.j <- as.matrix(SpSinv %*% t(S[idx.j,]))
			B <- A.i %*% (C %*% t(A.j))
			K <- K + 2*B
		}
	}

	return(K)
}
