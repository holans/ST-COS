# Computes 
#   K <- SpSinvSp %*% FullCovar %*% t(SpSinvSp)
# where
#   SpSinvSp <- solve(t(S) %*% S) %*% t(S)
#   FullCovar is the autocovariance of VAR(1), the first lag_max blocks
# We don't need to construct FullCovar in its entirety to do this operation.
#
# K minimizes || FullCovar - S X S^T ||_Frob over symmetric psd matrices X
sptcovar.vectautoreg <- function(Qinv, M, S, lag_max)
{
	n <- nrow(Qinv)
	d <- nrow(S)
	r <- ncol(S)

	SpS <- t(S) %*% S
	if (rankMatrix(SpS) < r) {
		warning("The matrix (S' S) is rank-deficient. Consider reducing the dimension of S")
	}
	SpSinv <- ginv(as.matrix(SpS))

	# Get all the the autocovariances we'll need
	logger("About to call covVAR1\n")
	G <- covVAR1(M, Qinv, lag_max = lag_max)
	Gamma <- function(h) {
		if (h >= 0) {
			G[,,h+1]
		} else {
			t(G[,,-h+1])
		}
	}

	C <- matrix(0, r, r)
	for (i in 1:lag_max) {
		for (j in 1:lag_max) {
			logger("Computing block (%d,%d)\n", i, j)
			idx.i <- seq(n*(i-1)+1, n*i)
			idx.j <- seq(n*(j-1)+1, n*j)
			C <- C + t(S[idx.i,]) %*% Gamma(i-j) %*% S[idx.j,]
		}
	}

	logger("Finished computing K\n")
	return(SpSinv %*% C %*% SpSinv)
}

sptcovar.randwalk <- function(Qinv, M, S, lag_max)
{
	n <- nrow(Qinv)
	d <- nrow(S)
	r <- ncol(S)
	SpS <- t(S) %*% S

	if (rankMatrix(SpS) < r) {
		warning("The matrix (S' S) is rank-deficient. Consider reducing the dimension of S")
	}
	SpSinv <- ginv(as.matrix(SpS))

	C <- matrix(0, r, r)
	for (i in 1:lag_max) {
		for (j in 1:lag_max) {
			logger("Computing block (%d,%d)\n", i, j)
			idx.i <- seq(n*(i-1)+1, n*i)
			idx.j <- seq(n*(j-1)+1, n*j)
			B <- min(i,j)*Qinv
			C <- C + t(S[idx.i,]) %*% B %*% S[idx.j,]
		}
	}

	logger("Finished computing K\n")
	return(SpSinv %*% C %*% SpSinv)
}
