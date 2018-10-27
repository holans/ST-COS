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
	SpSinv <- pinv(as.matrix(SpS))

	# Get all the the autocovariances we'll need
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
			idx.i <- seq(n*(i-1)+1, n*i)
			idx.j <- seq(n*(j-1)+1, n*j)
			C <- C + t(S[idx.i,]) %*% Gamma(i-j) %*% S[idx.j,]
		}
	}

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
	SpSinv <- pinv(as.matrix(SpS))

	C <- matrix(0, r, r)
	for (i in 1:lag_max) {
		for (j in 1:lag_max) {
			idx.i <- seq(n*(i-1)+1, n*i)
			idx.j <- seq(n*(j-1)+1, n*j)
			B <- min(i,j)*Qinv
			C <- C + t(S[idx.i,]) %*% B %*% S[idx.j,]
		}
	}

	return(SpSinv %*% C %*% SpSinv)
}

sptcovar.indep <- function(Qinv, S, lag_max)
{
	n <- nrow(Qinv)
	d <- nrow(S)
	r <- ncol(S)
	SpS <- t(S) %*% S
	
	if (rankMatrix(SpS) < r) {
		warning("The matrix (S' S) is rank-deficient. Consider reducing the dimension of S")
	}
	SpSinv <- pinv(as.matrix(SpS))
	
	C <- matrix(0, r, r)
	for (i in 1:lag_max) {
		idx <- seq(n*(i-1)+1, n*i)
		C <- C + t(S[idx,]) %*% Qinv %*% S[idx,]
	}

	return(SpSinv %*% C %*% SpSinv)
}
