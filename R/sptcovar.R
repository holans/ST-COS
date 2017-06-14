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

	# Get the autocovariance at lag 0
	# We can compute the other lags as we need them (without storing all of them)
	C <- matrix(0, r, r)

	logger("About to call covVAR1\n")
	Gamma <- covVAR1(M, Qinv, lag_max = 0)

	B <- Gamma[,,1]
	for (i in 1:lag_max) {
		logger("Computing block (%d,%d)\n", i, i)
		idx <- seq(n*(i-1)+1, n*i)
		S.i <- S[idx,]
		C <- C + t(S.i) %*% B %*% S.i
	}

	M.h <- diag(1,n)
	for (h in 1:(lag_max-1)) {
		M.h <- M %*% M.h
		B <- as.matrix(M.h %*% Gamma[,,1] %*% t(M.h))

		for (i in seq(1, lag_max-h)) {
			logger("Computing block (%d,%d)\n", i, h)
			j <- i + h
			idx.i <- seq(n*(i-1)+1, n*i)
			idx.j <- seq(n*(j-1)+1, n*j)
			S.i <- S[idx.i,]
			S.j <- S[idx.j,]
			C <- C + 2*(t(S.i) %*% B %*% S.j)
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
