# Sigma is covariance of error term
# M is AR(1) matrix
covVAR1 <- function(A, Sigma, lag_max)
{
	m <- nrow(Sigma)
	N <- m * (lag_max+1)
	Gamma <- array(NA, dim = c(m, m, lag_max+1))

	# Gamma0_vec <- solve(Diagonal(m^2,1) - (A %x% A), matrix(Sigma, m^2, 1))
	browser()
	out <- solve_Gamma0(as.matrix(A), as.matrix(Sigma))
	Gamma[,,1] <- matrix(out$x, m, m)

	for (h in seq_len(lag_max)) {
		Gamma[,,h+1] <- A %*% Gamma[,,h]
	}

	return(Gamma)
}

make_full_model_sptcovar_jon <- function(Sigma, M, lag_max)
{
	m <- nrow(Sigma)
	N <- m * (lag_max+1)
	K <- matrix(0, N, N)
	for (i in 1:lag_max) {
		for (j in 1:lag_max) {
			postmult = diag(m)
			if (i == j) {
				idx1 <- seq(m*(i-1)+1, m*i)
				idx2 <- seq(m*(i-1)+1, m*i)
				K[idx1, idx2] = 0.5 * Sigma
			}
			else if (j > i) {
				for (k in (i+1):j) {
					postmult = M %*% postmult 
				}
				idx1 <- seq(m*(i-1)+1, m*i)
				idx2 <- seq(m*(j-1)+1, m*j)
				K[idx1, idx2] = postmult %*% Sigma %*% t(postmult)
			}
		}
	}

	t(K) + K
}
