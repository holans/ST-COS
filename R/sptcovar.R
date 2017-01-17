# I'm not sure what' going on here; just ported Jon's code
make_full_model_sptcovar_9 <- function(Qinv, M, S, n1)
{
	FullCovar <- matrix(0, nrow(S), ncol(S))
	for (i in 1:9) {
		for (j in 1:9) {
			postmult <- diag(1,n1)
			if (i == j) {
				FullCovar(n1*(i-1)+1:n1*i, n1*(i-1)+1:n1*i) <- 0.5 %*% Qinv
			} else if (j > i) {
				for (k in (i+1):j) {
					postmult <- M %*% postmult
				}
				FullCovar[n1*(i-1)+1:n1*i, n1*(j-1)+1:n1*j] = (postmult %*% Qinv) %*% t(postmult)
		    }
		}
		print(i);
	}

	FullCovar <- t(FullCovar) + FullCovar
	SpS <- t(S) %*% S
	SpSinv <- ginv(SpS)
	SpSinvSp <- SpSinv %*% t(S)
	K <- (SpSinvSp %*% FullCovar) %*% t(SpSinvSp)
	return(K)
}
