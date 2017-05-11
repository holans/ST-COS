library(MTS)
library(Matrix)

m <- 3
Sigma <- c(1,2,1) %o% c(1,2,1) + c(1,1,1) %o% c(1,1,1) + c(3,2,1) %o% c(3,2,1)
M <- 7/16 * diag(3) + 1/16
Gamma <- make_full_model_sptcovar_jon(Sigma, M, lag_max = 9)

h <- 0
Gamma[1:m, 1:m + h*m]

h <- 1
Gamma[1:m, 1:m + h*m]

VARMAcov(Phi = M, Sigma = Sigma, lag = 1)
covVAR1(M, Sigma, lag_max = 1)

# ----- Now try to compute inv(S' S) S' F S inv(S' S)
S <- diag(c(1,2,3), m*3, 3)

Gamma <- covVAR1(M, Sigma, lag_max = 2)
F <- matrix(NA, m*3, m*3)
F[1:m,1:m] <- Gamma[,,1]
F[1:m + m,1:m + m] <- Gamma[,,1]
F[1:m + 2*m,1:m + 2*m] <- Gamma[,,1]
F[1:m,1:m + m] <- Gamma[,,2]
F[1:m,1:m + 2*m] <- Gamma[,,3]
F[1:m + m,1:m] <- t(Gamma[,,2])
F[1:m + m,1:m + 2*m] <- Gamma[,,2]
F[1:m + 2*m,1:m] <- t(Gamma[,,3])
F[1:m + 2*m,1:m + m] <- t(Gamma[,,2])

Sc <- solve(t(S) %*% S, t(S))
Sc %*% F %*% t(Sc)

sptcovar(Sigma, M, S, lag_max = 2)
