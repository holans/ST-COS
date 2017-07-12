# draw_spt_basis_mc_OLD <- function(R, domain, period_len, report.period = 100)
# {
# 	n <- nrow(domain)
# 	T <- period_len
# 	s1 <- array(NA, dim = c(T, n, R))
# 	s2 <- array(NA, dim = c(T, n, R))	
# 
# 	for (j in 1:n) {
# 		if (j %% report.period == 0) {
# 			logger("Drawing points for area %d of %d\n", j, n)
# 		}
# 		for (t in 1:T) {
# 			P <- rArea(R, domain[j,])
# 			s1[t,j,] <- P[,1]
# 			s2[t,j,] <- P[,2]
# 		}
# 	}
# 
# 	list(s1 = s1, s2 = s2)
# }
# 
# draw_sp_basis_mc_OLD <- function(R, domain, report.period = 100)
# {
# 	n <- nrow(domain)
# 	s1 <- matrix(NA, n, R)
# 	s2 <- matrix(NA, n, R)
# 
# 	for (j in 1:n) {
# 		if (j %% report.period == 0) {
# 			logger("Drawing points for area %d of %d\n", j, n)
# 		}
# 
# 		P <- rArea(R, domain[j,])
# 		s1[j,] <- P[,1]
# 		s2[j,] <- P[,2]
# 	}
# 
# 	list(s1 = s1, s2 = s2)
# }
# 
# compute_spt_basis_mc_OLD <- function(basis, domain, period, s1, s2, report.period = 100)
# {
# 	n <- nrow(domain)
# 	r <- basis$get_dim()
# 	R <- dim(s1)[3]
# 	stopifnot(R == dim(s2)[3])
# 	S <- Matrix(0, n, r)
# 	T <- length(period)
# 
# 	for (r in 1:R) {
# 		if (r %% report.period == 0) {
# 			logger("Computing basis for rep %d of %d\n", r, R)
# 		}
# 		for (t in 1:T) {
# 			S <- S + basis$compute(s1[t,,r], s2[t,,r], period[t])
# 		}
# 	}
# 
# 	return( S / (R*T) )
# }
# 
# compute_sp_basis_mc_OLD <- function(basis, domain, s1, s2, report.period = 100)
# {
# 	n <- nrow(domain)
# 	r <- basis$get_dim()
# 	R <- dim(s1)[2]
# 	stopifnot(R == dim(s2)[2])
# 	S <- Matrix(0, n, r)
# 
# 	for (r in 1:R) {
# 		if (r %% report.period == 0) {
# 			logger("Computing basis for rep %d of %d\n", r, R)
# 		}
# 		S <- S + basis$compute(s1[,r], s2[,r])
# 	}
# 
# 	return(S / R)
# }

compute_sp_basis_mc <- function(basis, domain, R, report.period = 100)
{
	n <- nrow(domain)
	r <- basis$get_dim()
	S <- Matrix(0, n, r)

	for (j in 1:n) {
		if (j %% report.period == 0) {
			logger("Computing basis for area %d of %d\n", j, n)
		}

		for (t in 1:T) {
			P <- rArea(R, domain[j,])
			S[j,] <- S[j,] + colSums(basis$compute(P[,1], P[,2]))
		}
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

		for (t in 1:T) {
			P <- rArea(R, domain[j,])
			S[j,] <- S[j,] + colSums(basis$compute(P[,1], P[,2], period[t]))
		}
	}

	return( S / (R*T) )
}
