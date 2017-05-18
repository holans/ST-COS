STCOSPrep <- R6Class("STCOSPrep",
	public = list(
		initialize = function(fine_domain, basis, basis_mc_reps = 500, report_period = 100) {
			stopifnot(inherits(fine_domain, "sf"))
			stopifnot(inherits(basis, "SpaceTimeBisquareBasis"))

			private$fine_domain <- fine_domain
			private$H_list <- list()
			private$S_list <- list()
			private$Z_list <- list()
			private$V_list <- list()
			private$N <- 0
			private$L <- 0
			private$basis_mc_reps <- basis_mc_reps
			private$basis <- basis
			private$basis_reduction <- identity
			private$report_period <- report_period
		}
	),
	private = list(
		fine_domain = NULL,
		H_list = NULL,
		S_list = NULL,
		Z_list = NULL,
		V_list = NULL,
		N = NULL,
		L = NULL,
		basis_mc_reps = NULL,
		basis = NULL,
		basis_reduction = NULL,
		report_period = NULL
	)
)

add_obs <- function(domain, time, period, estimate_name, variance_name)
{
	## Check the argument types
	stopifnot(inherits(domain, "sf"))
	stopifnot(class(estimate_name) == "character")
	stopifnot(class(variance_name) == "character")
	stopifnot(length(time) == 1)
	stopifnot(length(period) >= 1)

	logger("Begin adding observed space-time domain\n")
	n <- nrow(domain)

	logger("Extracting survey estimates from field '%s'", estimate_name)
	printf(" and variance estimates from field '%s'\n", variance_name)
	Z <- domain[[estimate_name]]
	V <- domain[[variance_name]]

	logger("Computing overlap matrix\n")
	H.prime <- compute.overlap(private$fine_domain, domain)
	H <- Matrix(apply(H.prime, 2, normalize))

	logger("Computing basis functions\n")
	draws.out <- private$draw_basis_mc(domain, length(period))
	S <- private$compute_spt_basis_mc(domain, period, draws.out$s1, draws.out$s2)

	private$N <- private$N + n
	private$L <- private$L + 1

	L <- private$L
	N <- private$N
	private$Z_list[[L]] <- Z
	private$V_list[[L]] <- V
	private$H_list[[L]] <- H
	private$S_list[[L]] <- S

	logger("Finished adding observed space-time domain\n")
}

draw_basis_mc <- function(domain, period_len)
{
	basis <- private$basis
	n <- nrow(domain)
	r <- basis$get_dim()
	R <- private$basis_mc_reps
	T <- period_len
	s1 <- array(NA, dim = c(T, n, R))
	s2 <- array(NA, dim = c(T, n, R))	

	for (j in 1:n) {
		if (j %% private$report_period == 0) {
			logger("Drawing points for area %d of %d\n", j, n)
		}
		for (t in 1:T) {
			P <- rArea(R, domain[j,])
			s1[t,j,] <- P[,1]
			s2[t,j,] <- P[,2]
		}
	}

	list(s1 = s1, s2 = s2)
}

compute_spt_basis_mc <- function(domain, period, s1, s2)
{
	basis <- private$basis
	n <- nrow(domain)
	r <- basis$get_dim()
	R <- dim(s1)[3]
	stopifnot(R == dim(s2)[3])
	S <- Matrix(0, n, r)
	T <- length(period)

	for (r in 1:R) {
		if (r %% private$report_period == 0) {
			logger("Computing basis for rep %d of %d\n", r, R)
		}
		for (t in 1:T) {
			S <- S + basis$compute(s1[t,,r], s2[t,,r], period[t])
		}
	}

	return( S / (R*T) )
}

compute_sp_basis_mc <- function(basis, domain, s1, s2)
{
	n <- nrow(domain)
	r <- basis$get_dim()
	R <- dim(s1)[2]
	stopifnot(R == dim(s2)[2])
	S <- Matrix(0, n, r)

	for (r in 1:R) {
		if (r %% private$report_period == 0) {
			logger("Computing basis for rep %d of %d\n", r, R)
		}
		S <- S + basis$compute(s1[,r], s2[,r])
	}

	return(S / R)
}

get_Z <- function()
{
	n <- nrow(private$fine_domain)
	Z <- numeric(private$N)
	L <- private$L
	cnt <- 0

	for (l in 1:L) {
		idx <- 1:length(private$Z_list[[l]]) + cnt
		Z[idx] <- private$Z_list[[l]]
		cnt <- cnt + length(idx)
	}

	return(Z)
}

get_V <- function()
{
	n <- nrow(private$fine_domain)
	V <- numeric(private$N)
	L <- private$L
	cnt <- 0

	for (l in 1:L) {
		idx <- 1:length(private$V_list[[l]]) + cnt
		V[idx] <- private$V_list[[l]]
		cnt <- cnt + length(idx)
	}

	return(V)
}

get_H <- function()
{
	n <- nrow(private$fine_domain)
	H <- Matrix(0, private$N, n)
	L <- private$L
	cnt <- 0

	for (l in 1:L) {
		idx <- 1:ncol(private$H_list[[l]]) + cnt
		H[idx,] <- t(private$H_list[[l]])
		cnt <- cnt + length(idx)
	}

	return(H)
}

get_S <- function()
{
	r <- basis$get_dim()
	S <- Matrix(0, private$N, r)
	L <- private$L
	cnt <- 0

	for (l in 1:L) {
		idx <- 1:nrow(private$S_list[[l]]) + cnt
		S[idx,] <- private$S_list[[l]]
		cnt <- cnt + length(idx)
	}

	return(S)
}

get_reduced_S <- function()
{
	private$basis_reduction(self$get_S())
}

set_basis_reduction <- function(f = identity)
{
	stopifnot(is.function(f))
	private$basis_reduction <- f
}

# TBD: Should we take target *times* or target *periods*???
# Should Cinv ever be based on 5-year ACS means, for example?
get_Cinv <- function(target.periods)
{
	T <- length(target.periods)
	n <- nrow(private$fine_domain)
	r <- private$basis$get_dim()

	Sconnector <- Matrix(0, 0, r)
	draws.out <- private$draw_basis_mc(private$fine_domain, 1)
	for (t in 1:T) {
		idx <- 1:n + (t-1)*n
		logger("Constructing S matrix for fine-scale at time %d of %d\n", t, T)
		S <- private$compute_spt_basis_mc(private$fine_domain,
			target.periods[t], draws.out$s1, draws.out$s2)
		# Sconnector[idx,] <- S

		# rbind usually slows performance, but here it's a lot faster
		# than doing Sconnector[idx,] <- S
		Sconnector <- rbind(Sconnector, S)
	}

	# Reduction should be same as the one used on S matrix for the observations
	Sconnectorf <- private$basis_reduction(Sconnector)

	# Compute adjacency matrix
	logger("Computing adjacency matrix\n")
	out <- st_touches(private$fine_domain, private$fine_domain)
	A <- adjList2Matrix(out)
	countAdj <- Matrix(0, nrow(A), ncol(A))
	s <- rowSums(A)
	for (j in 1:n) {
		if (s[j] > 0) {
			countAdj[j,] = A[j,] / s[j]
		}
	}
	Q <- Diagonal(n,1) - 0.9*countAdj
	Qinv <- solve(Q)

	# Moran's I Propagator
	# With this choice of B, M just is the identity matrix
	# B <- cbind(diag(N), diag(N))
	# P_perp = diag(nrow(B)) - B %*% MASS::ginv(t(B) %*% B) %*% t(B)
	# eig = eigen(P_perp)
	# M = Re(eig$vectors)
	#### M <- Diagonal(n,1)
	
	# Moran's I Propagator (Experimental)
	# Do a spatial-only basis expansion of fine-domain, and use this as the
	# design matrix to project away from
	warning("We're trying to use variable knots from outside environment. Pass this correctly if we keep it!")
	logger("Computing Moran's I Propagator\n")
	sp.basis <- SpatialBisquareBasis$new(knots[1:250,1], knots[1:250,2], w = 1)
	X <- private$compute_sp_basis_mc(sp.basis, private$fine_domain, draws.out$s1[1,,], draws.out$s2[1,,])
	licols.out <- licols(as.matrix(X))
	B <- Matrix(licols.out$Xsub)
	P_perp <- Diagonal(nrow(B),1) - B %*% solve(t(B) %*% B, t(B))
	eig <- eigen(P_perp)
	M <- Re(eig$vectors)

	# Target Covariance
	# TBD: Should we multiply by a constant to make the elements' magnitude less extreme?
	logger("Computing target covariance\n")
	C.unscaled <- sptcovar(Qinv, M, Sconnectorf, lag_max = T)
	C <- C.unscaled / max(abs(as.matrix(C.unscaled)))
	warning("We're scaling C by a constant. Make sure this is okay!")
	eig <- eigen(C)
	P <- Re(eig$vectors)
	D <- Re(eig$values)
	D[D < 0] <- 0
	Dinv <- D
	Dinv[D > 0] <- 1 / D[D > 0]
	Cinv.higham <- P %*% (Dinv * t(P))

	return(Cinv.higham)
}

STCOSPrep$set("private", "draw_basis_mc", draw_basis_mc)
STCOSPrep$set("private", "compute_sp_basis_mc", compute_sp_basis_mc)
STCOSPrep$set("private", "compute_spt_basis_mc", compute_spt_basis_mc)
STCOSPrep$set("public", "get_Z", get_Z)
STCOSPrep$set("public", "get_V", get_V)
STCOSPrep$set("public", "get_H", get_H)
STCOSPrep$set("public", "get_S", get_S)
STCOSPrep$set("public", "get_reduced_S", get_reduced_S)
STCOSPrep$set("public", "get_Cinv", get_Cinv)
STCOSPrep$set("public", "add_obs", add_obs)
STCOSPrep$set("public", "set_basis_reduction", set_basis_reduction)
STCOSPrep$lock()
