STCOSPrep <- R6Class("STCOSPrep",
	public = list(
		initialize = function(fine_domain, basis) {
			stopifnot("sf" %in% class(fine_domain))
			stopifnot("BisquareBasis" %in% class(basis))

			private$fine_domain <- fine_domain
			private$H_list <- list()
			private$S_list <- list()
			private$Z_list <- list()
			private$V_list <- list()
			private$N <- 0
			private$L <- 0
			private$basis <- basis
			private$basis_reduction <- identity
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
		basis = NULL,
		basis_reduction = NULL
	)
)

add_obs <- function(domain, time, period, estimate_name, variance_name)
{
	## Check the argument types
	stopifnot("sf" %in% class(domain))
	stopifnot(class(estimate_name) == "character")
	stopifnot(class(variance_name) == "character")
	stopifnot(length(time) == 1)
	stopifnot(length(period) >= 1)

	logger("Begin adding observed space-time domain\n")
	private$N <- private$N + nrow(domain)
	private$L <- private$L + 1

	L <- private$L
	N <- private$N

	logger("Extracting survey estimates from field '%s'", estimate_name)
	printf(" and variance estimates from field '%s'\n", variance_name)
	private$Z_list[[L]] <- domain[[estimate_name]]
	private$V_list[[L]] <- domain[[variance_name]]
	
	logger("Computing overlap matrix\n")
	H.prime <- compute.overlap.v3(private$fine_domain, domain)
	H <- Matrix(apply(H.prime, 2, normalize))
	private$H_list[[L]] <- H

	logger("Computing basis functions\n")
	message("Write me")
	# for each area in the domain,
	#	Generate either a grid or a random draw of points
	# for b = 1:B
	#	for time = 1:periods
	#		S <- S + basis$compute(domain)
	#	end
	# end
	# S <- S / (B*length(periods))
	r <- basis$get_dim()
	S <- Matrix(0, N, r)
	private$S_list[[L]] <- S

	logger("Finished adding observed space-time domain\n")
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

STCOSPrep$set("public", "get_Z", get_Z)
STCOSPrep$set("public", "get_V", get_V)
STCOSPrep$set("public", "get_H", get_H)
STCOSPrep$set("public", "get_S", get_S)
STCOSPrep$set("public", "get_reduced_S", get_reduced_S)
STCOSPrep$set("public", "add_obs", add_obs)
STCOSPrep$set("public", "set_basis_reduction", set_basis_reduction)
STCOSPrep$lock()
