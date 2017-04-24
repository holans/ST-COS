STCOSPrep <- R6Class("STCOSPrep",
	public = list(
		initialize = function(fine_domain) {
			private$fine_domain <- fine_domain
			private$H_list <- list()
			private$S_list <- list()
			private$Z_list <- list()
			private$V_list <- list()
			private$N <- 0
			private$L <- 0
		}
	),
	private = list(
		fine_domain = NULL,
		H_list = NULL,
		S_list = NULL,
		Z_list = NULL,
		V_list = NULL,
		N = NULL,
		L = NULL
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

	logger("Finished adding observed space-time domain\n")
}

del_obs <- function(i) {
	private$H_list[[i]] <- NULL
}

get_H <- function()
{
	n <- nrow(fine_domain)
	H <- matrix(NA, private$N, n)
	L <- private$L
	cnt <- 0

	for (l in 1:L) {
		idx <- 1:nrow(private$H_list[[l]]) + cnt
		H[idx,] <- private$H_list[[l]]
		cnt <- cnt + length(idx)
	}

	return(H)
}

STCOSPrep$set("public", "get_H", get_H)
STCOSPrep$set("public", "add_obs", add_obs)
STCOSPrep$set("public", "del_obs", del_obs)
STCOSPrep$lock()
