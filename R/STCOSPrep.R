# The following documentation format was suggested in
# https://stackoverflow.com/questions/45431845/documenting-r6-classes-and-methods-within-r-package-in-rstudio/45603005#

#' STCOS Preparation
#' 
#' An \code{\link{R6Class}} for preparing an STCOS analysis.
#' 
#' @section Usage:
#' \preformatted{
#' sp <- STCOSPrep$new(fine_domain, fine_domain_geo_name, basis, basis_mc_reps = 500,
#'                     report_period = 100)
#' 
#' sp$add_obs(domain, period, estimate_name, variance_name, geo_name)
#' 
#' S <- sp$get_S()
#'
#' sp$set_basis_reduction(f)
#' 
#' sp$get_reduced_S()
#'
#' sp$get_Kinv(times, X = NULL, autoreg = TRUE)
#' 
#' sp$get_obs(idx)
#' 
#' sp$get_basis()
#' 
#' sp$get_geo()
#' }
#' 
#' @section Arguments:
#' \itemize{
#' \item \code{fine_domain} An \code{sf} object; the fine-level support for the analysis.
#' \item \code{fine_domain_geo_name} The name of the field in \code{fine_domain} which represents
#' a unique geographical ID.
#' \item \code{basis} An object of type \code{SpaceTimeBisquareBasis}.
#' \item \code{basis_mc_reps} Number of Monte Carlo reps to use when computing
#' area-level basis function decomposition.
#' \item \code{report_period} Gibbs sampler will report progress each time this many
#' iterations are completed.
#' \item \code{domain} An \code{sf} object; a source support for the analysis.
#' \item \code{period} A vector of times from which the estimates were computed. For
#' example, to indicate 2015 ACS 5-year estimates, use \code{period = 2011:2015}.
#' \item \code{estimate_name} Name of the field which contains direct estimates.
#' \item \code{variance_name} Name of the field which contains direct variance estimates.
#' \item \code{geo_name} Name of the fiend which represents unique geographical ID.
#' \item \code{idx} A vector of indices.
#' \item \code{f} A function which performs a dimension reduction.
#' \item \code{times} A vector of times relevant to the analysis.
#' \item \code{X} A fixed covariate, if one is available.
#' \item \code{autoreg} A boolean; if TRUE, assume an autoregressive covariance
#' structure in time. Otherwise assume independence between times.
#' }
#' 
#' @section Methods:
#' \itemize{
#' \item \code{$new} Create a new \code{STCOSPrep} object, which does not yet
#' contain any source supports.
#' \item \code{$add_obs} Add a source support.
#' \item \code{$get_obs} Get a source support(s) which have already been added.
#' \item \code{$domain2model} Compute \code{H} and \code{S} matrix fragments
#' for a given source support.
#' \item \code{$get_z} Get vector of direct estimates from added source supports.
#' \item \code{$get_v} Get vector of direct variance estimates from added source
#' supports.
#' \item \code{$get_H} Get matrix of overlaps between fine-level support and added
#' source supports.
#' \item \code{$get_S} Get design matrix based on basis function decomposition, based
#' on added source supports.
#' \item \code{$get_reduced_S} Same as \code{get_S}, except first apply the dimension
#' reduction function (which is the identity function by default).
#' \item \code{$get_geo} Get vector of GEO IDs from added source supports.
#' \item \code{$set_basis_reduction} Set the dimension reduction function to be
#' applied to \code{S} and related matrices.
#' \item \code{$get_basis} Get the basis function which has been used to construct
#' this object.
#' \item \code{$get_Kinv} Compute the \code{K.inv} matrix.
#' }
#'
#' @name STCOSPrep
#' 
#' @examples
#' \dontrun{
#' sp <- STCOSPrep$new(dom.fine, "GEO_ID", basis, 500)
#' 
#' # Add source support data
#' sp$add_obs(acs1.2015, 2015, "DirectEst", "DirectVar", "GEO_ID")
#' ...
#' sp$add_obs(acs5.2015, 2011:2015, "DirectEst", "DirectVar", "GEO_ID")
#' 
#' # Reduce dimension of the design matrix with basis expansion
#' S <- sp$get_S()
#' eig <- eigen(t(S) %*% S)
#' idx.S <- which(cumsum(eig$values) / sum(eig$values) < 0.75)
#' Tx.S <- t(eig$vectors[idx.S,])
#' f <- function(S) { S %*% Tx.S }
#' sp$set_basis_reduction(f)
#' 
#' # Get quantities needed for MCMC
#' S.reduced <- sp$get_reduced_S()
#' z <- sp$get_z()
#' v <- sp$get_v()
#' H <- sp$get_H()
#'
#' # compute K.inv matrix using Random-Walk method
#' K.inv <- sp$get_Kinv(2011:2015)
#' 
#' # Retrieve some of the source supports
#' sp$get_obs(idx = c(2,4,7))
#' 
#' # Return the basis function object from which sp was constructed
#' sp$get_basis()
#' 
#' # Return the GEO IDs which have been processed so far
#' sp$get_geo()
#' }
NULL

#' @export
#' @docType class
STCOSPrep <- R6Class("STCOSPrep",
	lock_objects = TRUE,
	lock_class = TRUE,
	private = list(
		fine_domain = NULL,
		fine_domain_geo_name = NULL,
		H_list = NULL,
		S_list = NULL,
		z_list = NULL,
		v_list = NULL,
		geo_list = NULL,
		N = NULL,
		L = NULL,
		basis = NULL,
		basis_reduction = NULL,
		basis_mc_reps = NULL,
		report_period = NULL
	),
	public = list(
		initialize = function(fine_domain, fine_domain_geo_name, basis, basis_mc_reps = 500, report_period = 100) {
			stopifnot(inherits(fine_domain, "sf"))
			stopifnot(inherits(basis, "SpaceTimeBisquareBasis"))
			stopifnot(fine_domain_geo_name %in% colnames(fine_domain))

			private$fine_domain <- fine_domain
			private$H_list <- list()
			private$S_list <- list()
			private$z_list <- list()
			private$v_list <- list()
			private$geo_list <- list()
			private$N <- 0
			private$L <- 0
			private$basis_mc_reps <- basis_mc_reps
			private$basis <- basis
			private$basis_reduction <- identity
			private$report_period <- report_period
			private$fine_domain_geo_name <- fine_domain_geo_name
		},
		add_obs = function(domain, period, estimate_name, variance_name, geo_name)
		{
			stopifnot(class(estimate_name) == "character")
			stopifnot(class(variance_name) == "character")
			stopifnot(geo_name %in% colnames(domain))

			logger("Begin adding observed space-time domain\n")
			out <- self$domain2model(domain, period, geo_name)

			logger("Extracting direct estimates from field '%s'", estimate_name)
			printf(" and variance estimates from field '%s'\n", variance_name)
			z <- domain[[estimate_name]]
			v <- domain[[variance_name]]

			# Update internal state
			private$N <- private$N + nrow(domain)
			private$L <- private$L + 1
			private$z_list[[private$L]] <- z
			private$v_list[[private$L]] <- v
			private$H_list[[private$L]] <- out$H
			private$S_list[[private$L]] <- out$S
			geo <- out$geo
			geo$obs <- private$L
			private$geo_list[[private$L]] <- geo

			logger("Finished adding observed space-time domain\n")
		},
		get_obs = function(idx = 1:private$L)
		{
			z <- private$z_list[idx]
			v <- private$v_list[idx]
			H <- private$H_list[idx]
			S <- private$S_list[idx]
			geo <- private$geo_list[idx]
			S.reduced <- lapply(S, private$basis_reduction)

			list(z = z, v = v, H = H, S = S, geo = geo, S.reduced = S.reduced)
		},
		domain2model = function(domain, period, geo_name)
		{
			## Check the argument types
			stopifnot(inherits(domain, "sf"))
			stopifnot(length(period) >= 1)
			stopifnot(class(geo_name) == "character")

			logger("Computing overlap matrix using field '%s'\n", geo_name)
			H.prime <- compute.overlap(private$fine_domain, domain,
				geo.name.D = private$fine_domain_geo_name, geo.name.G = geo_name)
			H <- t(Matrix(apply(H.prime, 2, normalize)))

			logger("Computing basis functions\n")
			S <- compute_spt_basis_mc(basis = private$basis, domain = domain,
				R = private$basis_mc_reps, period = period,
				report.period = private$report_period)
			S.reduced <- private$basis_reduction(S)

			geo <- data.frame(row = 1:nrow(domain), geo_id = domain[[geo_name]])
			list(H = H, S = S, S.reduced = S.reduced, geo = geo)
		},
		get_z = function()
		{
			n <- nrow(private$fine_domain)
			z <- numeric(private$N)
			L <- private$L
			cnt <- 0

			for (l in 1:L) {
				idx <- 1:length(private$z_list[[l]]) + cnt
				z[idx] <- private$z_list[[l]]
				cnt <- cnt + length(idx)
			}

			return(z)
		},
		get_v = function()
		{
			n <- nrow(private$fine_domain)
			v <- numeric(private$N)
			L <- private$L
			cnt <- 0

			for (l in 1:L) {
				idx <- 1:length(private$v_list[[l]]) + cnt
				v[idx] <- private$v_list[[l]]
				cnt <- cnt + length(idx)
			}

			return(v)
		},
		get_H = function()
		{
			n <- nrow(private$fine_domain)
			H <- Matrix(0, 0, n)
			L <- private$L
			cnt <- 0

			for (l in 1:L) {
				H <- rbind(H, private$H_list[[l]])
				cnt <- cnt + ncol(private$H_list[[l]])
			}

			return(H)
		},
		get_S = function()
		{
			r <- private$basis$get_dim()
			S <- Matrix(0, 0, r)
			L <- private$L
			cnt <- 0

			for (l in 1:L) {
				S <- rbind(S, private$S_list[[l]])
				cnt <- cnt + nrow(private$S_list[[l]])
			}

			return(S)
		},
		get_reduced_S = function()
		{
			private$basis_reduction(self$get_S())
		},
		get_geo = function()
		{
			G <- private$geo_list[[1]]
			L <- private$L

			# Note: this may not be something that should be collapsed. Revisit...
			for (l in setdiff(1:L, 1)) {
				G <- rbind(G, private$geo_list[[l]])
			}

			return(G)
		},
		set_basis_reduction = function(f = identity)
		{
			stopifnot(is.function(f))
			private$basis_reduction <- f
		},
		get_basis = function()
		{
			return(private$basis)
		},
		get_basis_mc_reps = function()
		{
			return(private$basis_mc_reps)
		},
		get_report_period = function()
		{
			return(private$report_period)
		},
		get_basis_reduction = function()
		{
			return(private$basis_reduction)
		},
		get_Kinv = function(times, X = NULL, method = c("moran", "randomwalk", "car", "independence"))
		{
			ll <- max(times) - min(times) + 1
			if (length(times) != ll) {
				stop("times must contain consecutive integers")
			}

			T <- length(times)
			n <- nrow(private$fine_domain)
			r <- private$basis$get_dim()

			if (method == "independence") {
				return(Diagonal(n = ncol(self$get_reduced_S())))
			}

			# Compute basis function for fine-level process
			Sconnector <- Matrix(0, 0, r)
			for (t in 1:T) {
				idx <- 1:n + (t-1)*n
				logger("Constructing S matrix for fine-scale at time %d of %d\n", t, T)
				S <- compute_spt_basis_mc(basis = private$basis, domain = private$fine_domain,
					R = private$basis_mc_reps, period = times[t], report.period = private$report_period)

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

			# Target covariance
			logger("Computing target covariance\n")
			if (method == "car") {
				# Assume covariance structure without dependence over time
				K <- sptcovar.indep(Qinv, Sconnectorf, lag_max = T)
			} else if (method == "randomwalk") {
				# Assume covariance structure with M as identity matrix
				M <- Diagonal(n,1)
				K <- sptcovar.randwalk(Qinv, M, Sconnectorf, lag_max = T)
			} else if (method == "moran") {
				# Assume covariance structure with M computed via Moran's I basis
				if (is.null(X)) {
					knots.sp <- unique(private$basis$get_cutpoints()[,1:2])
					w.sp <- private$basis$get_ws()
					basis.sp <- SpatialBisquareBasis$new(knots.sp[,1], knots.sp[,2], w.sp)
					X <- compute_sp_basis_mc(basis = basis.sp, domain = private$fine_domain,
						R = private$basis_mc_reps, report.period = private$report_period)
				} else {
					stopifnot(nrow(X) == nrow(Sconnectorf))
				}

				P_perp <- Diagonal(nrow(X),1) - X %*% solve(t(X) %*% X, t(X))
				eig <- eigen(P_perp, symmetric = TRUE)
				M <- Re(eig$vectors)
				M <- (M + t(M)) / 2
				K <- sptcovar.vectautoreg(Qinv, M, Sconnectorf, lag_max = T)
			} else {
				stop("Invalid argument for method")
			}

			eig <- eigen(K)
			P <- Re(eig$vectors)
			D <- Re(eig$values)
			D[D < 0] <- 0
			Dinv <- D
			Dinv[D > 0] <- 1 / D[D > 0]
			Kinv <- P %*% (Dinv * t(P))

			return(Kinv)
		}
	)
)
