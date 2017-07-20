STCOSPrep <- R6Class("STCOSPrep",
	private = list(
		fine_domain = NULL,
		fine_domain_geo_name = NULL,
		H_list = NULL,
		S_list = NULL,
		Z_list = NULL,
		V_list = NULL,
		geo_list = NULL,
		N = NULL,
		L = NULL,
		basis = NULL
	),
	public = list(
		basis_reduction = NULL,
		basis_mc_reps = NULL,
		report_period = NULL,
		initialize = function(fine_domain, fine_domain_geo_name, basis, basis_mc_reps = 500, report_period = 100) {
			stopifnot(inherits(fine_domain, "sf"))
			stopifnot(inherits(basis, "SpaceTimeBisquareBasis"))

			private$fine_domain <- fine_domain
			private$H_list <- list()
			private$S_list <- list()
			private$Z_list <- list()
			private$V_list <- list()
			private$geo_list <- list()
			private$N <- 0
			private$L <- 0
			self$basis_mc_reps <- basis_mc_reps
			private$basis <- basis
			self$basis_reduction <- identity
			self$report_period <- report_period
			private$fine_domain_geo_name <- fine_domain_geo_name
		},
		add_obs = function(domain, period, estimate_name, variance_name, geo_name)
		{
			stopifnot(class(estimate_name) == "character")
			stopifnot(class(variance_name) == "character")

			private$N <- private$N + nrow(domain)
			private$L <- private$L + 1

			logger("Begin adding observed space-time domain\n")
			out <- self$domain2model(domain, period, geo_name)

			logger("Extracting survey estimates from field '%s'", estimate_name)
			printf(" and variance estimates from field '%s'\n", variance_name)
			Z <- domain[[estimate_name]]
			V <- domain[[variance_name]]

			L <- private$L
			N <- private$N
			private$Z_list[[L]] <- Z
			private$V_list[[L]] <- V
			private$H_list[[L]] <- out$H
			private$S_list[[L]] <- out$S
			private$geo_list[[L]] <- out$geo

			logger("Finished adding observed space-time domain\n")
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
				R = self$basis_mc_reps, period = period,
				report.period = self$report_period)
			S.reduced <- self$basis_reduction(S)

			geo <- data.frame(obs = private$L, row = 1:nrow(domain), geo_id = domain[[geo_name]])
			list(H = H, S = S, S.reduced = S.reduced, geo = geo)
		},
		get_Z = function()
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
		},
		get_V = function()
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
			r <- basis$get_dim()
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
			self$basis_reduction(self$get_S())
		},
		get_geo = function()
		{
			G <- private$geo_list[[1]]
			L <- private$L

			warning("THIS MAY NOT BE COLLAPSEABLE!!")
			for (l in seq_int_ordered(2, L)) {
				G <- rbind(G, private$geo_list[[l]])
			}

			return(G)
		},
		set_basis_reduction = function(f = identity)
		{
			stopifnot(is.function(f))
			self$basis_reduction <- f
		},
		get_Cinv = function(times, X = NULL, autoreg = TRUE)
		{
			ll <- max(times) - min(times) + 1
			if (length(times) != ll) {
				stop("times must contain consecutive integers")
			}

			T <- length(times)
			n <- nrow(private$fine_domain)
			r <- private$basis$get_dim()

			Sconnector <- Matrix(0, 0, r)
			for (t in 1:T) {
				idx <- 1:n + (t-1)*n
				logger("Constructing S matrix for fine-scale at time %d of %d\n", t, T)
				S <- compute_spt_basis_mc(basis = private$basis, domain = private$fine_domain,
					R = self$basis_mc_reps, period = times[t], report.period = self$report_period)

				# rbind usually slows performance, but here it's a lot faster
				# than doing Sconnector[idx,] <- S
				Sconnector <- rbind(Sconnector, S)
			}

			# Reduction should be same as the one used on S matrix for the observations
			Sconnectorf <- self$basis_reduction(Sconnector)

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
			if (!autoreg) {
				# Assume no autocovariance between spatial domains
				C <- sptcovar.indep(Qinv, Sconnectorf, lag_max = T)
			} else if (is.null(X)) {
				# Take X to be an identity matrix, which leads to M being an identity matrix
				M <- Diagonal(n,1)

				# Target Covariance
				C <- sptcovar.randwalk(Qinv, M, Sconnectorf, lag_max = T)
			} else {
				# Use the given X to compute M
				licols.out <- licols(as.matrix(X))
				B <- Matrix(licols.out$Xsub)
				P_perp <- Diagonal(nrow(B),1) - B %*% solve(t(B) %*% B, t(B))
				eig <- eigen(P_perp, symmetric = TRUE)
				M <- Re(eig$vectors)
				M <- (M + t(M)) / 2

				# Target Covariance
				logger("Computing target covariance\n")
				C <- sptcovar.vectautoreg(Qinv, M, Sconnectorf, lag_max = T)
			}

			eig <- eigen(C)
			P <- Re(eig$vectors)
			D <- Re(eig$values)
			D[D < 0] <- 0
			Dinv <- D
			Dinv[D > 0] <- 1 / D[D > 0]
			Cinv.higham <- P %*% (Dinv * t(P))

			return(Cinv.higham)
		}
	)
)

# STCOSPrep$set("public", "get_Z", get_Z)
# STCOSPrep$set("public", "get_V", get_V)
# STCOSPrep$set("public", "get_H", get_H)
# STCOSPrep$set("public", "get_S", get_S)
# STCOSPrep$set("public", "get_reduced_S", get_reduced_S)
# STCOSPrep$set("public", "get_geo", get_geo)
# STCOSPrep$set("public", "get_Cinv", get_Cinv)
# STCOSPrep$set("public", "add_obs", add_obs)
# STCOSPrep$set("public", "domain2model", domain2model)
# STCOSPrep$set("public", "set_basis_reduction", set_basis_reduction)
# STCOSPrep$lock()
