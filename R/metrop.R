metrop.stcos <- function(Z, S, sig2eps, C.inv, H, R,
	proposal, report.period = R+1, burn = 0, thin = 1,
	init = NULL, fixed = NULL, metrop.grp = c(1,1,1))
{
	stopifnot(R > burn)

	r <- ncol(S)
	n <- length(Z)
	n_mu <- ncol(H)

	# Initial values
	if (is.null(init)) { init <- list() }
	if (is.null(init$sig2mu)) { init$sig2mu <- 1 }
	if (is.null(init$sig2xi)) { init$sig2xi <- 1 }
	if (is.null(init$sig2K)) { init$sig2K <- 1 }

	timer <- list(phi = 0, beta = 0)

	# Log-posterior for [sig2mu, sig2xi, sig2K | Z]
	logpost <- function(par, Data) {
		sig2mu <- exp(par[1])
		sig2K <- exp(par[2])
		sig2xi <- exp(par[3])

		# sig2mu <- init$sig2mu
		# sig2K <- init$sig2K
		# sig2xi <- init$sig2xi

		# G <- Matrix(0, n_mu+r, n_mu+r)
		# G[cbind(1:n_mu, 1:n_mu)] <- sig2mu
		# G[1:r + n_mu, 1:r + n_mu] <- sig2K * solve(C.inv)

		# st <- Sys.time()
		Ginv <- matrix(0, Data$n_mu + Data$r, Data$n_mu + Data$r)
		Ginv[cbind(1:Data$n_mu, 1:Data$n_mu)] <- 1/sig2mu
		Ginv[1:Data$r + Data$n_mu, 1:Data$r + Data$n_mu] <- 1/sig2K * Data$C.inv
		# logger("Ginv took %f secs\n", as.numeric(Sys.time() - st, unit = "secs"))

		# st <- Sys.time()
		Omega <- 1/(Data$sig2eps + sig2xi)
		XtOmega <- t(Omega * Data$X)
		Gamma <- as.matrix(XtOmega %*% Data$X + Ginv)
		# logger("Gamma took %f secs\n", as.numeric(Sys.time() - st, unit = "secs"))

		# st <- Sys.time()
		nu <- solve(Gamma, XtOmega %*% Data$y)
		# logger("nu took %f secs\n", as.numeric(Sys.time() - st, unit = "secs"))

		# st <- Sys.time()
		logdet.Omega <- sum(log(Omega))
		# logger("logdet.Omega took %f secs\n", as.numeric(Sys.time() - st, unit = "secs"))

		# st <- Sys.time()
		# Use cholesky decomposition and grab the diagonals to get the determinant.
		# Seems to be much faster than directly using eigen.
		# See http://stackoverflow.com/questions/18117218/alternative-ways-to-calculate-the-determinant-of-a-matrix-in-r
		
		logdet.Gamma <- tryCatch({
			L <- chol(Gamma)
			2 * sum(log(diag(L)))
		}, error = function(e) {
			logger("chol(Gamma) failed, trying eigen...\n")
			eig <- eigen(Gamma)$values
			sum(log(eig[eig > 0]))
		})
		# logger("logdet.Gamma took %f secs\n", as.numeric(Sys.time() - st, unit = "secs"))

		# st <- Sys.time()
		logdet.G <- Data$n_mu*log(sig2mu) + Data$r*log(sig2K) - Data$logdet.Cinv
		# logger("logdet.G took %f secs\n", as.numeric(Sys.time() - st, unit = "secs"))

		# st <- Sys.time()
		logC <- 0.5 * logdet.Omega - 0.5 * logdet.G - 0.5 * logdet.Gamma -
			0.5 * t(Data$y) %*% (Omega * Data$y) +
			0.5 * t(nu) %*% Gamma %*% nu
		# logger("logC took %f secs\n", as.numeric(Sys.time() - st, unit = "secs"))

		lprior <- dinvgamma(sig2mu, 2, 1, log = TRUE) +
			dinvgamma(sig2K, 2, 1e-2, log = TRUE) +
			dinvgamma(sig2xi, 2, 1e8, log = TRUE)
		ltx <- log(sig2mu) + log(sig2K) + log(sig2xi)
		ret <- as.numeric(logC + lprior + ltx)
		printf("sig2mu=%e, sig2K=%e, sig2xi=%e, logpost=%e\n", sig2mu, sig2K, sig2xi, ret)
		if (is.infinite(ret)) browser()
		return(ret)
	}

	rBeta <- function(par, Data) {
		R <- nrow(par)
		Beta.draws <- matrix(NA, R, ncol(Data$X))

		for (s in 1:R) {
			if (s %% report.period == 0) { logger("Drawing beta for rep %d of %d\n", s, R) }
			sig2mu <- exp(par[s,1])
			sig2K <- exp(par[s,2])
			sig2xi <- exp(par[s,3])

			Ginv <- matrix(0, Data$n_mu+r, Data$n_mu+r)
			Ginv[cbind(1:Data$n_mu, 1:Data$n_mu)] <- 1/sig2mu
			Ginv[1:Data$r + Data$n_mu, 1:Data$r + Data$n_mu] <- 1/sig2K * Data$C.inv

			Omega <- 1/(Data$sig2eps + sig2xi)
			XtOmega <- t(Omega * Data$X)
			Gamma <- as.matrix(XtOmega %*% Data$X + Ginv)
			Gamma.inv <- solve(Gamma)
			# nu <- solve(Gamma, XtOmega %*% Data$y)
			nu <- Gamma.inv %*% (XtOmega %*% Data$y)

			zz <- rnorm(ncol(Data$X), 0, 1)
			Beta.draws[s,] <- as.numeric(chol(Gamma.inv) %*% zz + nu)
		}

		return(Beta.draws)
	} 

	# A workaround, since some of the eigvals in C.inv are negative
	# Perhaps because they were transferred from Matlab via CSV.
	eig <- eigen(C.inv)
	eigvals <- eig$values
	idx <- which(eigvals > 0)
	logdet.Cinv <- sum(log(eigvals[idx]))
	Data <- list(y = Z, X = cbind(H, S), sig2eps = sig2eps, C.inv = eig$vectors[,idx] %*% (eigvals[idx] * t(eig$vectors[,idx])),
		logdet.Cinv = logdet.Cinv, n = n, r = r, n_mu = n_mu)
	par.init <- log(c(init$sig2mu, init$sig2K, init$sig2xi))

	logger("Begin sampling sig2mu, sig2K, and sig2xi\n")
	st <- Sys.time()

	# We could try a laplace approximation first ...
	# laplace.out <- laplace(logpost, par.init, Data)
	# par.init2 <- laplace.out$mode

	rw.out <- rwmetrop(par.init = par.init, logpost = logpost, R = R,
		burn = burn, thin = thin, report.period = report.period, Data = Data,
		proposal = proposal, grp = metrop.grp)
	timer$phi <- as.numeric(Sys.time() - st, unit = "secs")
	par.hist <- rw.out$par
	phi.hist <- exp(par.hist)
	logger("Finished sampling sig2mu, sig2K, and sig2xi\n")

	if (FALSE) {
		# TBD: This part is slow because of the Gamma inversion
		# It would be easy to do in parallel across the draws
		logger("Begin sampling mu and eta\n")
		st <- Sys.time()
		Beta.hist <- rBeta(par.hist, Data)
		timer$beta <- as.numeric(Sys.time() - st, unit = "secs")
		logger("Finished sampling mu and eta\n")
	} else {
		# Just return NA for now
		Beta.hist <- matrix(NA, nrow(rw.out$par), ncol(Data$X))
	}

	# TBD: Remove this part after debugging
	if (FALSE) {
		print(rw.out$accept)
		par.mcmc <- mcmc(rw.out$par)
		phi.mcmc <- exp(par.mcmc)
		plot(par.mcmc)
		plot(phi.mcmc)
		summary(phi.mcmc)
		Beta.mcmc <- mcmc(Beta.hist)
		plot(Beta.mcmc[,1:3])
	}

	list(mu_B.hist = Beta.hist[,1:n_mu],
		 eta.hist = Beta.hist[,1:r + n_mu],
		 sig2mu.hist = phi.hist[,1],
		 sig2K.hist = phi.hist[,2],
		 sig2xi.hist = phi.hist[,3],
		 elapsed.sec = timer
	)
}
