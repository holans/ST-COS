metrop.stcos <- function(Z, S, sig2eps, C.inv, H, R,
	report.period = R+1, burn = 0, thin = 1,
	init = NULL, fixed = NULL)
{
	stopifnot(R > burn)

	r <- ncol(S)
	n <- length(Z)
	n_mu <- ncol(H)

	tt.keep <- 0
	R.keep <- ceiling((R - burn) / thin)
	mu_B.hist <- matrix(NA, R.keep, n_mu)
	eta.hist <- matrix(NA, R.keep, r)
	sig2mu.hist <- matrix(NA, R.keep, 1)
	sig2xi.hist <- matrix(NA, R.keep, 1)
	sig2K.hist <- matrix(NA, R.keep, 1)
	Y.hist <- matrix(NA, R.keep, n)

	# Initial values
	if (is.null(init)) { init <- list() }
	if (is.null(init$sig2mu)) { init$sig2mu <- 1 }
	if (is.null(init$sig2xi)) { init$sig2xi <- 1 }
	if (is.null(init$sig2K)) { init$sig2K <- 1 }
	sig2mu <- init$sig2mu
	sig2xi <- init$sig2xi
	sig2K <- init$sig2K

	logger("Begin sampling\n")


	# A workaround, since some of the eigvals in C.inv are negative
	# Perhaps because they were transferred from Matlab via CSV.
	eigvals <- eigen(C.inv)$values
	eigvals[eigvals < 0] <- -1*eigvals[eigvals < 0]
	logdet.Cinv <- sum(log(eigvals))

	Data <- list(y = Z, X = cbind(H, S), sig2eps = sig2eps, logdet.Cinv = logdet.Cinv)

	# Log-posterior for [sig2mu, sig2xi, sig2K | Z]
	logpost <- function(par, Data) {
		sig2mu <- exp(par[1])
		sig2K <- exp(par[2])
		sig2xi <- exp(par[3])

		# G <- Matrix(0, n_mu+r, n_mu+r)
		# G[cbind(1:n_mu, 1:n_mu)] <- sig2mu
		# G[1:r + n_mu, 1:r + n_mu] <- sig2K * solve(C.inv)

		Ginv <- Matrix(0, n_mu+r, n_mu+r)
		Ginv[cbind(1:n_mu, 1:n_mu)] <- 1/sig2mu
		Ginv[1:r + n_mu, 1:r + n_mu] <- 1/sig2K * C.inv

		Omega <- 1/Data$sig2eps + 1/sig2xi
		Gamma <- t(Data$X) %*% (Omega * Data$X) + Ginv
		nu <- solve(Gamma, t(Data$X) %*% (Omega * Data$y))

		# logdet.Omega <- sum(log(eigen(Omega)$values))
		logdet.Omega <- sum(log(Omega))

		# Computing eigen on Gamma is very slow.
		# It is much faster to use chol, then grab the diagonals
		# See http://stackoverflow.com/questions/18117218/alternative-ways-to-calculate-the-determinant-of-a-matrix-in-r
		# logdet.Gamma <- sum(log(eigen(Gamma)$values))
		L <- chol(Gamma)
		logdet.Gamma <- 2 * sum(log(diag(L)))

		logdet.G <- n_mu*log(sig2mu) + r*log(sig2K) - Data$logdet.Cinv

		logC <- 0.5 * logdet.Gamma + 0.5 * logdet.Omega - 0.5 * logdet.G -
			0.5 * t(Data$y) %*% (Omega * Data$y) +
			0.5 * t(nu) %*% Gamma %*% nu

		# TBD: Replace with our prior
		lprior <- dgamma(sig2mu, 1, 1, log = TRUE) +
			dgamma(sig2K, 1, 1, log = TRUE) +
			dgamma(sig2xi, 1, 1, log = TRUE)
		ljacobian <- log(sig2mu) + log(sig2K) + log(sig2xi)
		ret <- as.numeric(logC + lprior + ljacobian)
		if (is.infinite(ret)) browser()
		ret
	}

	rBeta <- function(par, Data) {
		R <- nrow(par)
		phi <- exp(par)
		Beta.draws <- matrix(NA, R, ncol(X))

		for (r in 1:R) {
			G <- diag(1, 2)
			Omega <- diag(1/phi[r,]^2, n)
			Gamma <- t(Data$X) %*% Omega %*% Data$X + solve(G)
			nu <- solve(Gamma, t(Data$X) %*% Omega %*% Data$y)
			Beta.draws[r,] <- rmvnorm(1, nu, solve(Gamma))
		}

		return(Beta.draws)
	} 

	browser()
	par.init <- log(c(init$sig2mu, init$sig2K, init$sig2xi))
	logpost(par.init, Data)

	rw.out <- rwmetrop(par.init = par.init, logpost = logpost, R = R,
		burn = burn, thin = thin, report.period = report.period, Data = Data,
		proposal = list(scale = 0.1, var = 1), use.cpp = FALSE)
	print(rw.out$accept)
	par.mcmc <- mcmc(rw.out$par)
	phi.mcmc <- exp(par.mcmc)
	plot(phi.mcmc)
	summary(phi.mcmc)

	Beta.mcmc <- mcmc(rBeta(par.mcmc, Data))
	plot(Beta.mcmc)

	logger("Finished sampling\n")

	list(mu_B.hist = mu_B.hist,
		 eta.hist = eta.hist,
		 sig2mu.hist = sig2mu.hist,
		 sig2xi.hist = sig2xi.hist,
		 sig2K.hist = sig2K.hist,
		 Y.hist = Y.hist,
		 elapsed.sec = timer
	)
}

