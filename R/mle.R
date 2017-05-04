mle.stcos <- function(Z, S, V, H, init = NULL,
	optim.control = list())
{
	n <- ncol(H)
	r <- ncol(S)

	# Initial values
	if (is.null(init)) { init <- list() }
	if (is.null(init$sig2xi)) { init$sig2xi <- 1 }

	loglik <- function(theta, Data) {
		sig2xi <- exp(theta)
		Omega <- 1/(Data$V + sig2xi)
		XtOmega <- t(Omega * Data$X)
		Beta <- solve(XtOmega %*% Data$X, XtOmega %*% Data$y)
		mean <- as.numeric(Data$X %*% Beta)
		sum(dnorm(Data$y, mean, sqrt(1/Omega), log = TRUE))
	}

	Data <- list(y=Z, X=cbind(H,S), V=V)
	optim.control$fnscale <- -1
	st <- Sys.time()
	res <- optim(par = log(init$sig2xi), loglik, method = "L-BFGS-B",
		control = optim.control, Data = Data)
	elapsed.sec <- as.numeric(Sys.time() - st, unit = "secs")

	sig2xi.hat <- exp(res$par)
	Omega.hat <- 1/(Data$V + sig2xi.hat)
	XtOmega <- t(Omega.hat * Data$X)
	Beta.hat <- solve(XtOmega %*% Data$X, XtOmega %*% Data$y)
	mu.hat <- Beta.hat[1:n]
	eta.hat <- Beta.hat[1:r + n]

	list(sig2xi.hat = sig2xi.hat, mu.hat = mu.hat, eta.hat = eta.hat,
		res = res, elapsed.sec = elapsed.sec)
}

