# Note: OverdispersionModelsInR is a package written by Andrew.
# Recent versions with rwmetrop were never made public.
library(OverdispersionModelsInR)
library(coda)

# ----- A simple test case -----
n <- 400
X <- cbind(1, x1 = rnorm(n))
Beta.true <- c(1,2)
sigma.true <- 0.75
y <- rnorm(n, X %*% Beta.true, sigma.true)

plot(X[,2], y); abline(Beta.true, col = "red", lty = 2, lwd = 2)

Data <- list(y = y, X = X)
logpost <- function(par, Data) {
	phi <- exp(par)
	G <- diag(1, 2)
	Omega <- diag(1/phi^2, n)
	Gamma <- t(Data$X) %*% Omega %*% Data$X + solve(G)
	nu <- solve(Gamma, t(Data$X) %*% Omega %*% Data$y)

	# logdet.Omega <- sum(log(eigen(Omega)$values))
	logdet.Omega <- n*log(1 / phi^2)
	logdet.Gamma <- sum(log(eigen(Gamma)$values))
	logdet.G <- sum(log(eigen(G)$values))

	logC <- 0.5 * logdet.Gamma + 0.5 * logdet.Omega - 0.5 * logdet.G -
		0.5 * t(Data$y) %*% Omega %*% Data$y +
		0.5 * t(nu) %*% Gamma %*% nu
	lprior <- dgamma(phi, 1, 1, log = TRUE)
	ljacobian <- log(phi)
	ret <- logC + lprior + ljacobian
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

rw.out <- rwmetrop(par.init = log(0.25), logpost = logpost, R = 3000,
	burn = 500, thin = 5, report.period = 100, Data = Data,
	proposal = list(scale = 0.1, var = 1))
print(rw.out$accept)
par.mcmc <- mcmc(rw.out$par)
phi.mcmc <- exp(par.mcmc)
plot(phi.mcmc)
summary(phi.mcmc)

Beta.mcmc <- mcmc(rBeta(par.mcmc, Data))
plot(Beta.mcmc)
