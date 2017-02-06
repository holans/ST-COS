# This code was copied from my OverdispersionModelsInR package,
# which was adapted from the LearnBayes package

rwmetrop <- function(par.init, logpost, Data, proposal, R,
	grp = rep(1, length(par.init)), burn = 0, thin = 1,
	report.period = 1000)
{
	stopifnot(R > burn)
	stopifnot(length(par.init) == nrow(proposal$var))
	stopifnot(length(par.init) == ncol(proposal$var))

	qq <- length(par.init)
	R.keep <- ceiling((R - burn) / thin)
	idx.sample <- 0
	par.hist <- matrix(NA, R.keep, qq)

	b <- par.init
	logfb <- logpost(b, Data)
	V.proposal.half <- chol(proposal$var)

	grp.list <- list()
	G <- length(unique(grp))
	for (g in 1:G) {
		grp.list[[g]] <- which(grp == unique(grp)[g])
	}
	accept.grp <- rep(0, G)

	for (r in 1:R) {
		bc <- b + (proposal$scale * t(V.proposal.half)) %*% rnorm(qq)
		for (g in 1:G) {
			idx.grp <- grp.list[[g]]
			b_ <- b
			b_[idx.grp] <- bc[idx.grp]

			log.alpha <- logpost(b_, Data) - logfb
			if (!is.na(log.alpha)) {
				if (log(runif(1)) < log.alpha) {
					b[idx.grp] <- b_[idx.grp]
					accept.grp[g] <- accept.grp[g] + 1
					logfb <- logpost(b, Data)
				}
			} else {
				warning("log.alpha = NA")
			}
		}

		if (r > burn && r %% thin == 0) {
			idx.sample <- idx.sample + 1
			par.hist[idx.sample,] <- b
		}

		if (r %% report.period == 0) {
			acc <- paste(sprintf("%0.02f", accept.grp / r * 100), collapse = ", ")
			logger("After %d rep, accept%% {%s}\n", r, acc)
		}
	}

	list(par = par.hist, accept = accept.grp / R)
}

laplace <- function (logpost, mode, Data, optim.control = list(), optim.method = "L-BFGS-B")
{
	optim.control$fnscale <- -1
	fit <- optim(mode, logpost, gr = NULL, Data, hessian = TRUE,
		method = optim.method, control = optim.control)

	mode <- fit$par
	H <- -solve(fit$hessian)
	p <- length(mode)
	int <- p/2 * log(2 * pi) + 0.5 * log(det(H)) + logpost(mode, Data)

	list(mode = mode, var = H, int = int, converge = (fit$convergence == 0), optim.out = fit)
}

