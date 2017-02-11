library(Matrix)
library(stcos)
library(coda)

setwd("~/Documents/simulations/ST-COS/")

Z <- as.matrix(read.csv("dat/Zagg.txt.gz", header = FALSE))
sig2eps <- as.matrix(read.csv("dat/sigmavar.txt.gz", header = FALSE))[,1]
C.inv <- as.matrix(read.csv("dat/Kinv.txt.gz", header = FALSE))
LamHpinvVH <- read.csv("dat/LamHpinvVH.txt.gz", header = FALSE)[,1]

H.el <- read.csv("dat/H_sparse.txt.gz", header = FALSE)
H <- sparseMatrix(i = H.el[,1], j = H.el[,2], x = H.el[,3])
rm(H.el)

S.el <- read.csv("dat/S1_sparse.txt.gz", header = FALSE)
S <- sparseMatrix(i = S.el[,1], j = S.el[,2], x = S.el[,3])
rm(S.el)

# First compute the MLE
init <- list(sig2xi = exp(19))
mle.out <- mle.stcos(Z, S, sig2eps, H, init = init,
        optim.control = list(trace = 6))
sig2xi.hat <- mle.out$sig2xi.hat

# Now do Bayesiam with sig2xi fixed to MLE
init <- list(sig2mu = 1, sig2K = 1e-2, sig2xi = sig2xi.hat)
fixed <- list(sig2xi = TRUE)
hyper <- list(a.sig2mu = 2, a.sig2K = 2, a.sig2xi = 2,
	b.sig2mu = 1e12, b.sig2K = 2, b.sig2xi = 1e8)
gibbs.stcos(Z, S, sig2eps, C.inv, H, R = 100,
	report.period = 1, burn = 0, thin = 1,
	init = init, fixed = fixed, hyper = hyper)

mu_B.mcmc <- mcmc(metrop.out$mu_B.hist)
eta.mcmc <- mcmc(metrop.out$eta.hist)
sig2mu.mcmc <- mcmc(metrop.out$sig2mu.hist)
sig2xi.mcmc <- mcmc(metrop.out$sig2xi.hist)
sig2K.mcmc <- mcmc(metrop.out$sig2K.hist)

save.image("results.Rdata")

