library(rstan)
library(shinystan)
library(mvtnorm)
library(Matrix)
library(coda)

set.seed(1234)

sig2eps <- as.matrix(read.csv("dat/sigmavar.txt.gz", header = FALSE))[,1]
C.inv <- as.matrix(read.csv("dat/Kinv.txt.gz", header = FALSE))
LamHpinvVH <- read.csv("dat/LamHpinvVH.txt.gz", header = FALSE)[,1]

H.el <- read.csv("dat/H_sparse.txt.gz", header = FALSE)
H <- sparseMatrix(i = H.el[,1], j = H.el[,2], x = H.el[,3])
rm(H.el)

S.el <- read.csv("dat/S1_sparse.txt.gz", header = FALSE)
S <- sparseMatrix(i = S.el[,1], j = S.el[,2], x = S.el[,3])
rm(S.el)

N <- nrow(S)
r <- ncol(S)
n <- ncol(H)

muB.true <- rep(0, n)
eta.true <- rep(0, r)
sig2xi.true <- 100^2
z <- rnorm(N, as.numeric(H %*% muB.true + S %*% eta.true), sqrt(sig2eps + sig2xi.true))

dat <- list(N = N, n = n, r = r, z = z, H = as.matrix(H), S = as.matrix(S), sig2eps = sig2eps)

# NUTS MCMC
fit <- stan(file = "model.stan", data = dat, 
	iter = 1000, chains = 1)
print(fit)
traceplot(fit)
my_sso <- launch_shinystan(fit)

# VB with ADVI
m <- stan_model(file = "model.stan")
fit.vb <- vb(m, data = dat, algorithm = "meanfield",
	elbo_samples = 1000)
print(fit.vb, pars = "mu")
traceplot(fit.vb, pars = "mu")
