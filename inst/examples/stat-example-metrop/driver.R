library(Matrix)
library(stcos)
library(coda)

source("../../../R/metrop.R")
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

init <- list(sig2mu = 1, sig2K = 1, sig2xi = 1e8)
metrop.out <- metrop.stcos(Z, S, sig2eps, C.inv, H, R = 100,
	report.period = 1, burn = 0, thin = 1, init = init)

mu_B.mcmc <- mcmc(metrop.out$mu_B.hist)
eta.mcmc <- mcmc(metrop.out$eta.hist)
sig2mu.mcmc <- mcmc(metrop.out$sig2mu.hist)
sig2xi.mcmc <- mcmc(metrop.out$sig2xi.hist)
sig2K.mcmc <- mcmc(metrop.out$sig2K.hist)

save.image("results.Rdata")