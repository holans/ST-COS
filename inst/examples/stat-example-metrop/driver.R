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

gibbs.out <- gibbs.stcos(Z, S, sig2eps, C.inv, H, R = 1000,
	report.period = 1, burn = 100, thin = 5)

mu_B.mcmc <- mcmc(gibbs.out$mu_B.hist)
xi.mcmc <- mcmc(gibbs.out$xi.hist)
eta.mcmc <- mcmc(gibbs.out$eta.hist)
sig2mu.mcmc <- mcmc(gibbs.out$sig2mu.hist)
sig2xi.mcmc <- mcmc(gibbs.out$sig2xi.hist)
sig2K.mcmc <- mcmc(gibbs.out$sig2K.hist)
Y.mcmc <- mcmc(gibbs.out$Y.hist)

save.image("results.Rdata")