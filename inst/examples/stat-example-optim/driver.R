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

init <- list(sig2xi = exp(34.17265))
mle.out <- mle.stcos(Z, S, sig2eps, H, init = init,
	optim.control = list(trace = 6))

mu_B.hat <- mle.out$mu.hat
eta.hat <- mle.out$eta.hat
sig2xi.hat <- mle.out$sig2xi.hat
Z.hat <- as.numeric(H %*% mu_B.hat + S %*% eta.hat)

save.image("results.Rdata")
q("no")
