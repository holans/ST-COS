library(Matrix)

Z <- as.matrix(read.csv("~/Documents/simulations/ST-COS/dat/Zagg.txt", header = FALSE))
sigmavar <- as.matrix(read.csv("~/Documents/simulations/ST-COS/dat/sigmavar.txt", header = FALSE))
Kinv <- as.matrix(read.csv("~/Documents/simulations/ST-COS/dat/Kinv.txt", header = FALSE))

H.el <- read.csv("~/Documents/simulations/ST-COS/dat/H_sparse.txt", header = FALSE)
H <- sparseMatrix(i = H.el[,1], j = H.el[,2], x = H.el[,3])
rm(H.el)

S.el <- read.csv("~/Documents/simulations/ST-COS/dat/S_sparse.txt", header = FALSE)
S <- sparseMatrix(i = S.el[,1], j = S.el[,2], x = S.el[,3])
rm(S.el)
