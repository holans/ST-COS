library(Matrix)

Z <- as.matrix(read.csv("Zagg.txt.gz", header = FALSE))
sig2eps <- as.matrix(read.csv("sigmavar.txt.gz", header = FALSE))[,1]
C.inv <- as.matrix(read.csv("Kinv.txt.gz", header = FALSE))

H.el <- read.csv("H_sparse.txt.gz", header = FALSE)
H <- sparseMatrix(i = H.el[,1], j = H.el[,2], x = H.el[,3])
rm(H.el)

S.el <- read.csv("S1_sparse.txt.gz", header = FALSE)
S <- sparseMatrix(i = S.el[,1], j = S.el[,2], x = S.el[,3])
rm(S.el)
