# [Xsub,idx]=

# Extract a linearly independent set of columns of a given matrix X
# Copied from a Matlab forum and ported to R
# https://www.mathworks.com/matlabcentral/answers/108835-how-to-get-only-linearly-independent-rows-in-a-matrix-or-to-remove-linear-dependency-b-w-rows-in-a-m
# Input:
#   X: The given input matrix
#   tol: A rank estimation tolerance. Default=1e-10
# Output:
#   Xsub: The extracted columns of X
#   idx:  The indices (into X) of the extracted columns
#
# Important Note: I get the same result as Jon if I take X to be a full
# (not sparse) matrix, but it takes a much much longer time to compute
# the QR decomposition than if X is sparse. The sparse QR in R seems
# to give different results than the dense or sparse QR in Matlab.
# Running non-sparse QR on the 32836 x 4750 S matrix takes about 15
# minutes in R; maybe not much less in Matlab.
licols <- function(X, tol = 1e-10)
{
	stopifnot(class(X) == "matrix")
	
	if (norm(X, type = "F") == 0) {
		# X has no non-zeros and hence no independent columns
		return(list(Xsub = matrix(NA,0,0), idx = numeric(0)))
	}

	system.time(qr.out <- qr(X, LAPACK = TRUE, tol = 1e-10))
	# diagr <- abs(diag(qr.R(qr.out)))
	diagr <- abs(diag(qr.R(qr.out)))
	mm <- max(diagr)
	# diagr <- diagr[diagr > 0]
	r <- sum(diagr >= tol*mm)

	r <- 1250
	idx <- qr.out$pivot[1:r]
	list(Xsub = X[,idx], idx = idx)
}