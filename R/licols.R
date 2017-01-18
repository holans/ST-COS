# [Xsub,idx]=

# Extract a linearly independent set of columns of a given matrix X
# Copied from a Matlab forum and ported to R
# https://www.mathworks.com/matlabcentral/answers/108835-how-to-get-only-linearly-independent-rows-in-a-matrix-or-to-remove-linear-dependency-b-w-rows-in-a-m
# Input:
#  X: The given input matrix
#  tol: A rank estimation tolerance. Default=1e-10
# Output:
#  Xsub: The extracted columns of X
#  idx:  The indices (into X) of the extracted columns

licols <- function(X, tol = 1e-10)
{
	if ~nnz(X) {
		# X has no non-zeros and hence no independent columns
		return(list(Xsub = matrix(NA,0,0), idx = numeric(0)))
	}

	[Q, R, E] = qr(X, 0)
	if (~isvector(R)) {
		diagr = abs(diag(R))
	} else {
		diagr = R(1)
	}

	# Rank estimation
	r = find(diagr >= tol*diagr(1), 1, 'last')
	idx = sort(E(1:r))
	list(Xsub = X[,idx], idx = idx)
}