#' licols
#' 
#' Extract a linearly independent set of columns of a matrix.
#'
#' @param X A matrix.
#' @param tol A tolerance for rank estimation. Default is 1e-10.
#' @param quiet logical; if FALSE, print a warning about computation time if \code{X} is large.
#'
#' @return \code{Xsub} contains the extracted columns of \code{X} and \code{idx}
#' contains the indices (into X) of those columns. The elapsed time is stored in
#' \code{elapsed.sec}.
#'
#' @details
#' An R version of a Matlab \code{licols} function given in
#' \href{https://www.mathworks.com/matlabcentral/answers/108835-how-to-get-only-linearly-independent-rows-in-a-matrix-or-to-remove-linear-dependency-b-w-rows-in-a-m#answer_117458}{this MathWorks forum post}.
#'
#' @examples
#' x = 0:19 %% 3 + 1
#' Z = model.matrix(~ as.factor(x) - 1)
#' X = cbind(1, Z)
#' licols(X)
#' 
#' @export
licols = function(X, tol = 1e-10, quiet = FALSE)
{
	# Note: It takes a much much longer time to compute
	# the QR decomposition than if X is sparse. The sparse QR in R seems
	# to give different results than the dense or sparse QR in Matlab.
	# Running non-sparse QR on the 32836 x 4750 S matrix takes about 15
	# minutes in R; maybe not much less in Matlab.
	#
	# Another note: The sparse QR algorithm in the matrix package doesn't
	# return the same kind of pivot. It appears to be a "non-rank-revealing
	# decomposition".

	stopifnot("matrix" %in% class(X))
	m = nrow(X)
	n = ncol(X)
	if (m*n > 1e6 && !quiet) {
		msg = sprintf("QR for %d x %d matrix might take a while...", m, n)
		message(msg)
	}

	if (norm(X, type = "F") == 0) {
		# X has no non-zeros and hence no independent columns
		return(list(Xsub = matrix(NA,0,0), idx = numeric(0)))
	}

	st = Sys.time()
	qr.out = qr(X, LAPACK = TRUE, tol = tol)
	diagr = abs(diag(qr.R(qr.out)))
	mm = max(diagr)
	r = sum(diagr >= tol*mm)
	idx = qr.out$pivot[1:r]

	elapsed.sec = as.numeric(Sys.time() - st, units = "secs")
	list(Xsub = X[,idx,drop=FALSE], idx = idx, elapsed.sec = elapsed.sec)
}
