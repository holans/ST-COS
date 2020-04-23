#' @export
car_precision = function(A, tau = 1, scale = FALSE)
{
	D = Diagonal(x = rowSums(A))
	Q = D - tau*A
	if (scale) {
		dd = rowSums(A) + (rowSums(A) == 0)
		return(1/dd * Q)
	} else {
		return(Q)
	}
}
