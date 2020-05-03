#' CAR Precision Matrix
#' 
#' A convenience function to compute the CAR precision matrix
#' based on a given adjacency matrix.
#' 
#' @param A An adjacency matrix.
#' @param tau The CAR dependency parameter \eqn{\tau \in [0,1]}.
#' See "Value". Default: \code{1}.
#' @param scale Whether to scale matrix entries. See "Value".
#' Default: \code{FALSE}.
#'
#' @return CAR precision matrix.
#' 
#' @details
#' Suppose \eqn{\bm{A}} is an \eqn{n \times n} adjacency matrix and
#' \eqn{
#'   \bm{D} = \textrm{Diag}(\bm{A} \bm{1})
#'   = \textrm{Diag}(a_{1+}, \ldots, a_{n+}),
#' }
#' and let \eqn{\bm{D}^{-}} be a diagonal matrix whose \eqn{(i,i)}th
#' entry is \eqn{1/a_{i+}} if \eqn{a_{i+} > 0} and \eqn{0} otherwise.
#' If \code{scale} is \code{FALSE}, return the CAR precision matrix
#' \deqn{
#'   \bm{Q} = \bm{D} - \tau \bm{A}.
#' }
#' If \code{scale} is \code{TRUE}, return a scaled version
#' \deqn{
#'   \tilde{\bm{Q}} = \bm{D}^{-} \bm{Q}.
#' }
#' Taking \eqn{\tau = 1} corresponds to the Intrinsic CAR
#' precision matrix.
#' 
#' Typically in a modeling context, the precision matrix will be
#' multiplied by a scaling parameter; e.g., a CAR model for
#' random effects \eqn{\bm{\phi}} could be
#' \deqn{
#'   f(\bm{\phi} \mid \alpha) \propto
#'   \alpha^{-q} \exp\left\{ -\frac{1}{2 \alpha^2}
#'   \bm{\phi}^\top \bm{Q} \bm{\phi} \right\}.
#' }
#' where \eqn{q = \textrm{rank}(Q)}.
#' 
#' @examples
#' data("acs_sf")
#' dom = acs5_2013[1:4,]
#' A = adjacency_matrix(dom)
#' Q = car_precision(A)
#' 
#' @export
car_precision = function(A, tau = 1, scale = FALSE)
{
	stopifnot(0 <= tau && tau <= 1)
	D = Diagonal(x = rowSums(A))
	Q = D - tau*A
	if (scale) {
		dd = rowSums(A) + (rowSums(A) == 0)
		return(1/dd * Q)
	} else {
		return(Q)
	}
}
