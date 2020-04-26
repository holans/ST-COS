#' Spatial Bisquare Basis
#' 
#' @description
#' Spatial bisquare basis on point data.
#' 
#' @param X An \eqn{n \times 2} numeric matrix with points to evaluate.
#' @param knots An \eqn{R \times 2} numeric matrix with knots for the basis.
#' @param w Radius for the basis.
#'
#' @return An \eqn{n \times R} sparse matrix where the \eqn{(i,j)}th entry
#' represents the \eqn{i}th row of \code{X} evaluated with the \eqn{j}th
#' row of \code{knots}.
#'   
#' @details
#' For each \eqn{\bm{u}} in the set of points to evaluate, compute the basis functions
#' \deqn{
#' \varphi_j(\bm{u},t) =
#' \left[ 1 - \frac{\Vert\bm{u} - \bm{c}_j \Vert^2}{w^2} \right]^2  \cdot
#' I(\Vert \bm{u} - \bm{c}_j \Vert \leq w).
#' }
#' for \eqn{j = 1, \ldots, R}.
#' 
#' If \code{knots} is interpreted as a matrix, the three columns correspond
#' to x-axis and y-axis coordinates. Here, it is assumed that \code{dom}
#' and \code{sf} are based on a common coordinate system.
#' 
#' @examples
#' set.seed(1234)
#' seq_x = seq(0, 1, length.out = 3)
#' seq_y = seq(0, 1, length.out = 3)
#' knots = merge(seq_x, seq_y)
#' x = runif(50)
#' y = runif(50)
#' rad = 0.5
#' 
#' spatial_bisquare(cbind(x,y), knots, rad)
#' 
#' # Plot the knots and the points at which we evaluated the basis
#' plot(knots[,1], knots[,2], pch = 4, cex = 1.5, col = "red")
#' points(x, y, cex = 0.5)
#' 
#' # Draw a circle representing the basis' radius around one of the knot points
#' tseq = seq(0, 2*pi, length=100) 
#' coords = cbind(rad * cos(tseq) + seq_x[2], rad * sin(tseq) + seq_y[2])
#' lines(coords, col = "red")
#' 
#' @export
spatial_bisquare = function(X, knots, w)
{
	stopifnot(ncol(X) == 2)
	stopifnot(ncol(knots) == 2)
	compute_basis_sp(as.matrix(X), as.matrix(knots), w)
}
