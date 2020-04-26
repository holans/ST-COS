#' Space-Time Bisquare Basis
#' 
#' @description
#' Space-time bisquare basis on point data.
#' 
#' @param dom An \eqn{n \times 3} numeric matrix with points \eqn{(\bm{u},v)}
#' to evaluate.
#' @param knots_s An \eqn{R \times 2} numeric matrix with spatial knots
#' \eqn{\bm{c}_1, \ldots, \bm{c}_R} for the basis.
#' @param knots_t A numeric vector with temporal knots \eqn{g_1, \ldots, g_T}
#' for the basis.
#' @param w_s Spatial radius for the basis.
#' @param w_t Temporal radius for the basis.
#'
#' @return A sparse \eqn{n \times RT} matrix whose \eqn{i}th row
#' represents the \eqn{i}th point \eqn{(\bm{u}_i,v_i)} evaluated at every
#' basis function for \eqn{j = 1, \ldots, R} and \eqn{t = 1, \ldots, T}.
#' 
#' @details
#' For each \eqn{(\bm{u},v)} in the set of points to evaluate,
#' compute the basis functions
#' \deqn{
#' \varphi_{jt}(\bm{u},v) =
#' \left[ 2 - \frac{\Vert \bm{u} - \bm{c}_j \Vert^2}{w_s^2}- \frac{|v - g_t|^2}{w_t^2} \right]^2  \cdot
#' I(\Vert \bm{u} - \bm{c}_j \Vert \leq w_s) \cdot
#' I(|v - g_t| \leq w_t).
#' }
#' for \eqn{j = 1, \ldots, R} and \eqn{t = 1, \ldots, T}.
#' 
#' If \code{knots} is interpreted as a matrix, the three columns correspond
#' to x-axis and y-axis coordinates. Here, it is assumed that \code{dom}
#' and \code{sf} are based on a common coordinate system.
#' 
#' @examples
#' set.seed(1234)
#' seq_x = seq(0, 1, length.out = 3)
#' seq_y = seq(0, 1, length.out = 3)
#' seq_t = seq(0, 1, length.out = 3)
#' knots = expand.grid(seq_x, seq_y, seq_t)
#' x = runif(50)
#' y = runif(50)
#' t = sample(1:3, size = 50, replace = TRUE)
#' rad = 0.5
#' 
#' spacetime_bisquare(cbind(x,y,t), knots, w_s = rad, w_t = 1)
#' 
#' # Plot the (spatial) knots and the points at which we evaluated the basis
#' plot(knots[,1], knots[,2], pch = 4, cex = 1.5, col = "red")
#' text(x, y, labels = t, cex = 0.75)
#' 
#' # Draw a circle representing the basis' radius around one of the knot points
#' tseq = seq(0, 2*pi, length=100) 
#' coords = cbind(rad * cos(tseq) + seq_x[2], rad * sin(tseq) + seq_y[2])
#' lines(coords, col = "red")
#' 
#' @export
spacetime_bisquare = function(dom, knots_s, knots_t, w_s, w_t)
{
	out = prepare_bisquare(dom, knots_s, knots_t, type = "point")
	compute_basis_spt(out$X, out$knot_mat, w_s, w_t)
}
