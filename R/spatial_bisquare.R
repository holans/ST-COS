#' Spatial Bisquare Basis
#' 
#' @description
#' Spatial bisquare basis on point data.
#' 
#' @param dom Points \eqn{\bm{u}_1, \ldots, \bm{u}_n} to evaluate. See
#' "Details".
#' @param knots Knots \eqn{\bm{c}_1, \ldots, \bm{c}_R} for the basis.
#' See "Details".
#' @param w Radius for the basis.
#'
#' @return A sparse \eqn{n \times R} matrix whose \eqn{i}th row
#' is
#' \eqn{
#' \bm{s}_i^\top =
#' \Big(
#' \varphi_1(\bm{u}_i), \ldots, \varphi_R(\bm{u}_i)
#' \Big).
#' }
#'   
#' @details
#' Both \code{dom} and \code{knots} may be provided as either \code{sf} or
#' \code{sfc} objects, or as matrices of points.
#' \itemize{
#' \item If an \code{sf} or \code{sfc} object is provided for \code{dom}, \eqn{n}
#'   two-dimensional \code{POINT} entries are expected in \code{st_geometry(dom)}.
#'   Otherwise, \code{dom} will be interpreted as an \eqn{n \times 2} numeric matrix.
#' \item If an \code{sf} or \code{sfc} object is provided for \code{knots}, \eqn{R}
#'   two-dimensional \code{POINT} entries are expected in \code{st_geometry(knots)}.
#'   Otherwise, \code{knots} will be interpreted as an \eqn{R \times 2} numeric matrix.
#' }
#' If both \code{dom} and \code{knots} are given as \code{sf} or \code{sfc} objects,
#' they will be checked to ensure a common coordinate system.
#' 
#' For each \eqn{\bm{u}_i}, compute the basis functions
#' \deqn{
#' \varphi_j(\bm{u}) =
#' \left[ 1 - \frac{\Vert\bm{u} - \bm{c}_j \Vert^2}{w^2} \right]^2  \cdot
#' I(\Vert \bm{u} - \bm{c}_j \Vert \leq w)
#' }
#' for \eqn{j = 1, \ldots, R}.
#' 
#' @examples
#' set.seed(1234)
#' 
#' # Create knot points
#' seq_x = seq(0, 1, length.out = 3)
#' seq_y = seq(0, 1, length.out = 3)
#' knots = expand.grid(x = seq_x, y = seq_y)
#' knots_sf = st_as_sf(knots, coords = c("x","y"), crs = NA, agr = "constant")
#' 
#' # Points to evaluate
#' x = runif(50)
#' y = runif(50)
#' pts = data.frame(x = x, y = y)
#' dom = st_as_sf(pts, coords = c("x","y"), crs = NA, agr = "constant")
#' 
#' rad = 0.5
#' spatial_bisquare(cbind(x,y), knots, rad)
#' spatial_bisquare(dom, knots, rad)
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
spatial_bisquare = function(dom, knots, w)
{
	out = prepare_bisquare(dom, knots, type = "point")
	compute_basis_sp(out$X, out$knot_mat, w)
}
