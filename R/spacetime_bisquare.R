#' Space-Time Bisquare Basis
#' 
#' @description
#' Space-time bisquare basis on point data.
#' 
#' @param dom Space-time points \eqn{(\bm{u}_1,v_1), \ldots, (\bm{u}_n,v_n)}
#' to evaluate. See "Details".
#' @param knots Spatio-temporal knots
#' \eqn{(\bm{c}_1,g_1), \ldots, (\bm{c}_r,g_r)}
#' for the basis. See "Details".
#' @param w_s Spatial radius for the basis.
#' @param w_t Temporal radius for the basis.
#'
#' @return A sparse \eqn{n \times r} matrix whose \eqn{i}th row
#' is
#' \deqn{
#' \bm{s}_i^\top =
#' \Big(
#' \psi_1(\bm{u}_i,v_i), \ldots, \psi_r(\bm{u}_i,v_i)
#' \Big).
#' }
#' 
#' @details
#' Notes about arguments:
#' \itemize{
#' \item Both \code{dom} and \code{knots} may be provided as either \code{sf} or
#'   \code{sfc} objects, or as matrices of points.
#' \item If an \code{sf} or \code{sfc} object is provided for \code{dom}, \eqn{n}
#'   three-dimensional \code{POINT} entries are expected in \code{st_geometry(dom)}.
#'   Otherwise, \code{dom} will be interpreted as an \eqn{n \times 3} numeric matrix.
#' \item If an \code{sf} or \code{sfc} object is provided for \code{knots}, \eqn{r}
#'   three-dimensional \code{POINT} entries are expected in \code{st_geometry(knots)}.
#'   Otherwise, \code{knots} will be interpreted as an \eqn{r \times 3} numeric matrix.
#' \item If both \code{dom} and \code{knots_s} are given as \code{sf} or \code{sfc} objects,
#'   they will be checked to ensure a common coordinate system.
#' }
#' 
#' For each \eqn{(\bm{u}_i,v_i)}, compute the basis functions
#' \deqn{
#' \psi_j(\bm{u},v) =
#' \left[ 2 - \frac{\Vert \bm{u} - \bm{c}_j \Vert^2}{w_s^2}- \frac{|v - g_j|^2}{w_t^2} \right]^2  \cdot
#' I(\Vert \bm{u} - \bm{c}_j \Vert \leq w_s) \cdot
#' I(|v - g_j| \leq w_t)
#' }
#' for \eqn{j = 1, \ldots, r}.
#' 
#' Due to the treatment of \eqn{\bm{u}_i} and \eqn{\bm{c}_j} as points in a
#' Euclidean space, this basis is more suitable for coordinates from a map
#' projection than coordinates based on a globe representation.
#' 
#' @examples
#' set.seed(1234)
#' 
#' # Create knot points
#' seq_x = seq(0, 1, length.out = 3)
#' seq_y = seq(0, 1, length.out = 3)
#' seq_t = seq(0, 1, length.out = 3)
#' knots = expand.grid(x = seq_x, y = seq_y, t = seq_t)
#' knots_sf = st_as_sf(knots, coords = c("x","y","t"), crs = NA, dim = "XYM", agr = "constant")
#' 
#' # Points to evaluate
#' x = runif(50)
#' y = runif(50)
#' t = sample(1:3, size = 50, replace = TRUE)
#' pts = data.frame(x = x, y = y, t = t)
#' dom = st_as_sf(pts, coords = c("x","y","t"), crs = NA, dim = "XYM", agr = "constant")
#' 
#' rad = 0.5
#' spacetime_bisquare(cbind(x,y,t), knots, w_s = rad, w_t = 1)
#' spacetime_bisquare(dom, knots_sf, w_s = rad, w_t = 1)
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
#' @family bisquare
#' @export
spacetime_bisquare = function(dom, knots, w_s, w_t)
{
	prep = prepare_bisquare(dom, knots, type = "point")
	out = compute_basis_spt(prep$X, prep$knot_mat, w_s, w_t)
	sparseMatrix(i = out$ind_row + 1, j = out$ind_col + 1, x = out$values,
		dims = out$dim)
}
