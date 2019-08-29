#' Spatial Bisquare Basis
#' 
#' An \code{R6Class} representing the spatial bisquare basis.
#' 
#' @section Usage:
#' \preformatted{
#' bs = SpatialBisquareBasis$new(knots_x, knots_y, w)
#' bs$compute(x, y, time)
#' bs$get_dim()
#' bs$get_cutpoints()
#' bs$get_w()
#' }
#' 
#' @section Arguments:
#' \itemize{
#' \item \code{knots_x} numeric vector; x-coordinates of knot points.
#' \item \code{knots_y} numeric vector; y-coordinates of knot points.
#' \item \code{w} numeric; radius for the basis.
#' \item \code{x} numeric vector; x-coordinates for points to evaluate.
#' \item \code{y} numeric vector; y-coordinates for points to evaluate.
#' }
#' 
#' @section Methods:
#' \itemize{
#' \item \code{new} Create a new \code{SpatialBisquareBasis} object.
#' \item \code{get_dim} Get the number of cutpoints used to construct this basis.
#' \item \code{get_cutpoints} Get the cutpoints used to construct this basis.
#' \item \code{get_w} Get the radius used to construct this basis. 
#' \item \code{compute} Evaluate this basis on specific points.
#' }
#'
#' @name SpatialBisquareBasis
#' 
#' @examples
#' set.seed(1234)
#' seq_x = seq(0, 1, length.out = 3)
#' seq_y = seq(0, 1, length.out = 3)
#' knots = merge(seq_x, seq_y)
#' x = runif(50)
#' y = runif(50)
#' 
#' bs = SpatialBisquareBasis$new(knots[,1], knots[,2], w = 0.5)
#' bs$compute(x, y)
#' bs$get_dim()
#' bs$get_cutpoints()
#' bs$get_w()
#' 
#' # Plot the knots and the points at which we evaluated the basis
#' plot(knots[,1], knots[,2], pch = 4, cex = 1.5, col = "red")
#' points(x, y, cex = 0.5)
#' 
#' # Draw a circle representing the basis' radius around one of the knot points
#' tseq = seq(0, 2*pi, length=100) 
#' rad = bs$get_w()
#' coords = cbind(rad * cos(tseq) + seq_x[2], rad * sin(tseq) + seq_y[2])
#' lines(coords, col = "red")
NULL

#' @export
#' @docType class
SpatialBisquareBasis = R6Class("SpatialBisquareBasis",
	lock_objects = TRUE,
	lock_class = TRUE,
	private = list(
		r = NULL,
		cutpoints = NULL,
		w = NULL
	),
	public = list(
		initialize = function(knots_x, knots_y, w) {
			r = length(knots_x)
			stopifnot(length(knots_y) == r)
			private$cutpoints = cbind(knots_x, knots_y)
			private$w = w
			private$r = r
		},
		get_dim = function() {
			private$r
		},
		get_cutpoints = function() {
			private$cutpoints
		},
		get_w = function() {
			private$w
		},
		compute = function(x, y) {
			X = cbind(x, y)
			S = compute_basis_sp(X, private$cutpoints, private$w)
			return(S)
		}
	)
)
