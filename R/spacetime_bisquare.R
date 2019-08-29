#' SpaceTime Bisquare Basis
#' 
#' An \code{R6Class} representing the space-time bisquare basis.
#' 
#' @section Usage:
#' \preformatted{
#' bs = SpaceTimeBisquareBasis$new(knots_x, knots_y, knots_t,
#'     w_s, w_t)
#' bs$compute(x, y, time)
#' bs$get_dim()
#' bs$get_cutpoints()
#' bs$get_ws()
#' bs$get_wt()
#' }
#' 
#' @section Arguments:
#' \itemize{
#' \item \code{knots_x} numeric vector; x-coordinates of knot points.
#' \item \code{knots_y} numeric vector; y-coordinates of knot points.
#' \item \code{knots_t} numeric vector; time coordinate of knot points.
#' \item \code{w_s} numeric; spatial radius for the basis.
#' \item \code{w_t} numeric; temporal radius for the basis.
#' \item \code{x} numeric vector; x-coordinates for points to evaluate.
#' \item \code{y} numeric vector; y-coordinates for points to evaluate.
#' \item \code{time} numeric vector; time coordinates for points to evaluate.
#' }
#' 
#' @section Methods:
#' \itemize{
#' \item \code{new} Create a new \code{SpaceTimeBisquareBasis} object.
#' \item \code{get_dim} Get the number of cutpoints used to construct this basis.
#' \item \code{get_cutpoints} Get the cutpoints used to construct this basis.
#' \item \code{get_ws} Get the spatial radius used to construct this basis.
#' \item \code{get_wt} Get the temporal radius used to construct this basis.
#' \item \code{compute} Evaluate this basis on specific points.
#' }
#'
#' @name SpaceTimeBisquareBasis
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
#' 
#' bs = SpaceTimeBisquareBasis$new(knots[,1], knots[,2], knots[,3], w_s = 0.5, w_t = 1)
#' bs$compute(x, y, t)
#' bs$get_dim()
#' bs$get_cutpoints()
#' bs$get_ws()
#' bs$get_wt()
#' 
#' # Plot the (spatial) knots and the points at which we evaluated the basis
#' plot(knots[,1], knots[,2], pch = 4, cex = 1.5, col = "red")
#' text(x, y, labels = t, cex = 0.75)
#' 
#' # Draw a circle representing the basis' radius around one of the knot points
#' tseq = seq(0, 2*pi, length=100) 
#' rad = bs$get_ws()
#' coords = cbind(rad * cos(tseq) + seq_x[2], rad * sin(tseq) + seq_y[2])
#' lines(coords, col = "red")
NULL

#' @export
#' @docType class
SpaceTimeBisquareBasis = R6Class("SpaceTimeBisquareBasis",
	lock_objects = TRUE,
	lock_class = TRUE,
	private = list(
		r = NULL,
		cutpoints = NULL,
		w_s = NULL,
		w_t = NULL
	),
	public = list(
		initialize = function(knots_x, knots_y, knots_t, w_s, w_t) {
			r = length(knots_x)
			stopifnot(length(knots_y) == r)
			stopifnot(length(knots_t) == r)
			private$cutpoints = cbind(knots_x, knots_y, knots_t)
			private$w_s = w_s
			private$w_t = w_t
			private$r = r
		},
		get_dim = function() {
			private$r
		},
		get_cutpoints = function() {
			private$cutpoints
		},
		get_ws = function() {
			private$w_s
		},
		get_wt = function() {
			private$w_t
		},
		compute = function(x, y, time) {
			X = cbind(x, y, time)
			S = compute_basis_spt(X, private$cutpoints, private$w_s, private$w_t)
			return(S)
		}
	)
)
