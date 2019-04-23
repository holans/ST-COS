#' Spacial Bisquare Basis
#' 
#' An \code{\link{R6Class}} representing the spacial bisquare basis.
#' 
#' @section Usage:
#' \preformatted{
#' basis <- SpaceTimeBisquareBasis$new(knots.x, knots.y, w)
#' basis$compute(x, y, time)
#' basis$get_dim()
#' basis$get_cutpoints()
#' basis$get_rl()
#' basis$get_w()
#' }
#' 
#' @section Arguments:
#' \itemize{
#' \item \code{knots.x} x-coordinate of knot points.
#' \item \code{knots.y} y-coordinate of knot points.
#' \item \code{w} (Original, before transformation) radius.
#' \item \code{x} Vector of x-coordinates for points on which to evaluate the basis.
#' \item \code{y} Vector of y-coordinates for points on which to evaluate the basis.
#' }
#' 
#' @section Methods:
#' \itemize{
#' \item \code{new} Create a new \code{SpatialBisquareBasis} object.
#' \item \code{get_dim} Get the number of cutpoints used to construct this basis.
#' \item \code{get_cutpoints} Get the cutpoints used to construct this basis.
#' \item \code{get_rl} Get the transformed spatial radius. The transformation is
#' based on quantiles of distances between knots.
#' \item \code{get_w} Get the original radius used to construct this
#' basis. This is transformed before it is applied to account for the geography
#' being used.
#' \item \code{compute} Evaluate this basis on specific points.
#' }
#'
#' @name SpatialBisquareBasis
#' 
#' @examples
#' set.seed(1234)
#' seq.x <- seq(0, 1, length.out = 3)
#' seq.y <- seq(0, 1, length.out = 3)
#' knots = merge(seq.x, seq.y)
#' x <- runif(50)
#' y <- runif(50)
#' 
#' basis <- SpatialBisquareBasis$new(knots[,1], knots[,2], w = 0.5)
#' basis$compute(x, y)
#' basis$get_dim()
#' basis$get_cutpoints()
#' basis$get_rl()
#' basis$get_w()
#' 
#' # Plot the knots and the points at which we evaluated the basis
#' plot(knots[,1], knots[,2], pch = 4, cex = 1.5, col = "red")
#' points(x, y, cex = 0.5)
#' 
#' # Draw a circle representing the basis' radius around one of the knot points
#' tseq <- seq(0, 2*pi, length=100) 
#' rad <- basis$get_rl()
#' coords <- cbind(rad * cos(tseq) + seq.x[2], rad * sin(tseq) + seq.y[2])
#' lines(coords, col = "red")
NULL

#' @export
#' @docType class
SpatialBisquareBasis <- R6Class("SpatialBisquareBasis",
	lock_objects = TRUE,
	lock_class = TRUE,
	private = list(
		r = NULL,
		cutpoints = NULL,
		w = NULL,
		rl = NULL
	),
	public = list(
		initialize = function(knots.x, knots.y, w) {
			r <- length(knots.x)
			stopifnot(length(knots.y) == r)
			private$cutpoints <- cbind(knots.x, knots.y)
			private$w <- w
			private$r <- r

			# Jon's code computes basis with rl instead of w.s
			# Use type 1 quantile algorithm to match Matlab 
			G <- dist(private$cutpoints)
			private$rl <- as.numeric(w * quantile(G[G > 0], prob = 0.05, type = 1))
		},
		get_dim = function() {
			private$r
		},
		get_cutpoints = function() {
			private$cutpoints
		},
		get_rl = function() {
			private$rl
		},
		get_w = function() {
			private$w
		},
		compute = function(x, y) {
			X <- cbind(x, y)
			S <- compute_basis_sp(X, private$cutpoints, private$rl)
			return(S)
		}
	)
)
