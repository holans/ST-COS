#' Spacial Bisquare Basis
#' 
#' An \code{\link{R6Class}} representing the spacial bisquare basis.
#' 
#' @section Usage:
#' \preformatted{
#' basis <- SpaceTimeBisquareBasis$new(cutpoints.x, cutpoints.y, w)
#' basis$compute(x, y, time)
#' basis$get_dim()
#' basis$get_cutpoints()
#' basis$get_rl()
#' basis$get_w()
#' }
#' 
#' @section Arguments:
#' \itemize{
#' \item \code{cutpoints.x} x-coordinate of knot points.
#' \item \code{cutpoints.y} y-coordinate of knot points.
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
#' \dontrun{
#' basis <- SpaceTimeBisquareBasis$new(cutpoints.x, cutpoints.y, w = 1)
#' basis$compute(x, y, time)
#' basis$get_dim()
#' basis$get_cutpoints()
#' basis$get_rl()
#' basis$get_w()
#' }
NULL

#' @export
#' @docType class
SpatialBisquareBasis <- R6Class("SpatialBisquareBasis",
	private = list(
		r = NULL,
		cutpoints = NULL,
		w = NULL,
		rl = NULL
	),
	public = list(
		initialize = function(cutpoints.x, cutpoints.y, w) {
			r <- length(cutpoints.x)
			stopifnot(length(cutpoints.y) == r)
			private$cutpoints <- cbind(cutpoints.x, cutpoints.y)
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

# SpatialBisquareBasis$set("public", "compute", compute)
# SpatialBisquareBasis$lock()

