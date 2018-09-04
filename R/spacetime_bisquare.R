#' SpaceTime Bisquare Basis
#' 
#' An \code{\link{R6Class}} representing the space-time bisquare basis.
#' 
#' @section Arguments:
#' \itemize{
#' \item \code{cutpoints.x} x-coordinate of knot points.
#' \item \code{cutpoints.y} y-coordinate of knot points.
#' \item \code{cutpoints.t} time coordinate of knot points.
#' \item \code{w.s} (Original, before transformation) spatial radius.
#' \item \code{w.t} Temporal radius.
#' \item \code{x} Vector of x-coordinates for points on which to evaluate the basis.
#' \item \code{y} Vector of y-coordinates for points on which to evaluate the basis.
#' \item \code{time} Vector of time coordinates for points on which to evaluate the basis.
#' }
#' 
#' @section Methods:
#' \itemize{
#' \item \code{new} Create a new \code{SpaceTimeBisquareBasis} object.
#' \item \code{get_dim} Get the number of cutpoints used to construct this basis.
#' \item \code{get_cutpoints} Get the cutpoints used to construct this basis.
#' \item \code{get_rl} Get the transformed spatial radius. The transformation is
#' based on quantiles of distances between knots.
#' \item \code{get_ws} Get the original spatial radius used to construct this
#' basis. This is transformed before it is applied to account for the geography
#' being used.
#' \item \code{get_wt} Get the temporal radius used to construct this basis.
#' \item \code{compute} Evaluate this basis on specific points.
#' }
#'
#' @name SpaceTimeBisquareBasis
#' 
#' @examples
#' \dontrun{
#' basis <- SpaceTimeBisquareBasis$new(cutpoints.x, cutpoints.y, cutpoints.t,
#'     w.s = 1, w.t = 1)
#' basis$compute(x, y, time)
#' basis$get_dim()
#' basis$get_cutpoints()
#' basis$get_rl()
#' basis$get_ws()
#' basis$get_wt()
#' }
NULL

#' @export
#' @docType class
SpaceTimeBisquareBasis <- R6Class("SpaceTimeBisquareBasis",
	private = list(
		r = NULL,
		cutpoints = NULL,
		w.s = NULL,
		w.t = NULL,
		rl = NULL
	),
	public = list(
		initialize = function(cutpoints.x, cutpoints.y, cutpoints.t, w.s, w.t) {
			r <- length(cutpoints.x)
			stopifnot(length(cutpoints.y) == r)
			stopifnot(length(cutpoints.t) == r)
			private$cutpoints <- cbind(cutpoints.x, cutpoints.y, cutpoints.t)
			private$w.s <- w.s
			private$w.t <- w.t
			private$r <- r

			# Jon's code computes basis with rl instead of w.s
			# Use type 1 quantile algorithm to match Matlab 
			G <- dist(private$cutpoints)
			private$rl <- as.numeric(w.s * quantile(G[G > 0], prob = 0.05, type = 1))
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
		get_ws = function() {
			private$w.s
		},
		get_wt = function() {
			private$w.t
		},
		compute = function(x, y, time) {
			X <- cbind(x, y, time)
			S <- compute_basis_spt(X, private$cutpoints, private$rl, private$w.t)
			return(S)
		}
	)
)

# SpaceTimeBisquareBasis$set("public", "compute", compute)
# SpaceTimeBisquareBasis$lock()
