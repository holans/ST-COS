#' Areal SpaceTime Bisquare Basis
#' 
#' An \code{\link{R6Class}} representing the space-time bisquare basis.
#' 
#' @section Usage:
#' \preformatted{
#' basis = SpaceTimeBisquareBasis$new(knots.x, knots.y, knots.t,
#'     w.s, w.t)
#' basis$compute(x, y, time)
#' basis$get_dim()
#' basis$get_cutpoints()
#' basis$get_rl()
#' basis$get_ws()
#' basis$get_wt()
#' }
#' 
#' @section Arguments:
#' \itemize{
#' \item \code{knots.x} x-coordinate of knot points.
#' \item \code{knots.y} y-coordinate of knot points.
#' \item \code{knots.t} time coordinate of knot points.
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
#' set.seed(1234)
#' seq_x = seq(0, 1, length.out = 3)
#' seq_y = seq(0, 1, length.out = 3)
#' seq_t = seq(0, 1, length.out = 3)
#' knots = expand.grid(seq_x, seq_y, seq_t)
#' x = runif(50)
#' y = runif(50)
#' t = sample(1:3, size = 50, replace = TRUE)
#' 
#' basis = SpaceTimeBisquareBasis$new(knots[,1], knots[,2], knots[,3], w.s = 0.5, w.t = 1)
#' basis$compute(x, y, t)
#' basis$get_dim()
#' basis$get_cutpoints()
#' basis$get_rl()
#' basis$get_ws()
#' basis$get_wt()
#' 
#' # Plot the (spatial) knots and the points at which we evaluated the basis
#' plot(knots[,1], knots[,2], pch = 4, cex = 1.5, col = "red")
#' text(x, y, labels = t, cex = 0.75)
#' 
#' # Draw a circle representing the basis' radius around one of the knot points
#' tseq = seq(0, 2*pi, length=100) 
#' rad = basis$get_rl()
#' coords = cbind(rad * cos(tseq) + seq_x[2], rad * sin(tseq) + seq_y[2])
#' lines(coords, col = "red")
NULL

#' @export
#' @docType class
ArealSpaceTimeBisquareBasis = R6Class("ArealSpaceTimeBisquareBasis",
	lock_objects = TRUE,
	lock_class = TRUE,
	private = list(
		mc_reps = NULL,
		report_period = NULL,
		basis_spt = NULL
	),
	public = list(
		initialize = function(knots_x, knots_y, knots_t, w_s, w_t, mc_reps, report_period = mc_reps + 1) {
			private$mc_reps = mc_reps
			private$report_period = report_period
			private$basis_spt = SpaceTimeBisquareBasis$new(knots_x, knots_y, knots_t, w_s = w_s, w_t = w_t)
		},
		get_mc_reps = function() {
			private$mc_reps
		},
		get_report_period = function() {
			private$report_period
		},
		get_basis_spt = function() {
			private$basis_spt
		},
		compute = function(dom, period) {
			# X = cbind(x, y, time)
			# S = compute_basis_spt(X, private$cutpoints, private$rl, private$w.t)
			# return(S)
			
			basis = private$basis_spt
			R = private$mc_reps

			n = nrow(dom)
			r = basis$get_dim()
			S = Matrix(0, n, r)
			T = length(period)
			report_period = private$report_period

			for (j in 1:n) {
				if (j %% report_period == 0) {
					logger("Computing basis for area %d of %d\n", j, n)
				}

				# Drawing samples from an area seems more time consuming than computing
				# basis function. Let's reuse samples over multiple lookbacks.
				# Request a few more samples than we'll need, to prevent the loop in rArea.
				P = rArea(R, dom[j,], blocksize = ceiling(1.2*R))

				for (t in 1:T) {
					S[j,] = S[j,] + colSums(basis$compute(P[,1], P[,2], period[t]))
				}
			}

			return( S / (R*T) )
		}
	)
)
