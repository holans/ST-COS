#' Areal SpaceTime Bisquare Basis
#' 
#' An \code{R6Class} representing the space-time bisquare basis.
#' 
#' @section Usage:
#' \preformatted{
#' bs = ArealSpaceTimeBisquareBasis$new(knots_x, knots_y, knots_t,
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
#' \item \code{w_s} (Original, before transformation) spatial radius.
#' \item \code{w_t} Temporal radius.
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
#' @name ArealSpaceTimeBisquareBasis
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
ArealSpaceTimeBisquareBasis = R6Class("ArealSpaceTimeBisquareBasis",
	lock_objects = TRUE,
	lock_class = TRUE,
	private = list(
		mc_reps = NULL,
		report_period = NULL,
		basis_spt = NULL
	),
	public = list(
		initialize = function(knots_x, knots_y, knots_t, w_s, w_t, mc_reps) {
			private$mc_reps = mc_reps
			private$basis_spt = SpaceTimeBisquareBasis$new(knots_x, knots_y, knots_t, w_s = w_s, w_t = w_t)
		},
		get_mc_reps = function() {
			private$mc_reps
		},
		get_basis_spt = function() {
			private$basis_spt
		},
		get_dim = function() {
			private$basis_sp$r
		},
		get_cutpoints = function() {
			private$basis_spt$cutpoints
		},
		get_ws = function() {
			private$basis_spt$w_s
		},
		get_wt = function() {
			private$basis_spt$w_t
		},
		compute = function(dom, period, report_period = nrow(dom) + 1) {
			basis = private$basis_spt
			R = private$mc_reps
			n = nrow(dom)
			r = basis$get_dim()
			S = Matrix(0, n, r)
			T = length(period)

			for (j in 1:n) {
				if (j %% report_period == 0) {
					logger("Computing basis for area %d of %d\n", j, n)
				}

				# Drawing samples from an area seems more time consuming than computing
				# basis function. Let's reuse samples over multiple lookbacks.
				# Request a few more samples than we'll need, to prevent the loop in rdomain.
				P = rdomain(R, dom[j,], blocksize = ceiling(1.2*R), itmax = R)

				for (t in 1:T) {
					S[j,] = S[j,] + colSums(basis$compute(P[,1], P[,2], period[t]))
				}
			}

			return( S / (R*T) )
		}
	)
)