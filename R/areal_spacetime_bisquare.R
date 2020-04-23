#' Areal SpaceTime Bisquare Basis
#' 
#' An \code{R6Class} representing the space-time bisquare basis.
#' 
#' @section Usage:
#' \preformatted{
#' bs = ArealSpaceTimeBisquareBasis$new(knots_x, knots_y, knots_t,
#'     w_s, w_t, mc_reps)
#' bs$compute(dom, period)
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
#' \item \code{dom} an \code{sf} object; areal units to evaluate.
#' \item \code{period} numeric vector; time coordinates for points to evaluate.
#' }
#' 
#' @section Methods:
#' \itemize{
#' \item \code{new} Create a new \code{SpaceTimeBisquareBasis} object.
#' \item \code{get_dim} Get the number of cutpoints used to construct this basis.
#' \item \code{get_cutpoints} Get the cutpoints used to construct this basis.
#' \item \code{get_ws} Get the spatial radius used to construct this basis.
#' \item \code{get_wt} Get the temporal radius used to construct this basis.
#' \item \code{compute} Evaluate this basis on specific areal units and periods.
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
#' 
#' # Create a simple domain from rectangles
#' shape1 = matrix(c(0.0,0.0, 0.5,0.0, 0.5,0.5, 0.0,0.5, 0.0,0.0), ncol=2, byrow=TRUE)
#' shape2 = shape1 + cbind(rep(0.5,5), rep(0.0,5))
#' shape3 = shape1 + cbind(rep(0.0,5), rep(0.5,5))
#' shape4 = shape1 + cbind(rep(0.5,5), rep(0.5,5))
#' sfc = st_sfc(
#'    st_polygon(list(shape1)),
#'    st_polygon(list(shape2)),
#'    st_polygon(list(shape3)),
#'    st_polygon(list(shape4))
#' )
#' dom = st_sf(data.frame(geoid = 1:length(sfc), geom = sfc))
#' period = c(0.4, 0.7)
#' 
#' bs = ArealSpaceTimeBisquareBasis$new(knots[,1], knots[,2], knots[,3],
#'     w_s = 0.5, w_t = 1, mc_reps = 200)
#' bs$compute(dom, period)
#' bs$get_dim()
#' bs$get_cutpoints()
#' bs$get_ws()
#' bs$get_wt()
#' 
#' # Plot the (spatial) knots and the (spatial) domain at which we evaluated
#' # the basis.
#' plot(knots[,1], knots[,2], pch = 4, cex = 1.5, col = "red")
#' plot(dom[,1], col = NA, add = TRUE)
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
	lock_class = FALSE,
	private = list(
		method = NULL,
		mc_reps = NULL,
		nx = NULL,
		ny = NULL,
		report_period = NULL,
		basis_spt = NULL
	),
	public = list(
		initialize = function(knots_x, knots_y, knots_t, w_s, w_t) {
			# Default is Monte Carlo with 1000 reps
			self$set_monte_carlo(1000)
			private$report_period = Inf
			private$basis_spt = SpaceTimeBisquareBasis$new(knots_x, knots_y, knots_t, w_s = w_s, w_t = w_t)
		},
		set_monte_carlo = function(reps) {
			private$method = "MonteCarlo"
			private$mc_reps = reps
			private$nx = NULL
			private$ny = NULL
		},
		set_quad = function(nx, ny) {
			private$method = "Quadrature"
			private$mc_reps = NULL
			private$nx = nx
			private$ny = ny
		},
		set_w = function(w_s, w_t) {
			knots = private$basis_spt$get_cutpoints()
			private$basis_spt = SpaceTimeBisquareBasis$new(knots[,1], knots[,2],
				knots[,3], w_s = w_s, w_t = w_t)
		},
		set_report_period = function(report_period) {
			private$report_period = report_period
		},
		get_method = function() {
			private$method
		},
		get_mc_reps = function() {
			private$mc_reps
		},
		get_nx = function() {
			private$nx
		},
		get_ny = function() {
			private$nx
		},
		get_basis_spt = function() {
			private$basis_spt
		},
		get_dim = function() {
			private$basis_spt$get_dim()
		},
		get_cutpoints = function() {
			private$basis_spt$get_cutpoints()
		},
		get_ws = function() {
			private$basis_spt$get_ws()
		},
		get_wt = function() {
			private$basis_spt$get_wt()
		},
		compute = function(dom, period) {
			if (private$method == "MonteCarlo") {
				S = private$compute_mc(dom, period)
			} else if (private$method == "Quadrature") {
				S = private$compute_quad(dom, period)
			} else {
				stop("Unknown method")
			}
			return(S)
		}
	)
)

ArealSpaceTimeBisquareBasis$set("private", "compute_mc", 
	function(dom, period) {
		bs = private$basis_spt
		R = private$mc_reps
		n = nrow(dom)
		r = bs$get_dim()
		S = Matrix(0, n, r)
		report_period = private$report_period

		for (j in 1:n) {
			if (j %% report_period == 0) {
				logger("Computing basis for area %d of %d\n", j, n)
			}

			# Monte Carlo approximation is done only over space and not period.
			# Therefore, it should be okay to reuse randomly drawn spatial points
			# for the same area across multiple periods.
			# 
			# The factor of 1.2 helps to achieve the desired sample size
			# in rdomain without repeating the loop there.
			P = rdomain(R, dom[j,], blocksize = ceiling(1.2*R), itmax = R)

			for (t in seq_along(period)) {
				B = bs$compute(P$x, P$y, period[t])
				S[j,] = S[j,] + colSums(B)
			}
		}

		return( S / (R * length(period)) )
	}
)

ArealSpaceTimeBisquareBasis$set("private", "compute_quad", 
	function(dom, period) {
		bs = private$basis_spt
		n = nrow(dom)
		r = bs$get_dim()
		S = Matrix(0, n, r)
		nx = private$nx
		ny = private$ny
		report_period = private$report_period

		for (j in 1:n) {
			if (j %% report_period == 0) {
				logger("Computing basis for area %d of %d\n", j, n)
			}

			out = make_grid(dom[j,], nx, ny)
			X = out$X
			dx = out$dx
			dy = out$dy

			area = as.numeric(st_area(dom[j,]))
			for (t in seq_along(period)) {
				B = bs$compute(X[,1], X[,2], period[t])
				S[j,] = S[j,] + colSums(B) * dx * dy / area
			}
		}
		return(S / length(period))			
	}
)

ArealSpaceTimeBisquareBasis$lock()
