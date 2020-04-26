#' Areal Spatial Bisquare Basis
#' 
#' @description
#' Spatial bisquare basis on areal data.
#' 
#' @param dom An \code{sf} or \code{sfc} object with areas
#' \eqn{A_1, \ldots, A_n} to evaluate.
#' @param knots Knots \eqn{\bm{c}_1, \ldots, \bm{c}_R} for the basis.
#' See "Details".
#' @param w Radius for the basis.
#' @param control A \code{list} of control arguments. See "Details".
#'
#' @return A sparse \eqn{n \times R} matrix whose \eqn{i}th row
#' is
#' \eqn{
#' \bm{s}_i^\top =
#' \Big(
#' \psi_1(A_i), \ldots, \psi_R(A_i)
#' \Big).
#' }
#' 
#' @details
#' \code{knots} may be provided as either an \code{sf} or \code{sfc} object, or as a
#' matrix of points.
#' \itemize{
#' \item If an \code{sf} or \code{sfc} object is provided for \code{knots}, \eqn{R}
#'   two-dimensional \code{POINT} entries are expected in \code{st_geometry(knots)}.
#'   Otherwise, \code{knots} will be interpreted as an \eqn{R \times 2} numeric matrix.
#' }
#' If \code{knots} is an \code{sf} or \code{sfc} object, it is checked
#' to ensure the coordinate system matches \code{dom}.
#' 
#' For each area \eqn{A} in the given domain, compute an approximation to
#' the basis functions
#' \deqn{
#' \psi_j(A) = \frac{1}{|A|} \int_A \varphi_j(\bm{u}) d\bm{u}
#' }
#' for \eqn{j = 1, \ldots, R}. Here, \eqn{\varphi_j(\bm{u})} represent
#' \link{spatial_bisquare} basis functions defined at the point level
#' using \eqn{\bm{c}_j} and \eqn{w}.
#' 
#' The \code{control} argument is a list which may provide any of the following:
#' \itemize{
#' \item \code{method} specifies computation method: use \code{"mc"} for
#' Monte Carlo or \code{"quad"} for quadrature. Default is \code{"mc"}.
#' \item \code{mc_reps} is number of repetitions to use for Monte Carlo.
#' Default is 1000.
#' \item \code{nx} is number of x-axis grid points to use for quadrature
#' method. Default is 50.
#' \item \code{ny} is number of y-axis grid points to use for quadrature
#' method. Default is 50.
#' \item \code{report_period} is an integer; print a message with progress each
#' time this many areas are processed. Default is \code{Inf} so that message
#' is suppressed.
#' \item \code{verbose} is a logical; if \code{TRUE} print descriptive
#' messages about the computation. Default is \code{FALSE}.
#' \item \code{mc_sampling_factor} is a positive number; an oversampling factor
#' used to compute \code{blocksize} in the \link{rdomain} function. I.e.,
#' \code{blocksize = ceiling(mc_sampling_factor * mc_reps)}. Default
#' is 1.2.
#' }
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
#' # Create a simple domain (of rectangles) to evaluate
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
#' 
#' rad = 0.5
#' areal_spatial_bisquare(dom, knots, rad)
#' areal_spatial_bisquare(dom, knots_sf, rad)
#' 
#' # Plot the knots and the points at which we evaluated the basis
#' plot(knots[,1], knots[,2], pch = 4, cex = 1.5, col = "red")
#' plot(dom[,1], col = NA, add = TRUE)
#' 
#' # Draw a circle representing the basis' radius around one of the knot points
#' tseq = seq(0, 2*pi, length=100) 
#' coords = cbind(rad * cos(tseq) + seq_x[2], rad * sin(tseq) + seq_y[2])
#' lines(coords, col = "red")
#' 
#' @export
areal_spatial_bisquare = function(dom, knots, w, control = NULL)
{
	out = prepare_bisquare(dom, knots, knots_t = NULL, type = "areal")
	X = out$X
	knot_mat = out$knot_mat
	R = out$R

	n = nrow(dom)
	S = Matrix(0, n, R)

	if (is.null(control)) { control = list() }
	if (is.null(control$mc_reps)) { control$mc_reps = 1000 }
	if (is.null(control$nx)) { control$nx = 50 }
	if (is.null(control$ny)) { control$ny = 50 }
	if (is.null(control$report_period)) { control$report_period = n + 1 }
	if (is.null(control$verbose)) { control$verbose = FALSE }
	if (is.null(control$method)) { control$method = "mc" }
	if (is.null(control$mc_sampling_factor)) { control$mc_sampling_factor = 1.2 }

	reps = control$mc_reps
	nx = control$nx
	ny = control$ny
	report_period = control$report_period
	verbose = control$verbose
	method = control$method
	blocksize = ceiling(control$mc_sampling_factor * reps)

	if (verbose && method == "mc") {
		printf("Using Monte Carlo method with %d reps\n", reps)
		printf("Sampling block size is %d\n", blocksize)
	} else if (verbose && method == "quad") {
		printf("Using quadrature method with %d x %d grid\n", nx, ny)
	}

	if (verbose) {
		printf("Computing %d areas using %d knots\n", n, R)
		printf("Radius w = %g\n", w)
	}

	for (j in 1:n) {
		if (j %% report_period == 0) {
			logger("Computing basis for area %d of %d\n", j, n)
		}

		if (method == "mc") {
			# The blocksize factor of helps to achieve the desired sample size
			# in rdomain without repeating the loop there.
			P = rdomain(reps, dom[j,], blocksize = blocksize, itmax = reps)
			B = compute_basis_sp(cbind(P$x,P$y), knot_mat, w)
			S[j,] = colSums(B) / reps
		} else if (method == "quad") {
			grid_out = make_grid(dom[j,], nx, ny)
			B = compute_basis_sp(grid_out$X, knot_mat, w)
			area = as.numeric(st_area(dom[j,]))
			S[j,] = colSums(B) * grid_out$dx * grid_out$dy / area
		} else {
			stop("Unrecongnized method. Use `mc` for Monte Carlo or 'quad' for quadrature")
		}
	}

	return(S)
}
