#' Areal Space-Time Bisquare Basis
#' 
#' @description
#' Space-Time bisquare basis on areal data.
#' 
#' @param dom An \code{sf} object with \eqn{n} areas to evaluate.
#' @param period A numeric vector of time periods \eqn{v_1, \ldots, v_m}
#' to evaluate for each area.
#' @param knots_s Spatial knots \eqn{\bm{c}_1, \ldots, \bm{c}_R} for the
#' basis. If an \code{sf} object is provided, it should contain \eqn{R}
#' \code{POINT} entries. Otherwise, it will be interpreted as an
#' \eqn{R \times 2} numeric matrix.
#' @param knots_t A numeric vector with temporal knots \eqn{g_1, \ldots, g_T}
#' for the basis.
#' @param w Radius for the basis.
#' @param control A \code{list} of control arguments. See "Details".
#'
#' @return A sparse \eqn{n \times RT} matrix whose \eqn{i}th row
#' represents the \eqn{i}th point \eqn{(\bm{u}_i,v_i)} evaluated at every
#' basis function for \eqn{j = 1, \ldots, R} and \eqn{t = 1, \ldots, T}.
#'   
#' @details
#' For each area \eqn{A} in the given domain and \eqn{v} in the given period,
#' compute an approximation to the basis functions
#' \deqn{
#' \psi_{jt}(A) = \frac{1}{m} \sum_{k=1}^m \int_A \varphi_{jt}(\bm{u},v_k) d\bm{u},
#' }
#' for \eqn{j = 1, \ldots, R} and \eqn{t = 1, \ldots, T}. Here,
#' \eqn{\varphi_{jt}{(\bm{u},v)}}
#' represent \link{spacetime_bisquare} basis functions defined at the point
#' level using \eqn{\bm{c}_j}, \eqn{g_t}, \eqn{w_s}, and \eqn{w_t}.
#' 
#' If \code{knots} is interpreted as a matrix, the three columns correspond
#' to x-axis and y-axis coordinates. Here, it is assumed that \code{dom}
#' and \code{sf} are based on a common coordinate system.
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
#' seq_x = seq(0, 1, length.out = 3)
#' seq_y = seq(0, 1, length.out = 3)
#' knots_s = expand.grid(x = seq_x, y = seq_y)#' 
#' knots_sf = st_as_sf(knots_s, coords = c("x","y"), crs = NA, agr = "constant")
#' knots_t = seq(0, 1, length.out = 3)

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
#' 
#' rad = 0.5
#' period = c(0.4, 0.7)
#' areal_spatial_bisquare(dom, period, knots, w_s = rad, w_t = 1)
#' areal_spatial_bisquare(dom, period, knots_sf, w_s = rad, w_t = 1)
#' 
#' # Plot the (spatial) knots and the (spatial) domain at which we evaluated
#' # the basis.
#' plot(knots[,1], knots[,2], pch = 4, cex = 1.5, col = "red")
#' plot(dom[,1], col = NA, add = TRUE)
#' 
#' # Draw a circle representing the basis' radius around one of the knot points
#' tseq = seq(0, 2*pi, length=100) 
#' coords = cbind(rad * cos(tseq) + seq_x[2], rad * sin(tseq) + seq_y[2])
#' lines(coords, col = "red")
#' 
#' @export
areal_spacetime_bisquare = function(dom, period, knots_s, knots_t, w_s, w_t, control = NULL)
{
	if ("sf" %in% class(knots_s)) {
		stopifnot(st_crs(dom) == st_crs(knots_s))
		knots_s = matrix(unlist(knots_s$geometry), ncol = 2, byrow = TRUE)
	} else {
		knots_s = as.matrix(knots_s)
	}

	n = nrow(dom)
	R = nrow(knots_s)
	T = length(knots_t)
	S = Matrix(0, n, R*T)
	m = length(period)

	knots = cbind(knots_s %x% matrix(1,T,1), matrix(1,R,1) %x% knots_t)

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
		period_str = ifelse(m > 4,
			paste(c(head(period,2), "...", tail(period,2)), collapse = ", "),
			paste(period, collapse = ", "))
		printf("Computing %d areas and %d periods\n", n, m)
		printf("Periods: %s\n", period_str)
		printf("Using %d spatial knots and %d temporal knots\n", R, T)
		printf("Spatial radius w_s = %g, temporal radius w_t = %g\n", w_s, w_t)
	}

	for (j in 1:n) {
		if (j %% report_period == 0) {
			logger("Computing basis for area %d of %d\n", j, n)
		}

		s_j = numeric(R*T)

		if (method == "mc") {
			# The blocksize factor of helps to achieve the desired sample size
			# in rdomain without repeating the loop there.
			P = rdomain(reps, dom[j,], blocksize = blocksize, itmax = reps)
			for (t in seq_along(period)) {
				X = cbind(P$x, P$y, period[t])
				B = compute_basis_spt(X, knots, w_s, w_t)
				s_j = s_j + colSums(B) / reps
			}
		} else if (method == "quad") {
			grid_out = make_grid(dom[j,], nx, ny)
			area = as.numeric(st_area(dom[j,]))
			for (t in seq_along(period)) {
				X = cbind(grid_out$X, period[t])
				B = compute_basis_spt(X, knots, w_s, w_t)
				s_j = s_j + colSums(B) * grid_out$dx * grid_out$dy / area
			}
		} else {
			stop("Unrecongnized method. Use `mc` for Monte Carlo or 'quad' for quadrature")
		}

		S[j,] = s_j
	}

	return(S / m)
}
