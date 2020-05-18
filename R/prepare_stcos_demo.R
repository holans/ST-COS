#' Prepare Demo Data for STCOS Model
#' 
#' Create demo data based on ACS example, making a few simple model choices.
#' The purpose of this function is to facilitate examples in other functions.
#' Uses functions in the package to create model terms from shapefiles.
#' 
#' @param num_knots_sp Number of spatial knots to use in areal space-time
#'     basis.
#' @param basis_mc_reps Number of monte carlo reps to use in areal space-time
#'     basis.
#' @param eigval_prop Proportion of variability to keep in dimension reduction
#'     of basis expansions.
#' @return A list containing the following:
#' \itemize{
#' \item \code{z} direct estimates.
#' \item \code{v} direct variance estimates.
#' \item \code{H} overlap matrix.
#' \item \code{S} design matrix of basis expansion.
#' \item \code{K} covariance matrix of the random effect.
#' }
#' 
#' @examples
#' \dontrun{
#' out = prepare_stcos_demo()
#' }
#' @export
prepare_stcos_demo = function(num_knots_sp = 200, basis_mc_reps = 200, eigval_prop = 0.65)
{
	logger("[1/7] Preparing data\n")
	myenv = new.env()
	utils::data("acs_sf", envir = myenv)
	acs5_2013 = myenv$acs5_2013
	acs5_2014 = myenv$acs5_2014
	acs5_2015 = myenv$acs5_2015
	acs5_2016 = myenv$acs5_2016
	acs5_2017 = myenv$acs5_2017
	dom_fine = acs5_2017
	n = nrow(dom_fine)

	logger("[2/7] Preparing areal space-time basis function\n")
	pts = st_sample(dom_fine, num_knots_sp, type = "hexagonal")
	knots_sp = st_coordinates(pts)
	knots_t = seq(2009, 2017, by = 0.5)
	knots = merge(knots_sp, knots_t)

	D = stats::dist(knots)
	ws_tx = quantile(D[D > 0], prob = 0.05, type = 1)
	bs_ctrl = list(mc_reps = basis_mc_reps)

	logger("[3/7] Preparing overlap matrix\n")
	H = rbind(
	    overlap_matrix(acs5_2013, dom_fine),
		overlap_matrix(acs5_2014, dom_fine),
		overlap_matrix(acs5_2015, dom_fine),
		overlap_matrix(acs5_2016, dom_fine),
	    overlap_matrix(acs5_2017, dom_fine)
	)

	logger("[4/7] Computing areal basis on source supports\n")
	S_full = rbind(
		areal_spacetime_bisquare(acs5_2013, 2009:2013, knots, ws_tx, w_t = 1, control = bs_ctrl),
		areal_spacetime_bisquare(acs5_2014, 2010:2014, knots, ws_tx, w_t = 1, control = bs_ctrl),
		areal_spacetime_bisquare(acs5_2015, 2011:2015, knots, ws_tx, w_t = 1, control = bs_ctrl),
		areal_spacetime_bisquare(acs5_2016, 2012:2016, knots, ws_tx, w_t = 1, control = bs_ctrl),
		areal_spacetime_bisquare(acs5_2017, 2013:2017, knots, ws_tx, w_t = 1, control = bs_ctrl)
	)

	logger("[5/7] Computing areal basis on fine-level support\n")
	S_fine_full = rbind(
		areal_spacetime_bisquare(dom_fine, 2009, knots, ws_tx, w_t = 1, control = bs_ctrl),
		areal_spacetime_bisquare(dom_fine, 2010, knots, ws_tx, w_t = 1, control = bs_ctrl),
		areal_spacetime_bisquare(dom_fine, 2011, knots, ws_tx, w_t = 1, control = bs_ctrl),
		areal_spacetime_bisquare(dom_fine, 2012, knots, ws_tx, w_t = 1, control = bs_ctrl),
		areal_spacetime_bisquare(dom_fine, 2013, knots, ws_tx, w_t = 1, control = bs_ctrl),
		areal_spacetime_bisquare(dom_fine, 2014, knots, ws_tx, w_t = 1, control = bs_ctrl),
		areal_spacetime_bisquare(dom_fine, 2015, knots, ws_tx, w_t = 1, control = bs_ctrl),
		areal_spacetime_bisquare(dom_fine, 2016, knots, ws_tx, w_t = 1, control = bs_ctrl),
		areal_spacetime_bisquare(dom_fine, 2017, knots, ws_tx, w_t = 1, control = bs_ctrl)
	)

	z = c(acs5_2013$DirectEst, acs5_2014$DirectEst, acs5_2015$DirectEst,
		acs5_2016$DirectEst, acs5_2017$DirectEst)
	v = c(acs5_2013$DirectVar, acs5_2014$DirectVar, acs5_2015$DirectVar,
		acs5_2016$DirectVar, acs5_2017$DirectVar)

	# A quick and dirty way to handle missing values so that the demo will run.
	# Not a good idea for a real data analysis...
	z[which(is.na(z))] = mean(z, na.rm = TRUE)
	v[which(is.na(v))] = mean(v, na.rm = TRUE)

	logger("[6/7] Dimension reduction\n")
	eig = eigen(t(S_full) %*% S_full)
	idx_S = which(cumsum(eig$values) / sum(eig$values) < eigval_prop)
	Ts = eig$vectors[,idx_S,drop=FALSE]

	S = S_full %*% Ts
	S_fine = S_fine_full %*% Ts

	logger("[7/7] Computing K matrix\n")
	W = adjacency_matrix(dom_fine)
	Q = car_precision(W, tau = 0.9, scale = TRUE)
	Qinv = solve(Q)

	K = cov_approx_blockdiag(Qinv, S_fine)
	logger("Done\n")

	list(z = z, v = v, H = H, S = S, K = K)
}
