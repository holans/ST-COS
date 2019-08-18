#' @export
printf = function(msg, ...)
{
	cat(sprintf(msg, ...))
}

#' @export
logger = function(msg, ...)
{
	sys.time = as.character(Sys.time())
	cat(sys.time, "-", sprintf(msg, ...))
}

#' @export
normalize = function(x)
{
	if (sum(x < 0) > 0) {
		stop("Cannot normalize a vector with negative elements")
	}

	s = sum(x)
	if (s > 0) {
		return(x / sum(x))
	}

	return(x)
}

# Copied from invgamma package. We may not need this...
dinvgamma = function (x, shape, rate, scale = 1/rate, log = FALSE) 
{
	if (missing(rate) && !missing(scale)) {
		rate = 1/scale
	}
	log_f = dgamma(1/x, shape, rate, log = TRUE) - 2 * log(x)
	if (log) {
		return(log_f)
	}
	exp(log_f)
}

# Code to compute an inverse for a non-pd matrix
# Do we still need this?
#my_inverse = function(Sigma, S) {
#	eig = eigen(Sigma)
#	P = Re(eig$vectors)
#	D = Re(eig$values)
#	D[D < 0] = 0
#	Dinv = D
#	Dinv[D > 0] = 1 / D[D > 0]
#	Kinv = P %*% (Dinv * t(P))
#	return(Kinv)
#}
