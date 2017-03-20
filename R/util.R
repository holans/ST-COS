printf <- function(msg, ...)
{
	cat(sprintf(msg, ...))
}

logger <- function(msg, ...)
{
	sys.time <- as.character(Sys.time())
	cat(sys.time, "-", sprintf(msg, ...))
}

normalize <- function(x)
{
	if (sum(x < 0) > 0) {
		stop("Cannot normalize a vector with negative elements")
	}

	s <- sum(x)
	if (s > 0) {
		return(x / sum(x))
	}

	return(x)
}

# Copied from invgamma package
dinvgamma <- function (x, shape, rate, scale = 1/rate, log = FALSE) 
{
	if (missing(rate) && !missing(scale)) 
		rate <- 1/scale
	log_f <- dgamma(1/x, shape, rate, log = TRUE) - 2 * log(x)
	if (log) 
		return(log_f)
	exp(log_f)
}
