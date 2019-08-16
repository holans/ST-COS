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

# Convert an adjacency list into a (sparse) matrix.
#' @export
adjList2Matrix = function(a)
{
	n.tuples = sum(unlist(Map(length, a)))
	i = integer(n.tuples)
	j = integer(n.tuples)
	x = rep(1, n.tuples)
	idx.last = 0
	L = length(a)
	for (h in 1:L) {
		idx.new = seq_len(length(a[[h]])) + idx.last
		i[idx.new] = h
		j[idx.new] = a[[h]]
		idx.last = idx.last + length(a[[h]])
	}
	sparseMatrix(i = i, j = j, x = 1, dims = c(L, L))
}

#' @export
intersection_matrix = function(dom1, dom2) {
	D = dom1 %>%
		mutate(D_row_num = row_number()) %>%
		st_set_agr("constant")
	G = dom2 %>%
		mutate(G_row_num = row_number()) %>%
		st_set_agr("constant")
	INT = st_intersection(D, G)
	sparseMatrix(i = INT$G_row_num, j = INT$D_row_num,
		x = as.numeric(st_area(INT)), dims = c(nrow(G), nrow(D)))
}

#' @export
covariance_approximant = function(Sigma, S) {
	eig = eigen(Sigma)
	P = Re(eig$vectors)
	D = Re(eig$values)
	D[D < 0] = 0
	Dinv = D
	Dinv[D > 0] = 1 / D[D > 0]
	Kinv = P %*% (Dinv * t(P))
	return(Kinv)
}

#' @export
precision_approximant = function(Sigma, S) {
	
}

