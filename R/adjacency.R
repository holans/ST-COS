# Convert an adjacency list into a (sparse) matrix.
#' @export
adjacency_matrix = function(dom)
{
	adj = st_touches(dom, dom)
	n_tuples = sum(unlist(Map(length, adj)))
	i = integer(n_tuples)
	j = integer(n_tuples)
	x = rep(1, n_tuples)
	idx.last = 0
	L = length(adj)
	for (h in 1:L) {
		idx.new = seq_len(length(adj[[h]])) + idx.last
		i[idx.new] = h
		j[idx.new] = adj[[h]]
		idx.last = idx.last + length(adj[[h]])
	}
	sparseMatrix(i = i, j = j, x = 1, dims = c(L, L))
}
