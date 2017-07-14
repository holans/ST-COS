rArea <- function(n, area, blocksize = n)
{
	N <- nrow(area)
	res <- matrix(NA, 0, 2)

	done <- FALSE
	while (!done) {
		u <- st_sample(area, size = blocksize)
		if (NROW(u) > 0) {
			M <- matrix(unlist(u), length(u), 2, byrow = TRUE)
			res <- rbind(res, M)
			done <- (nrow(res) >= n)			
		}
	}

	colnames(res) <- c("x", "y")
	matrix(res[1:n,], ncol = 2)
}

