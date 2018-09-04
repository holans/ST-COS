rArea <- function(n, area, blocksize = n)
{
	N <- nrow(area)
	res <- matrix(NA, 0, 2)
	bbox <- st_bbox(area)

	done <- FALSE
	while (!done) {
		## The sf package offers a simple way to draw a uniform sample on an area.
		## However, this currently seems to be slower than the way below.
		# u <- st_sample(area, size = blocksize)
		# if (NROW(u) > 0) {
		#	M <- matrix(unlist(u), length(u), 2, byrow = TRUE)
		#	res <- rbind(res, M)
		#	done <- (nrow(res) >= n)
		# }

		u1 <- runif(blocksize, bbox[1], bbox[3])
		u2 <- runif(blocksize, bbox[2], bbox[4])
		u_cand <- data.frame(u1 = u1, u2 = u2)
		DT_sf <- st_as_sf(u_cand, coords = c("u1", "u2"), crs = st_crs(area), agr = "constant")
		idx <- st_contains(area, DT_sf)[[1]]
		u <- u_cand[idx,]

		if (NROW(u) > 0) {
			res <- rbind(res, as.matrix(u))
			done <- (nrow(res) >= n)			
		}
	}

	colnames(res) <- c("x", "y")
	matrix(res[1:n,], ncol = 2)
}

