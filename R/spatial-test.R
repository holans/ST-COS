require(rgdal)
require(MCMCpack)

# Read the counties in Alabama
path <- "~/Documents/datasets/alabama2010_countiestl_2010_01_county10.shp"
# path <- "H:/work/projects/maf/20150225-revisions/mapping/counties/tl_2010_01_county10.shp"
counties <- readOGR(path, layer = "tl_2010_01_county10")
plot(counties)

# Select Autauga county and plot it
idx <- which(counties$NAME10 == 'Autauga')
plot(counties[idx,])

# How to sample uniformly from Autauga county?
#   It's not a convex shape, so we can't just find its convex hull.
#   We could put a convex envelope around it and sample from that!

chull(x, y = NULL)

a <- counties[idx,]
a@polygons

V <- matrix(NA, 4, 2)
V[1:2,] <- t(a@bbox)
V[3,] <- c(V[1,1], V[2,2])
V[4,] <- c(V[2,1], V[1,2])

points(V)

lambda <- rdirichlet(1000, rep(1,4))
u <- lambda %*% V
points(u, col = "blue")

u <- as.data.frame(u)
colnames(u) <- c("Longitude", "Latitude")
coordinates(u) <- ~ Longitude + Latitude
proj4string(u) <- proj4string(counties)

# Now reject the ones outside of the county...
idx <- which(over(u, counties)$NAME10 == "Autauga")
points(u[idx], col = "red")

