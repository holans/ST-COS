require(rgdal)
require(MCMCpack)

# Read the counties in Connecticut
path <- "../data/tl_2010_09_county10/tl_2010_09_county10.shp"
counties <- readOGR(path, layer = "tl_2010_09_county10")
plot(counties)

# Select a county and plot it
idx <- which(counties$NAME10 == 'Fairfield')
plot(counties[idx,])

# Sample uniformly from the county
# A bounding box around the county comes in the shapefile
# Sample from that uniformly and subset to the county

lon.low <- a@bbox[1,1]
lon.hi <- a@bbox[1,2]
lat.low <- a@bbox[2,1]
lat.hi <- a@bbox[2,2]

n <- 4000
u.lon <- runif(n, lon.low, lon.hi)
u.lat <- runif(n, lat.low, lat.hi)
u <- data.frame(Longitude = u.lon, Latitude = u.lat)
coordinates(u) <- ~ Longitude + Latitude
proj4string(u) <- proj4string(counties)
points(u)

# Now reject the ones outside of the county...
idx <- which(over(u, counties)$NAME10 == "Fairfield")
points(u[idx], col = "red")
