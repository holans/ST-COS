library(maptools)
library(raster)
library(shapefiles)

res$county@data[1,]
proj4string(res$county)
coordinates(res$county, proj4string = "")

coordinates(res$county[1,])

coordinates(spTransform(res$county[1,], CRS("+proj=longlat +datum=WGS84")))


county <- readOGR("shp", layer = "period1", verbose = FALSE)
head(county@data)

project2<-"+proj=eqdc +lat_0=0 +lon_0=0 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
data.shape <- readShapePoly("shp/period1", proj4string=CRS(project2), verbose=TRUE)
head(data.shape@data)
str(data.shape@data)

shp <- shapefile("shp/period1", verbose = TRUE, stringsAsFactors=TRUE)
shp

