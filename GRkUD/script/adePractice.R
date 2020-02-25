# This code is running through the examples in the documentation for the adehabitatHR package
# Project is part of a git repository for GR's honours (kernels)

# Load package
library(adehabitatHR)
# Creating relocations of an animal
xy <- matrix(runif(60), ncol=2)
head(xy)
# Convert xy relocations into a class of SpatialPoints
xysp <- SpatialPoints(xy)
# clusthr function implements single-linkage clustering algorithm - returns SpatialPolygonsDataFrame
clu <- clusthr(xysp)
class(clu)
plot(clu)
# Create relocations of a second animal 
xy2 <- matrix(runif(60), ncol=2)
# Bind the two animal matrices together 
xyt <- rbind(xy, xy2)
# Generate vector containing ID of animal for each relocation of xyt
id <- gl(2,30)
# Convert object id into SpatialPointsDataFrame 
idsp <- data.frame(id)
coordinates(idsp) <- xyt
class(idsp)
# Use clusthr on new object 
clu2 <- clusthr(idsp)
class(clu2)
clu2
# This creates a home range estimation for each animal with the objects being SpatialPolygonsDataFrame
# Class "MCHu" is a list of objects of SpatialPolygonsDataFrame with one element per animal
length(clu2)
class(clu2[[1]])
class(clu2[[2]])
plot(clu2)
# Load example dataset
data("puechabonsp")
names(puechabonsp)
head(as.data.frame(puechabonsp$relocs))
# Map of the elevation 
image(puechabonsp$map, col=grey(c(1:10)/10))
# Map of the relocations 
plot(puechabonsp$relocs, add=TRUE, col=as.data.frame(puechabonsp$relocs) [,1])


### Moving on to kernels ###

