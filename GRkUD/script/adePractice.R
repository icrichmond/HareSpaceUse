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

# Trying LSCV as the smoothing parameter as this is generally accepted as a good approach
# Load data
data("puechabonsp")
kud <- kernelUD(puechabonsp$relocs[,1], h="LSCV")
kud
image(kud)
# Values for smoothing parameters are stored in the slot "h" of each element in the list 
kud[[1]]@h
# important to look at LSCV minimization using plotLSCV & make sure minimization has occurred
# within the specified interval
plotLSCV(kud)
# the grid can be controlled using grid and extent parameters 
# same4all controls grid and uses the same grid for all individuals
kus <- kernelUD(puechabonsp$relocs[,1], same4all = TRUE)
image(kus)
# grid controls resolution of the grid and extent controls extent 
# manipulating individuals 
# brock
locs <- puechabonsp$relocs
firs <- locs[as.data.frame(locs)[,1]=="Brock",]
par(mar=c(0,0,2,0))
par(mfrow=c(2,2))
# Estimation of UD with grid = 20 and extent = 0.2
image(kernelUD(firs,grid=20,extent=0.2))
# When same4all is used, can make the object a SpatialPixelsDataFrame because all the UD are 
# on the same grid 
ii <- estUDm2spixdf(kus)
class(ii)
# Can make the grid equal to a different raster (e.g. environmental variables)
kudm <- kernelUD(puechabonsp$relocs[,1], grid = puechabonsp$map)


### Estimating home ranges from UD ###
# Estimate the 90% home range using LSCV algorithm in vector form
homerange <- getverticeshr(kud)
class(homerange)
plot(homerange, col=1:4)
# Can also estimate in raster mode from the UD 
vud <- getvolumeUD(kud)
vud
# Look at differences between the output of kernelUD and getvolumeUD 
par(mfrow=c(2,1))
par(mar=c(0,0,2,0))
# Output of kernelUD for first animal
image(kud[[1]])
title("Output of kernelUD")
# Convert into suitable data structure for the use of contour 
xyz <- as.image.SpatialGridDataFrame(kud[[1]])
contour(xyz, add=TRUE)
# Now for getvolumeUD
par(mar=c(0,0,2,0))
image(vud[[1]])
title("Output of getvolumeUD")
xyzv <- as.image.SpatialGridDataFrame(vud[[1]])
contour(xyzv, add=TRUE)
# store the volume under the UD of the first animal 
fud <- vud[[1]]
# Store the value of volume under UD in a vector hr95
hr95 <- as.data.frame(fud)[,1]
# If hr95 is <= 95 then the pixel belongs to the home range (takes value 1,0 otherwise)
hr95 <- as.numeric(hr95 <= 95)
# Convert into data frame 
hr95 <- data.frame(hr95)
# Convert to SpatialPixelsDataFrame
coordinates(hr95) <- coordinates(vud[[1]])
gridded(hr95) <- TRUE
image(hr95)
# Can also calculate the area of the home range - can do this for several probability levels
ii <- kernel.area(kud, percent = seq(50,95, by=5))
ii
plot(ii)
# Home range sizes returned this way are different than SpatialPolygonsDataFrame that is returned
# using getverticeshr. Different is smaller when grid resolution is finer 
# Can also take into account physical boundaries of the study area - skipping example as BL 
# does not have any physical barriers
