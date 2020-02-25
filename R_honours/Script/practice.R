# Gabby Riefesel, February 25, 2020

library(adehabitatHR)
# package provdes classes and methods for dealing with home range analysis

#----------------Basic Summary of the Functions of the Package------------------#
xy<-matrix(runif(60), ncol=2)
head(xy)
#sample of relocations dist. on plane

xysp<-SpatialPoints(xy)
#convert the locations into the class spatial points

#fnc clusthr, implements single-linkage custering algorithm, using the fnc clusthr on this object of class 'spatialpolygonsdataframe'
clu<-clusthr(xysp)
class(clu)

#show results
plot(clu)

#consider relocataions of another animal
xy2<-matrix(runif(60), ncol=2)
#bind matrices together
xyt<-rbind(xy, xy2)
#generate a vector containing the ID of the animals for each relocation of the object xyt
id<-gl(2,30)
#convert object id into spatialpointsdataframe
idsp<-data.frame(id)
coordinates(idsp)<-xyt
class(idsp)

# use same fnc clusthr on new object
clu2<-clusthr(idsp)
class(clu2)
clu2
#********** Multiple convex hull Home range of several Animals ************

#This object is a list with one component per animal.
#Each component is an object of class SpatialPolygonsDataFrame
#The home range has been estimated for the following animals:
# [1] "1" "2"

#an object of class MCHu is a list of objects of class 'spatialpolygonsdataframe, with one element per animal

length(clu2)
class(clu2[[1]])
class(clu2[[2]])

#the same HR estimating function can accept an object of class spatialpoints(one animal) or spatialpointsdataframe(several animals)
plot(clu2)

#-------------Example Dataset-------------------------#
#wildboar ex dataset, load the dataset puechabonsp
data("puechabonsp")
names(puechabonsp)

# the dataset is a list of 2 components,
# relocs is a spatialpointsdataframe containing the relocations of 4 wild boar
# following variables are available for each relocation

head(as.data.frame(puechabonsp$relocs))
#  Name Age Sex   Date      X       Y
#1 Brock   2   1 930701 699889 3161559
#2 Brock   2   1 930703 700046 3161541
#3 Brock   2   1 930706 698840 3161033
#4 Brock   2   1 930707 699809 3161496
#5 Brock   2   1 930708 698627 3160941
#6 Brock   2   1 930709 698719 3160989
# the first variable is the idenitity of the animals, use this for estimation process

#the component map is a spatialpixelsdataframe containing the maps of 4 environmental variables on the study area(elevation, slop, aspect, herbaceous cover
#Look at distribution of the distribution of the relocations on an elevation map

##map of elevation
image(puechabonsp$map, col=grey(c(1:10/10)))
##map of relocations
plot(puechabonsp$relocs, add=TRUE,col=as.data.frame(puechabonsp$relocs)[,1])

#------------------Kernel Estimation and the Utilization distribution---------#
#A) Utilization Distribution (UD)

#UD model, we consider that the animals use of space can be described as a bivariate probability density function, the UD, which gives the probability density to relocate teh animal at any place according to the coordinates xy of this place

#Input = relocations
#Output = the UD
#Methods: a bivariate kernel function is placed over each relocation, and the values of these functions are avergaed together
#bivariate kernel method default in this code
#variables
  # h is smoothing parameter,  can be computed by LSCV, the estimated value then minimizes the mean integrated square error (MISE) - teh difference in the vol bw the true UD and estimated UD, also the subjective visual choice for h based on sucessive trials is often sensible choice
  # n is number of relocations
  # Xi is the ith relocation of the sample

######## The fnc kernelUD: estimating the UD
# the fnc kernelUD, implements the methods ^^ to estimate a UD, the UD is estimated in each pixel of a grid superposed to the relocations. 
# the functions can either take an argument of spatialpoints object (1 animal) or spatialpointdatafram object (multiple animals)
# if h = href the reference bandwidth is used in the estimation, or alternatively it is possible to pass a numeric value for h

data("puechabonsp")
kud<-kernelUD(puechabonsp$relocs[,1], h="href")
kud
#********** Utilization distribution of several Animals ************
  
#  Type: probability density
#Smoothing parameter estimated with a  href smoothing parameter
#This object is a list with one component per animal.
#Each component is an object of class estUD
#See estUD-class for more information

#the resulting object is a list of class estUD this class extends the class spatialpixelsdataframe, now containing an additional attribute storing the info about h
#display results
image(kud)

# to get the h-value for the first animal
kud[[1]]@h

#-----------The least square cross validation LSCV-----------------#
#alt youu could set the h equal to LSCV, this algorith, searched for optimum h values in the interval specified by the parameter hlim
#ex est UD with a h parameter chosen with the LSCV algorithm for the puechbonsp dataste

kudl<-kernelUD(puechabonsp$relocs[,1], h="LSCV")
image(kudl)

#^^^^ the resulting UD fits the data more closely than href
## look at the results of the lscv minimization using the fnc plotLSCV, as the cv criterion cannot be minimized in some cases

plotLSCV(kudl)

#^^ here the minimization has been achieved in the specified interval, when the algoruthm does not converge toward a solution the estimate should not be used in further analysis
# you can also choose a numeric value for h, https://www.faunalia.it/dokuwiki/doku.php?id=public:animove

#----------------Controlling the Grid-----------------------------#
#The UD is estimated at the center of each pixel of a grid. To control the paaramters defining the grid , use the parmater grid functiion kernelUD

#-------------Passing a numeric value------------------------------#
#when a numeric value is passed to the parameter zgrid, the function kernelUD auto. computes the grid for the estimateion. the grid is computed in the following wat
  #the min and max coor. Xmin, Xmax, Ymin, Ymax of the relocations calc
  #the range of X and Y values is computed diff bw min and max value. Let Rx and Ry be these ranges
  #the coverage of the grid is defined thanks to the parameter 'extent' of the fnc 'kernelUD'. let 'e' be this parameter. The minimum and max Y coor. of the pixels of the grid are equal to the Ymin-(e x Ry) and Ymax+ (e X Ry) respectively
  #these extends are split into 'g' intervals where g is the value of the paramater grid
  #the paramter grid controls the resolution of the grid and the parameter extent contolrs its extent. 

## Example
#The relocations of "Brock"
locs<-puechabonsp$relocs
firs<-locs[as.data.frame(locs)[,1]=="Brock",]
#graphical parameters
par(mar=c(0,0,2,0))
par(mfrow=c(2,2))
#estimation of the UD with grid=20 and extent=0.2
image(kernelUD(firs, grid=20, extent=0.2))
title(main="grid=20, extent=o.2")
#estimation of the UD grid with grid=80 and extent=0.2
image(kernelUD(firs, grid=80, extent=0.2))
title(main="grid=20, extent=o.2")
#estimation of the UD grid with grid=20 and extent=3
image(kernelUD(firs, grid=20, extent=3))
title(main="grid=20, extent=o.2")
#estimation of the UD grid with grid=80 and extent=3
image(kernelUD(firs, grid=80, extent=3))
title(main="grid=20, extent=o.2")

#NOte that the paramter 'same4all' is an additional parameter allowing the control of the grid. If it is equal to TRUE, teh same grid is used for all animals....
# Example
kus<- kernelUD(puechabonsp$relocs[,1], same4all=TRUE)
image(kus)
#bc all the UD are est on the same grid, is its possible yo coerce the resulting objct as a 'spatialpixeldataframe'
ii<-estUDm2spixdf(kus)
class(ii)
#this object can then be handed with the functions of the package 'sp'


#----------------------Passing a SpatialPixelsDataFrame------------------#
#It can be useful to est the UD in each pixel of a raster map aready available. Esp. when the value of the UD is to be later related to the value of environmental varibales (habitat selection), it is possible to pass a raster map to the parameter grid of the function, ex) the data set puechabonsp contains an object of class spatialpixelsdataframe storing the environmental info for this varibale. This map can be passed to the parameter grid

kudm<-kernelUD(puechabonsp$relocs[,1], grid=puechabonsp$map)

# this object can be coerced into a spatialpixeldataframe using fnc estUDm2spixdf or to pass a list of spatialpixelsdataframe (one element/animal) to the paramter 'grid'

#---------------------Estimated the home range from the UD--------#
# the UD gives the prob density to relocate the animal at any given place. the home range deduced from the UD as the min area on which the prob. to relocate the animal is equal to a specified value Ex tehe 95% HR corresponds w the smallest area on which the probability to relocate the animal is equal to 0.95. The fnc 'getvolumeUD and getverticeshr provide utilities for hr estimation

#### HR in vector mode
#getverticeshr. to deduce the 90% HR from the UD estimated using the LSCV algorthm
homerange<-getverticeshr(kudl)
class(homerange)

# the resulting object ^ is of the class spatialpolygonsdataframe, as for teh MCP the functions of the package sp and maptools are avail to deal with this object or to export it toward a GIS
#to display hr
plot(homerange, COL=1:4)

#### hr in raster mode
#getvolumeUD useful to sestimate hr in raster mode from the UD, this fnc modifies the UD component of the object passed as argument so that the value of a pixel is equal to the % of the smallest hr containing this pixel
vud<-getvolumeUD(kudl)
vud
#********** Utilization distribution of several Animals ************
  
#  Type: volume under UD
#Smoothing parameter estimated with a  LSCV smoothing parameter
#This object is a list with one component per animal.
#Each component is an object of class estUD
#See estUD-class for more information

#to make clear difference bw the output of kernelUD and getvolumeUD look at teh values on the following counterplot....
  #set up graphical parameters
par(mfrow=c(2,1))
par(mar=c(0,0,2,0))
  #the output pg kernelUD for the 1st animal
image(kudl[[1]])
title("Output of KernelUD")
  #convert into suitable data structure for use of contour
xyz<-as.image.SpatialGridDataFrame(kudl[[1]])
contour(xyz, add=TRUE)
  # repeat for output of getvolumeUD
par(mar=c(0,0,2,0))
image(vud[[1]])
title("Output of getvolumeUD")
xyzv<-as.image.SpatialGridDataFrame(vud[[1]])
contour(xyzv, add=TRUE)

#the output kernelUD is the raw UD, the ouput of getvolumeUD can be used to compute the hr 
#store the volume under the UD (as computed by getvolumeUD)
# of the 1st animal in fud
fud<-vud[[1]]
#store the value of the volume under the UD in a vector hr95
hr95<- as.data.frame(fud)[,1]
#if hr95 is <=95 then the pxel belongs to the hr
#take the value 1,0 otherwise
hr95<- as.numeric(hr95<=95)
#convert into dataframe
hr95<- data.frame(hr95)
#covert to a spatialpixelsdataframe
coordinates(hr95)<-coordinates(vud[[1]])
gridded(hr95)<-TRUE
#display results
image (hr95)


#---------The hr size--------#
#kernel.area takes the UD as an argument. computes the hr size for several probabilities

ii<-kernel.area(kudl, percent=seq(50, 95, by=5))
ii
# resulting object is of the class hrsize. it can be plotted using th efnc plot
# the hr sizes retuened by this fnc are dslightly diff from the hr size stored in the spatial polygonsdataframe returned by the fnc getverticeshr. The former measured the area covered bu the raster hr(area covered by the set of pixels of the grid included in the hr with smoother contour), the difference bw the 2 estimates decrease as the resolution of the grid becomes finer

