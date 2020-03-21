# A script to estimate home range size of the radio collared snowshoe hares in 
# the Bloomfield Grid.
# 
# Should probably state in the manuscript (whenever it will happen) that we were 
# only interested in summer Home Ranges, hence we only did telemetry in 
# July-August and collected ~30 relocations per animal.
# This should give a good and conservative estimate of the home range for the 
# summer season, but likely won't tell us much about their winter home range
# So our inference is limited

#------------------------------------#
#                 Setup              #
#------------------------------------#

# source the Data Cleaning script
#source("Code/NEWCollarLocsEstimation.R")

#------------------------------------#
#       Home Range Estimation        #  
#------------------------------------#

# I will estimate the home range of these snowshoe hares using two widespread 
# methods: the Minimum Convex Polygon (MCP) and the Kernel Utilization 
# Distribution (KUD).

# Because of how the package adeHabitatHR works, it is much easier to work on 
# a single dataframe containing the data from each animal. The first step is, 
# then, to unify the 4 separate dataframes containing the results from the 
# sigloc() triangulation. Then, I will move on to estimate the HR using both MCP
# and KUD.

# 1. Bind the 4 datasets

# Despite multiple internet help pages stating the opposite, rbind() works just 
# fine to bind multiple SpatialPointsDataFrame objects together into a single, 
# massive one.
hares.triangd <- do.call(rbind, hares.list)

# Kernel UD
# Now let's move on to the Kernel Estimation. The following estimate uses Lest 
# Square Cross Validation (LSCV) to estimate the smoothing parameter

hares.kUD <- kernelUD(hares.triangd[,8], h = 'href', grid = 1000, extent = 1,
                      same4all = FALSE)

# let's visualize it
image(hares.kUD)

# to export it, let's first convert it to a SpatialPixelDataFrame
hareskUD_export <- estUDm2spixdf(hares.kUD)

# save the subset of data of interest as a raster
c452raster <- raster(hares.kUDexport[2])

# save it on disk
writeRaster(c452raster, "../Results/c452raster.tif", "GTiff")

# Now let's look at the home range size using a range of probability levels
kUD.hr.estimates <- kernel.area(hares.kUD, percent = seq(50, 100, 5), 
                                unin = "m", unout = "ha")

# now, let's extract the home range
hares.kUDhr <- getverticeshr(hares.kUD, percent = 90)

# and now save each home range as a separate shapefile for ue in GIS

c535.kUDhr <- hares.kUDhr[hares.kUDhr$id == "149.535",]
c673.kUDhr <- hares.kUDhr[hares.kUDhr$id == "149.673",]
c452.kUDhr <- hares.kUDhr[hares.kUDhr$id == "149.452",]
c394.kUDhr <- hares.kUDhr[hares.kUDhr$id == "149.394",]
c003.kUDhr <- hares.kUDhr[hares.kUDhr$id == "149.003",]
c053.kUDhr <- hares.kUDhr[hares.kUDhr$id == "149.053",]
c093.kUDhr <- hares.kUDhr[hares.kUDhr$id == "149.093",]
c173.kUDhr <- hares.kUDhr[hares.kUDhr$id == "149.173",]
c213.kUDhr <- hares.kUDhr[hares.kUDhr$id == "149.213",]
c274.kUDhr <- hares.kUDhr[hares.kUDhr$id == "149.274",]
c374.kUDhr <- hares.kUDhr[hares.kUDhr$id == "149.374",]
c474.kUDhr <- hares.kUDhr[hares.kUDhr$id == "149.474",]
c633.kUDhr <- hares.kUDhr[hares.kUDhr$id == "149.633",]
c653.kUDhr <- hares.kUDhr[hares.kUDhr$id == "149.653",]

shapefile(c535.kUDhr, "../Results/Coolar535_KUDhr.shp", overwrite = TRUE)
shapefile(c673.kUDhr, "../Results/Coolar673_KUDhr.shp", overwrite = TRUE)
shapefile(c452.kUDhr, "../Results/Coolar452_KUDhr.shp", overwrite = TRUE)
shapefile(c394.kUDhr, "../Results/Coolar394_KUDhr.shp", overwrite = TRUE)
shapefile(c003.kUDhr, "../Results/Coolar003_KUDhr.shp", overwrite = TRUE)
shapefile(c053.kUDhr, "../Results/Coolar093_KUDhr.shp", overwrite = TRUE)
shapefile(c093.kUDhr, "../Results/Coolar093_KUDhr.shp", overwrite = TRUE)
shapefile(c173.kUDhr, "../Results/Coolar173_KUDhr.shp", overwrite = TRUE)
shapefile(c213.kUDhr, "../Results/Coolar213_KUDhr.shp", overwrite = TRUE)
shapefile(c274.kUDhr, "../Results/Coolar274_KUDhr.shp", overwrite = TRUE)
shapefile(c374.kUDhr, "../Results/Coolar374_KUDhr.shp", overwrite = TRUE)
shapefile(c474.kUDhr, "../Results/Coolar474_KUDhr.shp", overwrite = TRUE)
shapefile(c633.kUDhr, "../Results/Coolar633_KUDhr.shp", overwrite = TRUE)
shapefile(c653.kUDhr, "../Results/Coolar653_KUDhr.shp", overwrite = TRUE)

# Now we can transform the kernel-estimated home range into a dataframe for 
# plotting with ggplot2
hares.kUDhr <- fortify(hares.kUDhr)

hares.kUDhr$id <- as.factor(hares.kUDhr$id)
hares.kUDhr.BL <- hares.kUDhr[! (hares.kUDhr$id == "149.374" | hares.kUDhr$id == "149.474"),]
hares.kUDhr.TN <- hares.kUDhr[(hares.kUDhr$id == "149.374"),]
hares.kUDhr.UNI <- hares.kUDhr[(hares.kUDhr$id == "149.474"),]

# finally, let's plot it as above using ggplot2
ggplot(hares.kUDhr.BL, aes(x = long, y = lat, col = id, 
                        fill = id, group = group)) +
  geom_polygon(alpha = 0.25) +
  # scale_fill_manual(values  = cbPalette, name = "Collar Frequency") +
  # scale_colour_manual(values  = cbPalette2, name = "Collar Frequency") +
  coord_equal() +
  theme_map()

ggplot(hares.kUDhr.TN, aes(x = long, y = lat)) +
  geom_polygon(alpha = 0.25) +
  coord_equal() +
  theme_map()

ggplot(hares.kUDhr.UNI, aes(x = long, y = lat)) +
  geom_polygon(alpha = 0.25) +
  coord_equal() +
  theme_map()

# Another option to visualize the Kernel UD is via raster. This methods 
# estimates the volume of the home range, and can give slightly different 
# results from the vector/area method used above. Let's see how it works out.

hares.vol <- getvolumeUD(hares.kUD)

# and let's visualize it 
image(hares.vol)

# To get the rasterized 95% hr of each animal, the procedure is a bit convoluted
# and needs to be done one individual at the time. Here is an example:

fud673 <- hares.vol[[4]] # separate the 1st individual's data

vol673hr95 <- as.data.frame(fud673)[,1]  # save it as a df
 
vol673hr95 <- as.numeric(vol673hr95 <= 95) # select points below 95%

vol673hr95 <- data.frame(vol673hr95) # save as df again

coordinates(vol673hr95) <- coordinates(fud673[1]) 
# covert to SpatialPixelsDataFrame

gridded(vol673hr95) <- TRUE # grid it

image(vol673hr95)

writeGDAL(vol673hr95, "../Results/c673volHR.tif", "GTiff") # save it as raster
