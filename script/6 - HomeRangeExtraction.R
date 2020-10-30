# Author: Isabella Richmond
# Last Edited: October 30, 2020

# This script is for extracting the food quality and kUD values at each habitat complexity
# sampling point.
# 3 - aKDERazimuth.R shows how the kUD and home range areas were calculated
# 5 - RiskOrdination.R shows the ordination of habitat complexity/predation risk values

# load required packages 
easypackages::packages("matrixStats", "tidyverse", "lubridate", "sf", "sp", "raster",
                       "tmap", "spatialEco", "grid", "osmdata")

# --------------------------------------- #
#             Data Preparation            #
# --------------------------------------- #
# import predation risk dataset with PCA axes data
# note - had to drop some of the complexity sampling points because 
# they had NA values (doesn't work in a PCA), there are 67 instead of 72
predrisk <- read_csv("output/predationriskpca.csv")
# load lowland blueberry C:N and C:P stoich rasters (as per JBF's and MR's results)
vaancn <- raster("input/VAAN_CN.tif")
vaancp <- raster("input/VAAN_CP.tif")
# clip raster to study area 
e <- extent(860000, 863000, 5383000, 5386000)
vaancnclip <- crop(vaancn, e)
image(vaancnclip)
vaancpclip <- crop(vaancp, e)
# reproject rasters so that they are using an up to date datum (WGS84 UTM Zone 22)
vaancnclip <- projectRaster(vaancnclip, crs="+init=epsg:32622")
vaancpclip <- projectRaster(vaancpclip, crs="+init=epsg:32622")
image(vaancnclip)
# load complexity sampling locations shapefile
bl_cs_pts <- read_sf("input/Mapping", layer = "cs_points")
bl_cs_pts <- sf::st_transform(bl_cs_pts, "+init=epsg:32622")
plot(bl_cs_pts, add=T)
# load home range area rasters - in a list (95% home range area in hectares)
kernel95 <- readRDS("large/akderasters.rds")
# normalize rasters so they are between 0,1
kernel95norm <- lapply(kernel95, spatialEco::raster.transformation)
# check that they are the same spatially
image(kernel95norm[[1]])
image(kernel95[[1]])
# create raster stack so extraction is possible
kernel95stack <- stack(kernel95norm)
names(kernel95stack) <- names(kernel95)

# --------------------------------------- #
#       Extract Data - Individuals        #
# --------------------------------------- #
# convert datasets into tibbles and code dates to make manipulation easier 
predrisk <- as_tibble(predrisk, .name_repair = "universal") %>%
  mutate(Date = lubridate::ymd(Date))
# make predation risk data spatially explicit by associating coordinate data 
predriskspatial <- inner_join(predrisk, bl_cs_pts, by = "Plot")
# convert VAAN rasters to dataframe to plot in ggplot
vaancnclip.df <- as.data.frame(vaancnclip, xy=TRUE)
vaancpclip.df <- as.data.frame(vaancpclip, xy=TRUE)
# plot the stoich layer with the sampling points
cncs <- tm_shape(vaancnclip)+
  tm_raster(title = "VAAN C:N", style = "cont", 
            palette = "-RdYlBu")+
  tm_scale_bar()+
  tm_layout(legend.bg.color = "white")+
tm_grid()+
tm_shape(predriskspatial$geometry)+
  tm_dots(size = 0.15)
cncs
# note how the grid is proportionally high food quality when compared to the surrounding area 
# same pattern for C:P is present 

# convert predriskspatial into a SpatialPointsDataFrame
predriskspatial <- st_as_sf(predriskspatial)
predriskspatial <- as_Spatial(predriskspatial)
# extract the stoich values for each home range
# create raster brick of kUD values and stoich values for easier extraction 
stoich <- brick(list(vaancnclip, vaancpclip))
# extract the stoich values at each complexity sampling point 
csstoich <- extract(stoich, predriskspatial)
csstoich <- as.data.frame(csstoich)
# extract the kUD values at each complexity sampling point 
cskud <- extract(kernel95stack, predriskspatial)
cskud <- as.data.frame(cskud)
# combine CS plot name, stoich values, and kUD values 
stoichkud <- cbind(cskud, csstoich)
stoichkud <- add_column(stoichkud, Plot = predrisk$Plot)
write.csv(stoichkud, "output/cs_stoich_kud.csv")

# use clamp to transform any values under 0.15 to Nas
kernel95normz <- lapply(kernel95norm, function(i) raster::clamp(i, lower=0.15, useValues=FALSE))
# save clamped raster
saveRDS(kernel95normz, "large/rasternormclamp.rds")