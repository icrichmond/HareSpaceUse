# Author: Isabella Richmond
# Last Edited: July 15, 2020

# This script is for extracting the food quality and predation risk of each individual's
# core area and home range.
# aKDERazimuth.R shows how the kUD and home range areas were calculated
# RiskOrdination.R shows the ordination of habitat complexity/predation risk values

# load required packages 
easypackages::packages("ggpubr", "patchwork", "AICcmodavg", "broom", 
                       "maptools", "sf", "matrixStats", "tidyverse", "lubridate", 
                       "dplyr", "ggplot2", "ggcorrplot", "sf", "raster", "lme4", 
                       "RCurl", "tibble", "rgdal","tmap", "spatialEco")

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
bl_cs_pts <- st_transform(bl_cs_pts, crs ="+init=epsg:32622")
plot(bl_cs_pts, add=T)
# load home range area rasters - in a list (95% home range area in hectares)
kernel95 <- readRDS("large/akderasters.rds")
# normalize rasters so they are between 0,1
# also change 0 to No Data so that mapping is easier
scale01 <- function(r) {
  rmin <- cellStats(kernel95[[1]], 'min')
  rmax <- cellStats(kernel95[[1]], 'max')
  (r - rmin) / (rmax - rmin)
}

kernel95norm <- lapply(kernel95, scale01)

image(kernel95norm[[1]])
image(kernel95[[1]])

# --------------------------------------- #
#               Extract Data              #
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
cskud <- extract(kernel95, predriskspatial)
cskud <- as.data.frame(cskud)

