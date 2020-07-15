# Author: Isabella Richmond (code and data shared between Matteo Rizzuto: github.com/matteorizzuto)
# Last Edited: June 17, 2020

# This script is for extracting the food quality and predation risk of each individual's
# core area and home range.
# aKDERazimuth.R shows how the kUD and home range areas were calculated
# RiskOrdination.R shows the ordination of habitat complexity/predation risk values

# load required packages 
easypackages::packages("ggpubr", "patchwork", "AICcmodavg", "broom", 
                       "maptools", "sf", "matrixStats", "tidyr", "lubridate", 
                       "dplyr", "ggplot2", "ggcorrplot", "sf", "raster", "lme4", 
                       "RCurl", "tibble", "rgdal")

# --------------------------------------- #
#             Data Preparation            #
# --------------------------------------- #
# import predation risk dataset with PCA axes data
predrisk <- read_csv("output/predationriskpca.csv")
# load lowland blueberry C:N and C:P stoich rasters (as per JBF's and MR's results)
vaancn <- raster("input/VAAN_CN.tif")
vaancp <- raster("input/VAAN_CP.tif")
# clip raster to study area 
e <- extent(860000, 863000, 5383000, 5386000)
vaancnclip <- crop(vaancn, e)
image(vaancnclip)
vaancpclip <- crop(vaancp, e)
# load complexity sampling locations shapefile
bl_cs_pts <- read_sf("input/Mapping", layer = "cs_points")
bl_cs_pts <- st_transform(bl_cs_pts, crs = st_crs(vaancn))
# load home range area raster brick (95% home range area in hectares)
# ratio in rangeuse refers to 50%:95% home range area (ha)
kernel95 <- raster::stack("output/Rasters/akde_homerange_pdf.tif")
# reproject kernels to match other CRS
crs(kernel95) <- crs(vaancn)

tmap::tm_shape(kernel95)+
  tmap::tm_raster(max.value=1)


