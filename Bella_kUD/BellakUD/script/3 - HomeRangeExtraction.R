# Author: Isabella Richmond (code and data shared between Matteo Rizzuto: github.com/matteorizzuto)
# Last Edited: June 10, 2020

# This script is for estimating the predation risk of snowshoe hare habitat using structural
# complexity of the environment

# KernelEstimation.R shows how the kUD and home range areas were calculated
# HabitatComplexity.R shows the calculation of range use ratio and ordination 
# of habitat complexity/predation risk values

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
# load vUD raster layer 
vUDBrick <- brick("output/vUDRaster.grd")
# load home range areas and polygons 
# ratio in rangeuse refers to 50%:95% home range area (ha)
rangeuse <- read_csv("output/rangeuseratio.csv")
kernel95 <- read_sf("output/hares.kudhr.95.shp")
kernel95 <- st_transform(kernel95, crs = st_crs(vaancn))
kernel50 <- read_sf("output/hares.kudhr.50.shp")
kernel50 <- st_transform(kernel50, crs = st_crs(vaancn))

# visualize home ranges 
ggplot(kernel95) +
  geom_sf(fill = "dark grey") +
  geom_sf(data = kernel50, aes(fill = "light grey"))+
  xlab("Longitude") + ylab("Latitude")

# convert datasets into tibbles and code dates to make manipulation easier 
predrisk <- as_tibble(predrisk, .name_repair = "universal") %>%
  mutate(Date = lubridate::ymd(Date))

rangeuse <- as_tibble(rangeuse)

# make predation risk data spatially explicit by associating coordinate data 
predriskspatial <- inner_join(predrisk, bl_cs_pts, by = "Plot")
# convert VAAN rasters to dataframe to plot in ggplot
vaancnclip.df <- as.data.frame(vaancnclip, xy=TRUE)
vaancpclip.df <- as.data.frame(vaancpclip, xy=TRUE)
# plot the stoich layer with the sampling points
cn95cs <- tm_shape(vaancnclip)+
  tm_raster(title = "VAAN C:N", style = "cont", 
            palette = "-RdYlBu")+
  tm_scale_bar()+
  tm_layout(legend.bg.color = "white")+
tm_grid()+
tm_shape(kernel95)+
  tm_polygons(alpha = 0.2)+
tm_shape(predriskspatial$geometry)+
  tm_dots(size = 0.15)
tmap_save(cn95cs, "graphics/CN_95_CS.png")
# convert predriskspatial into a SpatialPointsDataFrame
pts <- cbind(predriskspatial$POINT_X_x, predriskspatial$POINT_Y_y)
predriskspatial <- SpatialPointsDataFrame(pts, predriskspatial, proj4string = CRS("+proj=tmerc +lat_0=0 +lon_0=-61.5 +k=0.9999 +x_0=304800 +y_0=0 +ellps=GRS80 +units=m +no_defs"))




# extract the stoich and kUD values at each sampling point
# create raster brick of kUD values and stoich values for easier extraction 
stoichkUD <- brick(list(vaancnclip, vUDBrick))
csValues <- raster::extract(stoichkUD, cchc.spatial, df = TRUE)
# want to combine cchc and csValues dataframes so that there are all values for 
# each sampling plot 
cchc <- as_tibble(cchc) %>%
  mutate(Plot = as.character(Plot))
csValues <- as_tibble(csValues) %>%
  mutate(ID = as.character(ID))
# use bind and not join because the plots are in the same order 
cscchc <- bind_cols(csValues, cchc)
# get median value of all collars for each sampling plot 
cscchc <- as_tibble(cscchc) %>%
  rowwise() %>% 
  mutate(kUDmedian =  median(c(!!! rlang::syms(grep('X', names(.), value=TRUE)))))
