# Author: Isabella Richmond (code and data shared between Matteo Rizzuto: github.com/matteorizzuto)
# Last Edited: June 11, 2020

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

# --------------------------------------- #
#               Extract Data              #
# --------------------------------------- #
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

# extract the stoich values for each home range
# create raster brick of kUD values and stoich values for easier extraction 
stoich <- brick(list(vaancnclip, vaancpclip))
# extract the mean and the individual measurements for each home range 
stoichhr <- raster::extract(stoich, kernel95, df = TRUE,na.rm =  TRUE)
stoichhrmean <- stoichhr %>% drop_na() %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise(MeanCN = mean(VAAN_CN), MeanCP = mean(VAAN_CP), NumberCells = n())
# extract the predation risk values for each home range 
# need pred risk to be SpatialPoints and kernel95 to be SpatialPolygons
predrisksf <- st_as_sf(predriskspatial)
predriskhr <- sf::st_intersection(kernel95, predrisksf)
# want to calculate mean for each home range and also calculate the number of measurements per home range
predriskhr <- dplyr::rename(predriskhr, CollarFrequency = id)
predriskhr <- predriskhr %>% dplyr::group_by(CollarFrequency) %>%
  dplyr::summarise(MeanPC1 = mean(PC1), MeanPC2 = mean(PC2), NumberPoints = n_distinct(Plot))
# join predation risk and stoich data together 
# rasters extract in the same order as the polygon ID - sort predation risk by collar and then join
predriskhr <- dplyr::arrange(predriskhr, CollarFrequency)
prstoichhr <- bind_cols(predriskhr, stoichhrmean)
# join home range area and ratio data 
rangeuse <- dplyr::arrange(rangeuse, CollarFrequency)
finaldata <- bind_cols(prstoichhr, rangeuse)
write_csv(finaldata, "output/RangeStoichRisk.csv")

# now extract the stoich and predation risk values for each core area 
