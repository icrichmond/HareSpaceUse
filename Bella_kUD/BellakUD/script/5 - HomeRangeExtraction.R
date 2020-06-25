# Author: Isabella Richmond (code and data shared between Matteo Rizzuto: github.com/matteorizzuto)
# Last Edited: June 17, 2020

# This script is for extracting the food quality and predation risk of each individual's
# core area and home range.
# KernelEstimation.R shows how the kUD and home range areas were calculated
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
# load home range areas and polygons 
# ratio in rangeuse refers to 50%:95% home range area (ha)
rangeuse <- read_csv("output/rangeuseratio.csv")
kernel90 <- read_sf("output/hares.kudhr.90.shp")
kernel90 <- st_transform(kernel90, crs = st_crs(vaancn))
kernel50 <- read_sf("output/hares.kudhr.50.shp")
kernel50 <- st_transform(kernel50, crs = st_crs(vaancn))

# --------------------------------------- #
#               Extract Data              #
# --------------------------------------- #
# visualize home ranges 
ggplot(kernel90) +
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
cn90cs <- tm_shape(vaancnclip)+
  tm_raster(title = "VAAN C:N", style = "cont", 
            palette = "-RdYlBu")+
  tm_scale_bar()+
  tm_layout(legend.bg.color = "white")+
tm_grid()+
tm_shape(kernel90)+
  tm_polygons(alpha = 0.2)+
tm_shape(predriskspatial$geometry)+
  tm_dots(size = 0.15)
cn90cs
tmap_save(cn90cs, "graphics/CN_90_CS.png")
# convert predriskspatial into a SpatialPointsDataFrame
pts <- cbind(predriskspatial$POINT_X_x, predriskspatial$POINT_Y_y)
predriskspatial <- SpatialPointsDataFrame(pts, predriskspatial, proj4string = CRS("+proj=tmerc +lat_0=0 +lon_0=-61.5 +k=0.9999 +x_0=304800 +y_0=0 +ellps=GRS80 +units=m +no_defs"))

# extract the stoich values for each home range
# create raster brick of kUD values and stoich values for easier extraction 
stoich <- brick(list(vaancnclip, vaancpclip))
# extract the mean and the individual measurements for each home range 
stoichhr <- raster::extract(stoich, kernel90, df = TRUE,na.rm =  TRUE)
stoichhrmean <- stoichhr %>% drop_na() %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise(MeanCN_Home = mean(VAAN_CN), MeanCP_Home = mean(VAAN_CP), NumberCells_Home = n())
# extract the predation risk values for each home range 
# need pred risk to be SpatialPoints and kernel95 to be SpatialPolygons
predrisksf <- st_as_sf(predriskspatial)
predriskhr <- st_intersection(kernel90, predrisksf)
# want to calculate mean for each home range and also calculate the number of measurements per home range
predriskhr <- dplyr::rename(predriskhr, CollarFrequency = id)
predriskhr <- predriskhr %>% dplyr::group_by(CollarFrequency) %>%
  dplyr::summarise(MeanOverPCA_Home = mean(overPCA), MeanUnderPCA_Home = mean(underPCA), NumberPoints_Home = n_distinct(Plot))
# join predation risk and stoich data together 
# rasters extract in the same order as the polygon ID - sort predation risk by collar and then join
predriskhr <- dplyr::arrange(predriskhr, CollarFrequency)
prstoichhr <- bind_cols(predriskhr, stoichhrmean)
# join home range area and ratio data 
rangeuse <- dplyr::arrange(rangeuse, CollarFrequency)
homerangedata <- bind_cols(prstoichhr, rangeuse)

# now extract the stoich and predation risk values for each core area 
# extract the mean and the individual measurements for each core area 
stoichca <- raster::extract(stoich, kernel50, df = TRUE,na.rm =  TRUE)
stoichcamean <- stoichca %>% drop_na() %>%
  dplyr::group_by(ID) %>%
  dplyr::summarise(MeanCN_Core = mean(VAAN_CN), MeanCP_Core = mean(VAAN_CP), NumberCells_Core = n())
# extract the predation risk values for each core area 
# need pred risk to be SpatialPoints and kernel95 to be SpatialPolygons
predriskca <- st_join(kernel50, predrisksf)
predriskcamean <- predriskca %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(MeanOverPCA_Core = mean(overPCA), MeanUnderPCA_Core = mean(underPCA), NumberPoints_Core = n_distinct(Plot))
predriskcamean <- dplyr::rename(predriskcamean, CollarFrequency=id)
# join predation risk and stoich data together 
# rasters extract in the same order as the polygon ID - sort predation risk by collar and then join
predriskca <- dplyr::arrange(predriskcamean, CollarFrequency)
prstoichca <- bind_cols(predriskcamean, stoichcamean)
# join core area data to home range and ratio data 
prstoichca <- as_tibble(prstoichca)
homerangedata <- as_tibble(homerangedata)
finaldata <- inner_join(prstoichca, homerangedata, by="CollarFrequency")
write_csv(finaldata, "output/RangeStoichRisk.csv")
