# This script is for estimating the space use using kernels for snowshoe hares 
# over three years (2016-2019). Relocations only taken in summer season.

# Author: Isabella Richmond (code and data shared between Matteo Rizzuto: github.com/matteorizzuto)
# Last Edited: June 25, 2020

easypackages::packages("chron", "ctmm", "sp", "sf", "maptools", "tmap", "tmaptools", "SDMTools", 
                       "adehabitatHR", "adehabitatHS", "adehabitatLT", "ellipse", "ggplot2",
                       "nleqslv", "adehabitatMA", "adehabitatHR","dplyr", "gdtools", "ggmap",  
                       "ggrepel", "ggsci", "ggthemes", "maps", "raster", "spatial", "XML", 
                       "tidyr", "readr","rgdal", "rgeos", "reshape2", "dismo", "tibble")

# --------------------------------------- #
#           Data Preparation              #
# --------------------------------------- #
# This code taken directly from Matteo Rizzuto's HomeRangeEstimation.R script
# Find Matteo's repository at github.com/matteorizzuto/Chapter_2

# load triangulated hare data
hares.triangd <- read_csv("input/harestriangd.csv")
# load stoich data that is going to be used - Lowland blueberry (Vaccinium angustifolium) C:N
vaancn <- raster("input/VAAN_CN.tif")
image(vaancn)
# clip raster to study area 
e <- extent(860000, 863000, 5383000, 5386000)
vaancnclip <- crop(vaancn, e)
image(vaancnclip)
# convert clipped raster to SpatialPixelsDataFrame so it can be used in kUD analysis 
vaancnASC <- asc.from.raster(vaancnclip)
vaanCN <- asc2spixdf(vaancnASC)
class(vaanCN)

# --------------------------------------- #
#            Calculate kUDS               #
# --------------------------------------- #

# before estimating the kUD, remove collars 149.555, 150.132 due to too few relocations 
# available for estimating a reliable kernel Utilization Distribution
hares.triangd <- subset(hares.triangd, hares.triangd$Frequency != "149.555" 
                        & hares.triangd$Frequency != "150.132")

# Let's estimate the kernel Utilization Distribution using the ad hoc method and 
# a grid that is set to the same size as the stoich grid, that can adapt to the general 
# geographic area used by each animal
hares.kUD <- kernelUD(hares.triangd[,8], h = 'href', grid = vaanCN, extent = 1 ,same4all = FALSE)
# NOTE for future: stoich resolution and predation risk resolution (30mx30m) is used 
# for kUD calculations so that all layers are the same resolution

# If reverting back to using LSCV to estimate h, double-check that minimization
# of the cross-validation criteria is successful using:
#par(mar = c(1,1,1,1), mfrow = c(1,1))
#plotLSCV(hares.kUD)

# estimate home range in raster form using getvolumeUD
hares.vUD <- getvolumeUD(hares.kUD)
image(hares.vUD)

image(hares.vUD[[1]])
xyzv <- as.image.SpatialGridDataFrame(hares.vUD[[1]])
contour(xyzv, add=TRUE)
# estimate kUD area using a range of percent levels
kUD.hr.estimates <- kernel.area(hares.kUD, percent = seq(20, 95, 5), 
                                unout = "ha")
kUD.hr.estimates
plot(kUD.hr.estimates)

# and extract values to be used in later modelling
hrArea <- kUD.hr.estimates[1:16, ]
hrArea <- rownames_to_column(hrArea)
hrArea <- rename(hrArea, Kernel = rowname)
write_csv(hrArea, "output/homerangeareas.csv")
# calculate range use ratio with 50:95 home range areas 
rangeuse <- pivot_longer(hrArea, cols = 2:ncol(hrArea), names_to = "CollarFrequency", values_to = "HomeRangeArea")
rangeuse <- pivot_wider(rangeuse, names_from = "Kernel", values_from = "HomeRangeArea")
# calculating a 50:95 range use ratio (Webber et al., 2020)
rangeuse <- add_column(rangeuse, ratio = rangeuse$`50`/rangeuse$`95`)
# test if the 50% and 95% home range areas are correlated 
ggplot(data = rangeuse, aes(x = rangeuse$`20`, y = rangeuse$`70`))+geom_point()
# highly correlated, don't use ratio going forward. Use core area as response variable
write_csv(rangeuse, "output/rangeuseratio.csv")


# --------------------------------------- #
#         Extract Home Ranges             #
# --------------------------------------- #

# now, let's extract the home range
hares.kUDhr.90 <- getverticeshr(hares.kUD, percent = 90)
hares.kUDhr.50 <- getverticeshr(hares.kUD, percent = 50)
writeOGR(hares.kUDhr.90, "output/Shapefiles", "hares.kudhr.90", overwrite = TRUE, driver = "ESRI Shapefile")
writeOGR(hares.kUDhr.50, "output/Shapefiles", "hares.kudhr.50", overwrite = TRUE, driver = "ESRI Shapefile")

# loop through estUDm data and make each collar a raster layer and combine into a raster brick
vUDBrick <- brick(lapply(hares.vUD, function(x) try(raster(x))))
writeRaster(vUDBrick,"output/vUDRaster.grd", format="raster", overwrite = TRUE)

# plot the kernel area percent home range area ~ isopleth volume to see the type of curve
# as per Vander Wal &  Rodgers (2012), should be exponential
# melt the hrArea dataset 
hrAreamelt <- melt(hrArea, id.vars = "Kernel", variable.name = "CollarFrequency",value.name = "Area")
ggplot(hrAreamelt, aes(x=Kernel, y=Area))+
  geom_point()
# standardize the home range area proportional to the total area covered by the UD 
# and display as a percentage
hrAreastandard <- hrArea %>%
  mutate_at(vars(starts_with("X")), function(x){x/x[16]})
hrAreastandardmelt <- melt(hrAreastandard, id.vars = "Kernel", variable.name = "CollarFrequency", value.name = "StandardizedArea")
ggplot(hrAreastandardmelt, aes(x=Kernel, y=StandardizedArea))+
  geom_point()
# average the value area at each volume and see relationship
hrAreamean <- hrAreastandard %>%
  mutate(AreaMean = rowMeans(dplyr::select(hrAreastandard,starts_with("X"))))
ggplot(hrAreamean, aes(x=Kernel, y=AreaMean))+
  geom_point()
# relationship follows the exponential curve that you would expect as per Vander Wal & Rodgers (2012)
# curve is not very distinct 

# going to try an AKDE method that is better for small sample sizes
# output telemetry data so it can be uploaded to MoveBank
# convert UTM to lat/long for MoveBank
harestriangulated <- spTransform(hares.triangd, CRS("+proj=longlat +datum=WGS84"))
harestriangulated.df <- as.data.frame(harestriangulated)
harestriangulated.df$Time <- times(harestriangulated.df$Time)
write.csv(harestriangulated.df, "output/harestriangulated.csv", fileEncoding = "UTF-8")

# --------------------------------------- #
#           Visualize Kernels             #
# --------------------------------------- #

# loop through estUDm and make each collar a data frame to plot in ggplot
v_df <- lapply(hares.vUD, function(x) try(as.data.frame.estUD(x)))
# bind list of data frame 
v_df <- bind_rows(v_df, .id = "column_label")

# convert triangulated data to data frame to plot 
hares.triangd.df <- as_tibble(hares.triangd)
hares.triangd.df <- hares.triangd.df %>% mutate(Frequency = as.factor(Frequency))

# load complexity sampling points
bl_cs_pts <- read_sf("input/Mapping", layer = "cs_points")
bl_cs_pts <- st_transform(bl_cs_pts, crs = st_crs(vaancn))

bl_grid_pts <- read_sf("input/Mapping", layer = "bl_grid_points")
bl_grid_pts <- st_transform(bl_grid_pts, crs = st_crs(vaancn))
gridcoords <- as_tibble(st_coordinates(bl_grid_pts))
bl_grid_pts <- add_column(bl_grid_pts, X = gridcoords$X, Y = gridcoords$Y)

# plot contours
ggplot(data = hares.triangd.df, aes(x = X, y = Y)) +
  stat_density_2d(aes(group = Frequency, fill = stat(nlevel)), geom = "polygon", alpha = 0.15) + 
  scale_fill_viridis_c() +
  geom_point(aes(x = POINT_X_x, y = POINT_Y_y), bl_cs_pts)+
  geom_point(aes(x = X, y = Y),shape = 2, bl_grid_pts)+
  theme(
    legend.title = element_text(size = 8),
    panel.border = element_rect(size = 1, fill = NA),
    panel.background = element_rect(fill = "white"),)+
  labs(fill = "Probability")
ggsave("graphics/heatmapindividualsgrid.png")

png("graphics/heatmapalldata.png", width = 3000, height = 3000, units = "px", res = 600)
ggplot(data = hares.triangd.df, aes(x = X, y = Y)) +
  stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon", alpha = 0.25) + 
  scale_fill_viridis_c() +
  geom_point(aes(x = POINT_X_x, y = POINT_Y_y), bl_cs_pts)+
  theme(legend.position =  'none',
        panel.border = element_rect(size = 1, fill = NA),
        panel.background = element_rect(fill = "white"),
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        axis.title = element_blank())
dev.off()
