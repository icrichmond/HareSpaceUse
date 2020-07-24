



# read in raster list
kernel95 <- readRDS("large/akderasters.rds")
# removing 149.555 because it was sampled across multiple years
kernel95[[12]] <- NULL
# transform raster values so they are between 0,1
kernel95norm <- lapply(kernel95, spatialEco::raster.transformation)
# use clamp to transform any values under 0.15 to Nas
kernel95normz <- lapply(kernel95norm, function(i) raster::clamp(i, lower=0.15, useValues=FALSE))
# save clamped raster
saveRDS(kernel95normz, "large/rasternormclamp.rds")
clam <- readRDS("large/rasternormclamp.rds")
# stack clamped rasters and merge them - can see individuals when plotted 
clams <- raster::stack(clam)
kernel95normm <- raster::merge(clams)
# stack original raster list and overlay, summing values where there 
# is more than one raster - shows collective when plotted
kernel95norms <- raster::stack(kernel95norm)
over <- raster::overlay(kernel95norms, fun=sum)
# read in complexity sampling points
bl_cs_pts <- read_sf("input/Mapping", layer = "cs_points")
bl_cs_pts <- sf::st_transform(bl_cs_pts, "+init=epsg:32622")
# read in trap sampling points
bl_grid_pts <- read_sf("input/Mapping", layer = "bl_grid_points")
bl_grid_pts <- sf::st_transform(bl_grid_pts, "+init=epsg:32622")
# set extent to include all home ranges but not be too large
e <- raster::extent(278000,280000,5359000,5360500)
# plot merged data - individuals 
t <- tm_shape(kernel95normm, bbox=e)+
  tm_raster(title = "kUD", style = "cont", 
            palette = "-RdYlBu", alpha=0.9)+
  tm_scale_bar()+
  tm_grid()+
  tm_xlab("Easting")+
  tm_ylab("Northing")+
tm_shape(bl_cs_pts)+
  tm_dots(size = 0.15)+
tm_shape(bl_grid_pts)+
  tm_dots(size = 0.15, shape = 4)+
tm_legend()+
tm_layout(legend.bg.color = "white")
t
tmap_save(t, "graphics/kUD_raster_grids.png")

plot(bl_cs_pts$geometry)
plot(bl_grid_pts$geometry, add=T)

# plot overlay data - populations
tm_shape(over, bbox=e)+
  tm_raster(title = "kUD", style = "cont", 
            palette = "-RdYlBu", alpha=0.9)+
  tm_scale_bar()+
  tm_layout(legend.bg.color = "white")+
  tm_grid()+
tm_shape(bl_cs_pts)+
  tm_dots(size = 0.15)


ade <- read.csv("output/homerangeareas.csv")
ade <- filter(ade, Kernel == 95)
hr_df <- hr_df[-c(12,22),]
all <- bind_cols(ade, hr_df)

ggplot(all, aes(x=value, y=homerange))+
  geom_point()+
  labs(x = "AdeHabitat Area (ha)", y = "ctmm Area (ha)", title = "95% Home Ranges")
cor(all$value, all$homerange)
