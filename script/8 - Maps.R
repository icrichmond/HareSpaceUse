# Author: Isabella Richmond
# Last Edited: October 30, 2020

# This script is for creating figures/maps for this dataset

# load required packages 
easypackages::packages("matrixStats", "tidyverse", "lubridate", "sf", "sp", "raster",
                       "tmap", "spatialEco", "grid", "osmdata", "tmaptools")
# --------------------------------------- #
#              Home Ranges                #
# --------------------------------------- #
# read in home range rasters that were clamped in script 6 to set values 
# under 0.15 to NAs so that they could all be visualized together
clam <- readRDS("large/rasternormclamp.rds")
# stack clamped rasters and merge them - can see individuals when plotted 
clams <- raster::stack(clam)
kernel95normm <- raster::merge(clams)
# stack original raster list and overlay, summing values where there 
# is more than one raster - shows collective when plotted
kernel95norms <- raster::stack(kernel95normm)
over <- raster::overlay(kernel95norms, fun=sum)
# read in complexity sampling points
bl_cs_pts <- read_sf("input/Mapping", layer = "cs_points")
bl_cs_pts <- sf::st_transform(bl_cs_pts, "EPSG:32622")
# read in trap sampling points
bl_grid_pts <- read_sf("input/Mapping", layer = "bl_grid_points")
bl_grid_pts <- sf::st_transform(bl_grid_pts, "+init=epsg:32622")
# read in stoich data 
vaancn <- raster("input/VAAN_CN.tif")
vaancp <- raster("input/VAAN_CP.tif")
# set extent to include all home ranges but not be too large
e <- raster::extent(860000, 863000, 5383000, 5386000)
# crop rasters
vaancnclip <- crop(vaancn, e)
vaancpclip <- crop(vaancp, e)
# reproject rasters 
vaanextent <- projectExtent(vaancn, crs="+init=epsg:32622")
vaancnclip <- projectRaster(vaancnclip, vaanextent)
vaancpclip <- projectRaster(vaancpclip, vaanextent)

# NL boundary layer (code from https://github.com/wildlifeevoeco/study-area-figs)
zz <- opq(getbb('Newfoundland')) %>%
  add_osm_feature(key = 'place', value = 'island') %>%
  osmdata_sf()
lns <- zz$osm_lines
# union -> polygonize -> cast lines = geometry set
castpolys <- st_cast(st_union(lns))
# Combine geometries and cast as sf
nl <- st_as_sf(castpolys)
nl <- st_set_crs(nl, "+init=epsg:4326")
# reproject
utmNL <- sf::st_transform(nl, "+init=epsg:32622")
saveRDS(utmNL, "large/NLboundaryUTM.rds")

# get the min x and y and use them to calc aspect ratio
xy <- st_bbox(utmNL)
asp <- (xy$xmax - xy$xmin)/(xy$ymax - xy$ymin)

# plot merged data - individuals 
p <- get_brewer_pal("Greys", n = 10, contrast = c(0.3, 1))
e <- raster::extent(278000,280000,5359000,5360500)

t <- tm_shape(kernel95normm, bbox=e)+
  tm_raster(title = "kUD", style = "cont", 
            palette = p, alpha=0.9)+
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

# create map t with inset of Newfoundland 
# inset map
inset <- tm_shape(utmNL)+
  tm_lines()+
  tm_layout(outer.margins=c(0,0,0,0))+
  tm_shape(bl_grid_pts)+
  tm_dots(size=0.3, col="red")
w <- 0.29
h <- asp*w

print(inset, vp=viewport(x=0.89,y=0.82,width=w,height=h))
# save
tmap_save(t, insets_tm = inset, insets_vp = viewport(x=0.89,y=0.83,width=w,height=h), filename="graphics/kUD_raster_grids.pdf")




# plot overlay data - populations
tm_shape(over, bbox=e)+
  tm_raster(title = "kUD", style = "cont", 
            palette = "-RdYlBu", alpha=0.9)+
  tm_scale_bar()+
  tm_layout(legend.bg.color = "white")+
  tm_grid()+
  tm_shape(bl_cs_pts)+
  tm_dots(size = 0.15)


# --------------------------------------- #
#   Food Quality & Predation Risk         #
# --------------------------------------- #
# plot stoich and predation risk variation 
# make predrisk spatial
# import predation risk dataset with PCA axes data
# note - had to drop some of the complexity sampling points because 
# they had NA values (doesn't work in a PCA), there are 67 instead of 72
predrisk <- read_csv("output/predationriskpca.csv")
# make predation risk data spatially explicit by associating coordinate data 
predriskspatial <- inner_join(predrisk, bl_cs_pts, by = "Plot")
predriskspatial <- st_as_sf(predriskspatial)
predriskspatial <- as_Spatial(predriskspatial)

# plot
cnover <- tm_shape(vaancnclip, bbox=e)+
  tm_raster(palette = "BrBG", title="C:N")+
  tm_scale_bar()+
  tm_layout(legend.bg.color = "white", legend.title.size = 1, legend.outside = FALSE, 
            legend.position = c("left", "top"), legend.stack = "horizontal")+
  tm_grid()+
  tm_shape(predriskspatial)+
  tm_bubbles(size=0.2, col="overPCA", palette="Greys", title.col="Overstory", midpoint=NA)
#tmap_save(cnover, "graphics/OverstoryComplexity_CN.png")

cpunder <- tm_shape(vaancpclip, bbox=e)+
  tm_raster(palette = "BrBG", title = "C:P")+
  tm_scale_bar()+
  tm_layout(legend.bg.color = "white", legend.title.size=1, legend.outside=FALSE,
            legend.position = c("left", "top"), legend.stack="horizontal")+
  tm_grid()+
  tm_shape(predriskspatial)+
  tm_bubbles(size=0.2, col="underPCA", palette="Greys", title.col="Understory", midpoint=NA)
#tmap_save(cpunder, "graphics/UnderstoryComplexity_CP.png")

pan <- tmap_arrange(cnover, cpunder)
tmap_save(pan, "graphics/CNCPoverunder.png")