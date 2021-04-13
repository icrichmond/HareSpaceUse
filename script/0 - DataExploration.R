# This script is for examining the variation levels in C:N and C:P on Bloomfield
# where the study took place compared to the other 3 grids in our study area 

# Author: Isabella Richmond 

#### Packages ####
easypackages::packages("tidyverse", "ggpubr", "sf", "raster")

#### Data #### 
points <- read_sf("large/SamplePoints_StDMs.shp")
cn <- raster("input/VAAN_CN.tif")
cp <- raster("input/VAAN_CP.tif")

#### Extract #### 
# extract CN and CP values for all grid sites 
pointscn <- extract(cn, points, df = TRUE)
pointscn$ID <- points$PlotName
pointscp <- extract(cp, points, df = TRUE)
pointscp$ID <- points$PlotName
# join dataframes together 
cncp <- inner_join(pointscn, pointscp, by = "ID")
cncp <- rename(cncp, PlotName = ID)
full <- inner_join(cncp, points, by = "PlotName")
full <- drop_na(full)
full$SiteName <- as.factor(full$SiteName)

#### Variation #### 
# see if there are significant differences in variation across groups 
cngg <- ggboxplot(full, x = "SiteName", y = "VAAN_CN", add = "jitter", ylab = "C:N", xlab = FALSE)+
  rremove("x.text") + 
  rremove("x.ticks")
cpgg <- ggboxplot(full, x = "SiteName", y = "VAAN_CP", add = "jitter", ylab = "C:P", xlab = FALSE)
ggarrange(cngg, cpgg, ncol = 1, nrow = 2, common.legend = T)
ggsave("graphics/GridVariation.png", dpi = 400)

aovcn <- aov(VAAN_CN ~ SiteName, data = full)
TukeyHSD(aovcn)

aovcp <- aov(VAAN_CP ~ SiteName, data = full)
TukeyHSD(aovcp)

# Bloomfield is significantly different than the other grids