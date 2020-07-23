# plot complexity sampling and home ranges 
png("graphics/mapoverlap.png", width = 3000, height = 3000, units = "px", res = 600)
ggplot(hares.kUDhr.50.df, aes(long, lat))+
  geom_polygon(aes(group = id, fill = id), alpha = 0.30)+
  geom_point(aes(x = POINT_X_x, y = POINT_Y_y),bl_cs_pts)+
  coord_equal()+
  scale_fill_hue() +
  labs(y = "Latitude", x = "Longitude")+
  theme(legend.position = "none")

# make the diagnostic plots
# first plot is a residual vs fitted plot
p1 <- ggplot(overlapmodel, aes(.fitted, .resid))+geom_point()
p1 <- p1+stat_smooth(method="loess")+geom_hline(yintercept=0, col="red", linetype="dashed")
p1 <- p1+xlab("Fitted values")+ylab("Residuals")
p1 <- p1+ggtitle("Residual vs Fitted Plot")+theme_bw()+theme(plot.title = element_text(size=9),axis.title = element_text(size=8))
# second plot is a qqplot to test normality
p2 <- ggplot(overlapmodel, aes(sample=.stdresid))+stat_qq()
p2 <- p2+geom_qq_line(col='red')+xlab("Theoretical Quantiles")+ylab("Stdized Residuals")
p2 <- p2+ggtitle("Normal Q-Q")+theme_bw()+theme(plot.title = element_text(size=9),axis.title = element_text(size=8))
# third plot is a Cook's distance plot to assess outliers 
p3 <- ggplot(overlapmodel, aes(seq_along(.cooksd), .cooksd))+geom_bar(stat="identity", position="identity")
p3 <- p3+xlab("Obs. Number")+ylab("Cook's distance")
p3 <- p3+ggtitle("Cook's distance")+theme_bw()+theme(plot.title = element_text(size=9),axis.title = element_text(size=8))
# last plot is a check of influential data points, similar to Cook's distance
p4 <- ggplot(overlapmodel, aes(.hat, .stdresid))+geom_point(aes(size=.cooksd), na.rm=TRUE)
p4 <- p4+stat_smooth(method="loess", na.rm=TRUE)
p4 <- p4+xlab("Leverage")+ylab("Stdized Residuals")
p4 <- p4+ggtitle("Residual vs Leverage Plot")
p4 <- p4+scale_size_continuous("Cook's Distance", range=c(1,5))
p4 <- p4+theme_bw()+theme(plot.title = element_text(size=9),axis.title = element_text(size=8), legend.position="bottom")
# plot with patchwork and save
png("graphics/overlapdiagnostics.png", width = 3000, height = 3000, units = "px", res = 600)
(p1 | p2) /
  (p3 | p4) +
  patchwork::plot_annotation(title = paste("Diagnostic plots", "Overlap Model"))
dev.off()

# plot relationships
# make formulas of each relationship 
cover <- 
  
  png("graphics/modelrelationships.png", width = 3000, height = 5000, units = "px", res = 600)
cover <- ggplot(overlapdata, aes(x = CoverValue, y = overlap))+geom_point(color=rgb(35,77,32, maxColorValue = 255))+
  geom_smooth(method="lm", color = rgb(201,223,138, maxColorValue=255))+
  stat_cor(aes(label = paste(..rr.label..)), color = rgb(35,77,32, maxColorValue = 255))+
  theme(plot.background = element_rect(fill = rgb(119,171,89, maxColorValue = 255)))+
  theme(panel.background = element_rect(fill = rgb(240,247,218, maxColorValue = 255),colour =  rgb(240,247,218, maxColorValue = 255)))+
  theme(axis.title = element_text(colour = rgb(240,247,218, maxColorValue = 255)))+
  theme(axis.text = element_text(colour = rgb(240,247,218, maxColorValue = 255)))+
  labs(x = "Canopy Closure", y = "")
stoich <- ggplot(overlapdata, aes(x = VAAN_CN, y = overlap))+geom_point(color=rgb(35,77,32, maxColorValue = 255))+
  geom_smooth(method="lm", color = rgb(201,223,138, maxColorValue=255))+
  stat_cor(aes(label = paste(..rr.label..)), color = rgb(35,77,32, maxColorValue = 255))+
  theme(plot.background = element_rect(fill = rgb(119,171,89, maxColorValue = 255)))+
  theme(panel.background = element_rect(fill = rgb(240,247,218, maxColorValue = 255),colour =  rgb(240,247,218, maxColorValue = 255)))+
  theme(axis.title = element_text(colour = rgb(240,247,218, maxColorValue = 255)))+
  theme(axis.text = element_text(colour = rgb(240,247,218, maxColorValue = 255)))+
  labs(x = "Lowland Blueberry C:N", y = "")
hc <- ggplot(overlapdata, aes(x = meanhc, y = overlap))+geom_point(color=rgb(35,77,32, maxColorValue = 255))+
  geom_smooth(method="lm", color = rgb(201,223,138, maxColorValue=255))+
  stat_cor(aes(label = paste(..rr.label..)), color = rgb(35,77,32, maxColorValue = 255))+
  theme(plot.background = element_rect(fill = rgb(119,171,89, maxColorValue = 255)))+
  theme(axis.title = element_text(colour = rgb(240,247,218, maxColorValue = 255)))+
  theme(panel.background = element_rect(fill = rgb(240,247,218, maxColorValue = 255),colour =  rgb(240,247,218, maxColorValue = 255)))+
  theme(axis.title.y = element_text(size = 13))+
  theme(axis.text = element_text(colour = rgb(240,247,218, maxColorValue = 255)))+
  labs(x = "Horizontal Complexity", y = "Core Area Overlap")
(cover/hc/stoich)
dev.off()



# read in raster list
kernel95 <- readRDS("large/akderasters.rds")
# removing 149.394 because it is too large to be biologically relevant
kernel95[[7]] <- NULL
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
# set extent to include all home ranges but not be too large
e <- raster::extent(278000,280000,5359000,5360500)
# plot merged data - individuals 
tm_shape(kernel95normm, bbox=e)+
  tm_raster(title = "kUD", style = "cont", 
            palette = "-RdYlBu", alpha=0.9)+
  tm_scale_bar()+
  tm_layout(legend.bg.color = "white")+
  tm_grid()+
tm_shape(bl_cs_pts)+
  tm_dots(size = 0.15)
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
