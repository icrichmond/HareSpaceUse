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





# replace with your folder name:
dir <- "output/Shapefiles/aKDE_Home"
ff <- list.files(dir, pattern="\\.shp$", full.names=TRUE)
x <- lapply(ff, shapefile)
lapply(1:length(x), function(i) crs(x[[i]]) <- crs(vaancn))

crs(x[[1]]) <- crs(vaancn)
