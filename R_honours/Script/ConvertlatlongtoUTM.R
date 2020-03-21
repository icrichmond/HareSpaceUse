# import alpha diversity ground cover data set to R

# change coordinate system from lat long to UTM to work in arcmap
coordinates(Alpha_D) <- c("X_dd","Y_dd")


proj4string(Alpha_D) <- CRS("+proj=longlat +datum=WGS84")

Alpha_D <- spTransform(Alpha_D, CRS("+proj=utm +zone=22 ellps=WGS84"))
Alpha_D <- as.data.frame(Alpha_D)
Alpha_D$UTMZone <- as.integer(22)


write.csv(Alpha_D, "output/Alpha_D_utm.csv")
