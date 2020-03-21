# import alpha d ground cover diversity for only selected forage spp data set to R

# change coordinate system from lat long to UTM to work in arcmap
coordinates(SH_Alpha_D) <- c("X_dd","Y_dd")


proj4string(SH_Alpha_D) <- CRS("+proj=longlat +datum=WGS84")

SH_Alpha_D <- spTransform(SH_Alpha_D, CRS("+proj=utm +zone=22 ellps=WGS84"))
SH_Alpha_D<- as.data.frame(SH_Alpha_D)
SH_Alpha_D$UTMZone <- as.integer(22)


write.csv(Gamma_D, "output/SH_Alpha_D_utm.csv")
