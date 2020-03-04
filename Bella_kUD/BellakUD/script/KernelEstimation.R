# This script is for estimating the space use of snowshoe hare collared in Bloomfield, NL 
# over three years (2016-2019). Relocations only taken in summer season.

# Author: Isabella Richmond (code and data shared between Matteo Rizzuto: github.com/matteorizzuto)
# Last Edited: March 3, 2020

# load required packages 
easypackages::packages("sp", "sf", "maptools", "tmap", "tmaptools", 
              "adehabitatHR", "adehabitatHS", "adehabitatLT", "ellipse", 
              "nleqslv", "adehabitatMA", "adehabitatHR","dplyr", "gdtools", "ggmap", "ggplot2", 
              "ggrepel", "ggsci", "ggthemes", "maps", "raster", "spatial", 
              "tidyr", "readr","rgdal", "rgeos", "reshape2", "dismo", "tibble")


### NOTE: one package (sigloc) is no longer maintained - need to manaully load package ###
# Package can still be found on CRAN archive: https://cran.r-project.org/src/contrib/Archive/sigloc/


# --------------------------------------- #
#             Data Preparation            #
# --------------------------------------- #

# import telemetry datasetss
VHF2019 <- read.csv("input/TelemetryPoints_VHF_2019.csv")
View(VHF2019)
VHF2018 <- read.csv("input/VHF_CleanData_2018.csv")
View(VHF2018)
VHF2017 <- read.csv("input/VHF_CleanData_2017.csv")
View(VHF2017)

# 2019 uses lat/long instead of UTM
# convert 2019 to UTM for simplicity 

coordinates(VHF2019) <- c("Easting", "Northing")
proj4string(VHF2019) <- CRS("+proj=longlat +datum=WGS84")

VHF2019 <- spTransform(VHF2019, CRS("+proj=utm +zone=22 ellps=WGS84"))
VHF2019 <- as.data.frame(VHF2019)
VHF2019$UTMZone <- as.integer(22)

# following section developed from Matteo Rizzuto's (github.com/matteorizzuto) code
# process the dataset so it is ready for home range estimation.
# Matteo Rizzuto also using this data for his manuscript. Code can be found in his repository on GitHub

# isolate live collars to make working on them easier
Data17 <- subset(VHF2017, VHF2017$Frequency == "149.535" | 
                   VHF2017$Frequency == "149.673" |
                   VHF2017$Frequency == "149.394" |
                   VHF2017$Frequency == "149.452")

Data18 <- subset(VHF2018, VHF2018$Frequency == "149.003" |
                   VHF2018$Frequency == "149.053" |
                   VHF2018$Frequency == "149.093" |
                   VHF2018$Frequency == "149.173" |
                   VHF2018$Frequency == "149.213" |
                   VHF2018$Frequency == "149.274" |
                   VHF2018$Frequency == "149.374" |
                   VHF2018$Frequency == "149.474" |
                   VHF2018$Frequency == "149.613" |
                   VHF2018$Frequency == "149.633" |
                   VHF2018$Frequency == "149.653")

# drop data from 18-06-2018 for collar 149.653 
Data18 <- Data18[!(Data18$Frequency == "149.653" & Data18$Date == "2018-06-18"),]


Data19 <- subset(VHF2019, VHF2019$ï..Frequency == "149.124" |
                   VHF2019$ï..Frequency == "149.233" | 
                   VHF2019$ï..Frequency == "149.294" | 
                   VHF2019$ï..Frequency == "149.423" |
                   VHF2019$ï..Frequency == "149.513" |
                   VHF2019$ï..Frequency == "149.555" |
                   VHF2019$ï..Frequency == "149.594" |
                   VHF2019$ï..Frequency == "150.032" |
                   VHF2019$ï..Frequency == "150.052" |
                   VHF2019$ï..Frequency == "150.072" |
                   VHF2019$ï..Frequency == "150.091" |
                   VHF2019$ï..Frequency == "150.111" |
                   VHF2019$ï..Frequency == "150.132" |
                   VHF2019$ï..Frequency == "150.154" | 
                   VHF2019$ï..Frequency == "150.173" | 
                   VHF2019$ï..Frequency == "150.191" | 
                   VHF2019$ï..Frequency == "150.232" | 
                   VHF2019$ï..Frequency == "150.273" |
                   VHF2019$ï..Frequency == "150.314" |
                   VHF2019$ï..Frequency == "150.332" |
                   VHF2019$ï..Frequency == "150.373" |
                   VHF2019$ï..Frequency == "150.392")

# remove NAs
Data17 <- drop_na(Data17, Azimuth)
Data18 <- drop_na(Data18, Azimuth)
Data19 <- drop_na(Data19, Azimuth)

# add a year variable 
Data17 <- add_column(Data17, Year="2017")
Data18 <- add_column(Data18, Year="2018")
Data19 <- add_column(Data19, Year="2019")

# Categorize times as AM or PM in Data19

for (i in 1:(nrow(Data19) - 1)) {
  Data19$AM_PM <- ifelse(Data19$TimeCategory == "Morning","AM", "PM")
}

# rename collar heading in Data19 to match Data18 and Data17
Data19 <- dplyr::rename(Data19, Frequency = ï..Frequency)

# remove unused columns from the dataframes
Data17 <- dplyr::select(Data17, -c(Notes, EarTag, Wind, Rain, SampleTimeCat, SampleTimeGeneral, Clouds, Temp, Line, Alive, Time_O))
Data18 <- dplyr::select(Data18, -c(Notes, EarTag, Wind, Rain, SampleTimeCat, SampleTimeGeneral, Clouds, Temp, Line, Alive, Time_O))
Data19 <- dplyr::select(Data19, -c(FixLocation, Time_O, TimeCategory, TempC, Wind, Rain))

# combine all three years into one dataframe 
liveData <- rbind(Data17, Data18, Data19)
liveData$Frequency <- as.factor(liveData$Frequency)
liveData$Year <- as.factor(liveData$Year)
# remove NAs
liveData <- drop_na(liveData, Azimuth)
# remove second observer, who only did 1 day of triangulations 
# and no error reduction trials
liveData <- subset(liveData, liveData$Observer != "TH")

# The package sigloc requires identification of each individual set of 
# triangulation azimuths via a Group Identifier (GID). The following code 
# assigns a unique GID to each unique combination of Frequency and Date.

IDs <- paste(liveData[, "Frequency"])   # isolate the IDs of each data entry

UniqIDs <- unique(IDs) # extract unique ID values

OutputFields <- c(paste(colnames(liveData)), "GID") # create column names for new df

liveData$GID <- NA    # add GID column to original data
liveData$GID <- as.factor(liveData$GID)

vhfData <- liveData[FALSE,] # create new df

liveData$GID <- NULL # remove useless GID column from original dataset

# The following loop scans the data for groups of azimuths 
# obtained for each collar on each day of triangulation, and assigns a GID 
# based on the date of collection.

for (i in 1:length(UniqIDs)){
  # first, create a temporary dataset to host each Frequency's data in turn
  TmpInds <- which(IDs == UniqIDs[i]) 
  TmpData <- liveData[TmpInds,]   
  TmpData[1, "GID"] <- 1  # assign a starting GID
  for (j in 1:(nrow(TmpData) - 1)){
    if (TmpData$Year[j] == "2017" | TmpData$Grid[j] == "UNI"){
      # if the Date stays the same, we want to assign the same GID 
      if (TmpData[j + 1, "Date"] == TmpData[j, "Date"] &&
          TmpData[j + 1, "AM_PM"] == TmpData[j, "AM_PM"]){
        TmpData[j + 1, "GID"] <- TmpData[j, "GID"]
      }
      else {
        # if the Date changes, we want the GID to increase accordingly
        TmpData[j + 1, "GID"] <- TmpData[j, "GID"] + 1
      }
    }
    else {
      if (TmpData[j + 1, "Date"] == TmpData[j, "Date"]){
        TmpData[j + 1, "GID"] <- TmpData[j, "GID"]
      }
      else {
        # if the Date changes, we want the GID to increase accordingly
        TmpData[j + 1, "GID"] <- TmpData[j, "GID"] + 1
      }
    }
  }
  # browser()
  vhfData <- rbind(vhfData, TmpData) # bind each temp dataset into a new one
}

vhfData$GID <- as.factor(vhfData$GID)
vhfData$Time <- as.numeric(vhfData$Time)



# --------------------------------------- #
#       Collar Location Estimation        #
# --------------------------------------- #
# This code taken directly from Matteo Rizzuto's NewCollarLocsEstimation.R script
# Find Matteo's repository at github.com/matteorizzuto/Chapter_2

vhfData <- vhfData[with(vhfData, order(Year, Frequency)),]

# allocate an empty list to store each hare's triangulation results
hares.list <- vector("list", length(UniqIDs))

# run the triangulation loop
for (i in 1:length(UniqIDs)) {
  # we want to keep track of what R is doing, so let's have it print a message
  # telling us exactly which collar it is processing
  print(paste("Starting to triangulate collar", UniqIDs[i]))
  
  # Following a suggestion from sigloc's author, Simon Berg, we run the locate()
  # and plot() functions one GID at the time, for each collar.The for loop below 
  # takes care of this.
  
  # create a working dataset containing only data for the current collar
  test <- subset(vhfData, vhfData$Frequency == UniqIDs[i])
  
  # allocate vector along whose length to run the nested loop
  test.index <- unique(test$GID) 
  
  # allocate an empty dataframe to store the results of the triangulation
  workingCollar.loc <- data.frame("X" = numeric(), "Y" = numeric(), 
                                  "Badpoint" = integer(), "Var_X" = numeric(), 
                                  "Var_Y" = numeric(), "Cov_XY" = numeric(), 
                                  "AngleDiff" = numeric(), "Date" = character(), 
                                  "Time" = integer(), stringsAsFactors = FALSE)
  
  # set up plotting
  # frame()
  par(mar = c(2, 2, 2, 2), mfrow = c(6,6))
  
  for (j in 1:length(test.index)) {
    # separate a day's set of azimuths from the rest of the current working dataset
    test.dat <- subset(test, test$GID == j)
    
    # Prevents the looped code from looking for earlier GIDs
    test.dat$GID <- 1 
    
    # save the triangulation date as a vector for later use
    date <- as.character(unique(test.dat$Date))
    
    # set the current test dataframe as a receiver for triangulation
    test.rec <- as.receiver(test.dat)
    
    # find the intersection points of the three azimuths in the current test df
    test.int <- findintersects(test.rec)
    
    # visualize triangulation
    plot(test.rec, bearings = TRUE, xlab = "Easting", ylab = "Northing", asp = 1, 
         ylim = c(5359000, 5360000), 
         xlim = c(278500, 279900))
    title(main = unique(test.rec$Date), sub = unique(test.dat$Frequency))
    
    # locate the transmitting collar using a Maximum Likelihood Estimator
    test.loc <- locate(test.rec)
    
    # add the date to the localized collar
    test.loc$Date <- as.character(date) 
    
    # plot the localized collar on the same plot as above;
    # if the point is red, it means the MLE failed to retun a point and the 
    # triangulation was complete by taking the midpoint of the azimuths' 
    # intersections to estimate the location of the collar
    plot(test.loc, add = TRUE, badcolor = TRUE, errors = TRUE)
    
    # save subsequent triangulations in the current workingCollar.loc dataframe 
    workingCollar.loc <- rbind(workingCollar.loc, test.loc)
    
    # make sure the date is store as a character rather than a number
    workingCollar.loc$Date <- as.character(workingCollar.loc$Date)
    # browser()
    # Sys.sleep(1)
  }
  
  # add Frequency identifier to the current workingCollar.loc dataframe
  workingCollar.loc$Frequency <- as.character(UniqIDs[i])
  # don't have EarTag for 2019, change this line of code to Frequency
  workingCollar.loc$Frequency <- as.character(unique(test$Frequency))
  # browser()
  
  # Let's take a look at these points! In order to plot them as spatial points, 
  # we first need to convert them to a Spatial Object. To do so, first we create
  # a vector containing just the lat and long coordinates of the traingulated 
  # points
  
  workingCollar.coords <- workingCollar.loc[,c(1,2)]
  
  # then we remove the corresponding columns from the dataframe
  
  workingCollar.loc[,c(1,2)] <- NULL
  
  # finally, using package sp, we combine the two to obtain a spatial object
  
  workingCollar.spatial <- SpatialPointsDataFrame(coords = workingCollar.coords, 
                                                  data = workingCollar.loc, 
                                                  proj4string = CRS("+proj=utm +zone=22 +datum=NAD83"))
  
  # plotting
  # frame()
  par(mfrow=c(1,1))
  # very simple plot, just to check everything works
  plot(workingCollar.spatial, main = unique(workingCollar.spatial$Frequency))  
  # add identifiers for each relocation
  text(workingCollar.spatial@coords, labels = workingCollar.spatial$Date, pos = 2, cex = 0.7)
  
  # store the newly produced set of triangulation for collar i into the 
  # previously allocated hares.list
  hares.list[[i]] <- workingCollar.spatial
  
  if (is.na(UniqIDs[i+1])==TRUE) {
    print(paste("Finished triangulating collar", UniqIDs[i], ";", "Done"))
  } else {print(paste("Finished triangulating collar", UniqIDs[i],";", "moving on to collar", UniqIDs[i+1]))
  }
  
  # browser()
}

# convert the list of triangulated collars dataframes into a single dataframe 
# for further analyses
hares.triangd <- do.call("rbind", hares.list)

# manually went through and removed problematic points from the dataset using Excel 
# check with Matteo on this 
# input/HomeRanges.txt file describing which points are removed
# check with Matteo that datasets are the same and all problem points have been removed

# remove the two collars that are not in Bloomfield
hares.triangd <- subset(hares.triangd, hares.triangd$Frequency != "149.374" &
                          hares.triangd$Frequency != "149.474")


# --------------------------------------- #
#         Load Background Layers          #
# --------------------------------------- #
# This code taken directly from Matteo Rizzuto's HomeRangeEstimation.R script
# Find Matteo's repository at github.com/matteorizzuto/Chapter_2

# read the forest shapefile from the Land Cover geodatabase
forest_gdb <- read_sf("input/Mapping", layer = "Forest.shp")
# and set it to the same crs as the triangulation data
forest_gdb <- st_transform(forest_gdb, crs = st_crs(hares.triangd))

# load grid trap locations shapefile
bl_grid_pts <- st_transform(bl_grid_points, crs = st_crs(hares.triangd))

# now, let's intersect the forest
forest_clip <- st_intersection(bl.mcp.all500, forest_gdb)
rm(forest_gdb) # remove the large FRI shapefile to make spatial R go faster
tm_shape(forest_clip) + tm_polygons(alpha = 0) + tm_shape(bl_grid_pts) + 
  tm_dots(size = 0.3, shape = 3, col = "red") + tm_compass() + tm_scale_bar()


# --------------------------------------- #
#           Create Collar List            #
# --------------------------------------- #
# This code taken directly from Matteo Rizzuto's HomeRangeEstimation.R script
# Find Matteo's repository at github.com/matteorizzuto/Chapter_2

# create list of names of collars from Bloomfield study site for mapping each
# BL set of points and future mapping/listing use
UniqCIDs <- UniqIDs[UniqIDs!="149.374"]
UniqCIDs <- UniqCIDs[UniqCIDs!="149.474"]

# allocate empty list to store tmap objects for each map created below
hares.dots.maps <- vector("list", length(UniqCIDs))

# loop to produce a map for each hare showing a heatmap of the kUD and the 
# 50% and 90% isopleths
for (i in 1:length(UniqCIDs)) {
  # store contours for the current collar in a temporary linear object
  dots.temp <- subset(hares.triangd, hares.triangd$Frequency == UniqCIDs[i]) 
  
  # add a crs to the temporary linear object
  crs(dots.temp) <- crs(forest_clip) 
  
  bbox_temp <- st_bbox(dots.temp) # current bounding box
  
  xrange <- bbox_temp$xmax - bbox_temp$xmin # range of x values
  yrange <- bbox_temp$ymax - bbox_temp$ymin # range of y values
  
  bbox_temp[1] <- bbox_temp[1] - (0.05 * xrange) # xmin - left
  bbox_temp[3] <- bbox_temp[3] + (0.05 * xrange) # xmax - right
  bbox_temp[2] <- bbox_temp[2] - (0.05 * yrange) # ymin - bottom
  bbox_temp[4] <- bbox_temp[4] + (0.05 * yrange) # ymax - top
  
  bbox_temp <- bbox_temp %>%  # take the bounding box ...
    st_as_sfc() # ... and make it a sf polygon
  
  # save the number of levels in Date for later plotting
  col.lvls <- nlevels(as.factor(dots.temp$Date)) 
  
  # create a map using tmap and store it in its own slots in the list created
  # above
  hares.dots.maps[[i]] <- tm_shape(forest_clip) + 
    tm_borders(col = "grey") +
    tm_shape(bl_grid_pts) + tm_dots(size = 0.3, 
                                    shape = 3, col = "red") +
    tm_shape(dots.temp, bbox = bbox_temp) + 
    tm_dots(size = 0.5, shape = 21, col = "Date",
            palette = "viridis", 
            border.col = "grey") + 
    # tm_text("Date", auto.placement = TRUE, size = 0.75, xmod = 1, ymod = 0.5) +
    tm_compass(position = c("left", "top"), text.size = 1) +
    tm_scale_bar(text.size = 0.8,
                 position = c("right", "bottom")) + 
    tm_layout(main.title = paste("Collar Frequency:", UniqCIDs[i]),
              title = paste("Ear tag:", dots.temp$EarTag),
              legend.outside = TRUE, 
              legend.outside.position = "right",
              legend.text.size = 1.1,
              legend.title.size = 1.5,
              asp = 1,
              frame = FALSE) +
    tmap_options(max.categories = col.lvls)
  
  # Sys.sleep(1)
  # browser()
  
  print(hares.dots.maps[[i]])
  
  rm(dots.temp) # remove the temp dots object to avoid errors
}

# tmap_arrange(hares.dots.maps)



# --------------------------------------- #
#         Home Range Estimation           #
# --------------------------------------- #
# This code taken directly from Matteo Rizzuto's HomeRangeEstimation.R script
# Find Matteo's repository at github.com/matteorizzuto/Chapter_2

# before estimating the kUD, remove collars 149.555, 150.132 due to too few relocations 
# available for estimating a reliable kernel Utilization Distribution
# CHECK WITH MATTEO - WHAT DID HE REMOVE?
hares.triangd <- subset(hares.triangd, hares.triangd$Frequency != "149.555" 
                         & hares.triangd$Frequency != "150.132")

# Let's estimate the kernel Utilization Distribution using the ad hoc method and 
# a grid of a 1000 points that can adapt to the general geographic area used by
# each animal
hares.kUD <- kernelUD(hares.triangd[,8], h = 'LSCV', grid = 1000, extent = 1, same4all = FALSE)

# If reverting back to using LSCV to estimate h, double-check that minimization
# of the cross-validation criteria is successful using:
par(mar = c(1,1,1,1), mfrow = c(1,1))
plotLSCV(hares.kUD)

# Estimate kUD area using a range of percent levels
kUD.hr.estimates <- kernel.area(hares.kUD, percent = seq(50, 95, 5), 
                                unout = "ha")

# and extract values only for the core area to be used in later modelling
hrArea.50 <- kUD.hr.estimates[1, ]
hrArea.50 <- tidyr::pivot_longer(hrArea.50, cols = 1:ncol(hrArea.50), names_to = "CollarFrequency", values_to = "HRsize", names_prefix = "X")

# now, let's extract the home range
hares.kUDhr.90 <- getverticeshr(hares.kUD, percent = 90)
hares.kUDhr.50 <- getverticeshr(hares.kUD, percent = 50)

plot(hares.kUDhr.90, col=1:34)
plot(hares.kUDhr.50, col=1:34)

# NOTE for future: need to load in stoich raster, match grid to raster, estimate kUD 
# in raster format, improve plotting