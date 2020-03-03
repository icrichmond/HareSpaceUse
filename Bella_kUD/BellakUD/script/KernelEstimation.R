# This script is for estimating the space use of snowshoe hare collared in Bloomfield, NL 
# over three years (2016-2019). Relocations only taken in summer season.

# Author: Isabella Richmond (code and data shared between Matteo Rizzuto: github.com/matteorizzuto)
# Last Edited: March 3, 2020

# load required packages 
easypackages::packages("sp", "maptools", "tmap", "tmaptools", 
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
#          kUD Estimation                 #
# --------------------------------------- #
# This code taken directly from Matteo Rizzuto's NewCollarLocsEstimation.R script
# Find Matteo's repository at github.com/matteorizzuto/Chapter_2


