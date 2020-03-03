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

hares.list <- vector("list", length(UniqIDs))

for (i in 1:length(UniqIDs)) {
  print(paste("Starting to triangulate collar", UniqIDs[i]))
  
  # Following a suggestion from sigloc's author, Simon Berg, I will be
  # running the locate and plotting process one GID at the time, for each collar
  # (credit to S. Berg for the loop)
  # The for loop below takes care of this. Note that azimuths = 90 degrees had to
  # be changed to 91 degrees to avoid having a mathematical angle of 0, which would
  # cause the Maximum Likelihood Estimator used by locate() to return an error
  
  test <- subset(vhfData, vhfData$Frequency == UniqIDs[i])
  
  test.index <- unique(test$GID)
  
  # browser()
  
  workingCollar.loc <- data.frame("X" = numeric(), "Y" = numeric(), "Badpoint" = integer(), "Var_X" = numeric(), "Var_Y" = numeric(), "Cov_XY" = numeric(), "AngleDiff" = numeric(), "Date" = character(), "Time" = integer(), stringsAsFactors = FALSE)
  
  # str(workingCollar.loc)
  
  # browser()
  frame()
  par(mar = c(2, 2, 2, 2), mfrow = c(6,6))
  
  for (j in 1:length(test.index)) {
    
    test.dat <- subset(test, test$GID == j)
    
    # browser()
    # test.dat$GID_o <- vhfData$GID[i]
    
    test.dat$GID <- 1 # Prevents the looped code from looking for earlier GIDs
    
    date <- as.character(unique(test.dat$Date))
    
    test.rec <- as.receiver(test.dat)
    
    test.int <- findintersects(test.rec)
    
    plot(test.rec, bearings = TRUE, xlab = "Easting", ylab = "Northing", asp = 1, 
         ylim = c(5359000, 5360000), 
         xlim = c(278500, 279900))
    title(main = unique(test.rec$Date), sub = unique(test.dat$Frequency))
    
    test.loc <- locate(test.rec)
    
    test.loc$Date <- as.character(date) 
    
    plot(test.loc, add = TRUE, badcolor = TRUE, errors = TRUE)
    
    workingCollar.loc <- rbind(workingCollar.loc, test.loc)
    
    workingCollar.loc$Date <- as.character(workingCollar.loc$Date)
    # browser()
    # Sys.sleep(1)
  }
  
  workingCollar.loc$Frequency <- as.character(UniqIDs[i])
  # browser()
  
  # Let's take a look at these points! In order to plot them as spatial points, 
  # I first need to convert them to a Spatial Object. To do so, first I create
  # a vector containing just the lat and long coordinates of the points
  
  workingCollar.coords <- workingCollar.loc[,c(1,2)]
  
  # then I remove the corresponding columns from the dataframe
  
  workingCollar.loc[,c(1,2)] <- NULL
  
  # finally, using package sp, I combine the two to obtain a spatial object
  
  workingCollar.spatial <- SpatialPointsDataFrame(coords = workingCollar.coords, 
                                                  data = workingCollar.loc, 
                                                  proj4string = CRS("+proj=utm +zone=22 +datum=NAD83"))
  
  # plotting
  frame()
  par(mfrow=c(1,1))
  plot(workingCollar.spatial, main = unique(workingCollar.spatial$Frequency))  # very simple, just to check everything works
  text(workingCollar.spatial@coords, labels = workingCollar.spatial$Date, pos = 3)
  # more complex plotting (with ggplot2 etc) will take place below
  
  hares.list[[i]] <- workingCollar.spatial
  
  print(paste("Finished triangulating collar", UniqIDs[i], "moving on to collar", UniqIDs[i+1]))
  
  # browser()
}

# convert the list of triangulated collars dataframes into a single dataframe 
# for further analyses
hares.triangd <- do.call("rbind", hares.list)

# Clean up outlying relocations 
# remove triangulations 13, 20, 26, 28 June, 29 July, 29, 30 August 2019 
# for 149.233
toBeRemoved <- which(hares.triangd$Frequency =="149.233" &
                       hares.triangd$Date == "2019-06-13" | 
                       hares.triangd$Date == "2019-06-26" | 
                       hares.triangd$Date == "2019-06-28" |
                       hares.triangd$Date == "2019-07-29" |
                       hares.triangd$Date == "2019-08-29" |
                       hares.triangd$Date == "2019-08-30" | 
                       hares.triangd$Date == "2019-06-20")

hares.triangd <- hares.triangd[-toBeRemoved, ]

# remove triangulations 28 May, 20, 24, 27 June 2019 for 149.294
toBeRemoved <- which(hares.triangd$Frequency == "149.294" & 
                       hares.triangd$Date == "2019-05-28" | 
                       hares.triangd$Date == "2019-06-20" | 
                       hares.triangd$Date == "2019-06-24" |
                       hares.triangd$Date == "2019-06-27")

hares.triangd <- hares.triangd[-toBeRemoved, ]

# remove triangulations from 16, 18, 27 May, 18, 20 June 2019 for 149.423
toBeRemoved <- which(hares.triangd$Frequency == "149.423" & 
                       hares.triangd$Date == "2019-05-16" |
                       hares.triangd$Date == "2019-05-18" |
                       hares.triangd$Date == "2019-05-27" |
                       hares.triangd$Date == "2019-06-18" |
                       hares.triangd$Date == "2019-06-20")

hares.triangd <- hares.triangd[-toBeRemoved, ]

# remove triangulations from 25 May 2019 for 149.513
toBeRemoved <- which (hares.triangd$Frequency == "149.513" & 
                        hares.triangd$Date == "2019-05-25")


# remove triangulations from 19, 22, 28 May, 18, 22, 25, 26 June 2019 
# for 150.173
toBeRemoved <- which(hares.triangd$Frequency == "150.173"&
                       hares.triangd$Date == "2019-05-19" |
                       hares.triangd$Date == "2019-05-22" |
                       hares.triangd$Date == "2019-05-28" |
                       hares.triangd$Date == "2019-06-18" |
                       hares.triangd$Date == "2019-06-22" |
                       hares.triangd$Date == "2019-06-25" |
                       hares.triangd$Date == "2019-06-26")

hares.triangd <- hares.triangd[-toBeRemoved, ]

# remove collars 149.513, 149.555, 150.032, 150.052, 150.154 due to too few relocations 
# available for estimating a reliable kernel Utilization Distribution
hares.triangd <- subset(hares.triangd, 
                        hares.triangd$Frequency != "149.513" &
                          hares.triangd$Frequency != "149.555" &
                          hares.triangd$Frequency != "150.032" &
                          hares.triangd$Frequency != "150.052" &
                          hares.triangd$Frequency != "150.154")

# finally, remove the two collars that are not in Bloomfield
hares.triangd <- subset(hares.triangd, hares.triangd$Frequency != "149.374" &
                          hares.triangd$Frequency != "149.474")
