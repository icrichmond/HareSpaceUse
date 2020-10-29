# This script is for cleaning the VHF collar bearings taken for snowshoe hares 
# over three years (2016-2019). Relocations only taken in summer season.

# Author: Isabella Richmond (code and data shared between 
# Matteo Rizzuto: github.com/matteorizzuto and Alec Robitaille: github.com/robitalec)

# Last Edited: July 8, 2020


devtools::install_github("cppeck/razimuth")
easypackages::packages("razimuth", "stringr", "data.table", "tidyverse", "chron", "sf", "sp")

# --------------------------------------- #
#             Data Preparation            #
# --------------------------------------- #

# import telemetry datasets
VHF2019 <- read.csv("input/TelemetryPoints_VHF_2019.csv")
head(VHF2019)
VHF2018 <- read.csv("input/VHF_CleanData_2018.csv")
head(VHF2018)
VHF2017 <- read.csv("input/VHF_CleanData_2017.csv")
head(VHF2017)
# remove NAs
VHF2019 <- drop_na(VHF2019, Azimuth)
VHF2018 <- drop_na(VHF2018, Azimuth)
VHF2017 <- drop_na(VHF2017, Azimuth)

# 2019 uses lat/long and 2017-2018 use UTM
# convert all of them to meters to match stoich raster CRS
VHF2017 <- drop_na(VHF2017, Easting)
coordinates(VHF2017) <- c("Easting", "Northing")
proj4string(VHF2017) <- CRS("+init=epsg:32622")

VHF2017 <- spTransform(VHF2017, CRS("+proj=tmerc +lat_0=0 +lon_0=-61.5 +k=0.9999 +x_0=304800 +y_0=0 +ellps=GRS80 +units=m
                                    +no_defs"))
VHF2017 <- as.data.frame(VHF2017)
VHF2017$UTMZone <- as.integer(22)


coordinates(VHF2018) <- c("Easting", "Northing")
proj4string(VHF2018) <- CRS("+init=epsg:32622")

VHF2018 <- spTransform(VHF2018, CRS("+proj=tmerc +lat_0=0 +lon_0=-61.5 +k=0.9999 +x_0=304800 +y_0=0 +ellps=GRS80 +units=m
                                   +no_defs"))
VHF2018 <- as.data.frame(VHF2018)
VHF2018$UTMZone <- as.integer(22)


coordinates(VHF2019) <- c("Easting", "Northing")
proj4string(VHF2019) <- CRS("+proj=longlat +datum=WGS84")

VHF2019 <- spTransform(VHF2019, CRS("+proj=tmerc +lat_0=0 +lon_0=-61.5 +k=0.9999 +x_0=304800 +y_0=0 +ellps=GRS80 +units=m
                                    +no_defs"))
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
Data17 <- dplyr::select(Data17, -c(Notes, EarTag, Wind, Rain, SampleTimeCat, SampleTimeGeneral, Clouds, Temp, Line, Alive))
Data18 <- dplyr::select(Data18, -c(Notes, EarTag, Wind, Rain, SampleTimeCat, SampleTimeGeneral, Clouds, Temp, Line, Alive))
Data19 <- dplyr::select(Data19, -c(FixLocation,TimeCategory, TempC, Wind, Rain))

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
vhfData$Time_O <- times(vhfData$Time_O)
# add a column that has date and time for razimuth package
vhfData <- add_column(vhfData, datetime = as.POSIXct(paste(vhfData$Date, vhfData$Time_O), format="%Y-%m-%d %H:%M:%S", tz = "America/St_Johns", usetz = TRUE)) 
head(vhfData)
# select only those columns necessary for ATM calculation
vhfData <- vhfData %>%
  dplyr::select(c(-Point,-UTMZone,-Date,-Time_O,-Time,-AM_PM,-Observer,-Grid,-Year))
# rename columns to appropriate headers for ATM calculation
vhfData <- vhfData %>%
  rename(indiv = Frequency, obs_id = GID, azimuth = Azimuth, date = datetime, utm_x = Easting, utm_y = Northing)
# add column of prior_r. This is the maximum distance (m) between observer and transmitter prior to taking the 
# azimuth. If this was not obtained during data collection, still need to have an upper bound for each azimuth 
# BL grid is 500 m wide and buns do not generally move more than ~250m off grid. So max distance away would be 
# approximately 750m. Setting this for every row.
vhfData <- vhfData %>%
  add_column(prior_r = 750)

# clean up environment
rm(Data17, Data18, Data19, liveData, TmpData, VHF2017, VHF2018, VHF2019)

# --------------------------------------- #
#         Individual Collar               #
# --------------------------------------- #
# AR helped with this code - github.com/robitalec
# investigating sample data and functionality of package

# Check example data structure
grouse <- sim_atd(n_loc = 100, sq_km = 4, n_azimuth = 3, dist_vec = c(200,300),
                  kappa_sim = 50, prior_r = 1000, plot = FALSE)$sim_df

# For one individual
oneid <- '149.673'
onedf <- vhfData[vhfData$indiv=='149.673', ]
oneatm <- convert_atm(df = vhfData[vhfData$indiv == oneid,])
visualize_atm(atm_df = oneatm, obs_id = 1, add_prior = T)
atm_fit <- atm_mcmc(atm_df = oneatm, n_mcmc = 500, n_burn = 100)
which_pid(oneatm, obs_id=1)

# visualize one with error 
pdf('graphics/final_triangulation.pdf')
for(i in 1:length(unique(onedf$obs_id))){
  id_tmp <- which_pid(atm_df = oneatm, obs_id = unique(onedf$obs_id)[i])
  plot.new()
  visualize_atm(atm_df = oneatm, obs_id = unique(onedf$obs_id)[i], add_prior = T)
  p_isopleth(df = atm_fit$mu_ls[[id_tmp]]$pdraws, prob_lvls = 0.9, range_extend = 0,
            kde_n=50, col_vec = c(4,4))
  points(matrix(atm_fit$pmode[id_tmp, 2:3], ncol = 2), pch = 21, bg = 4)
  legend("topleft", c("Posterior Mode"), pch = 21, pt.bg = 4, bty = "n")

}
dev.off()

# --------------------------------------- #
#           Multiple Collars              #
# --------------------------------------- #
# AR helped with this code - github.com/robitalec
# For all individuals 
# Unique individuals
uids <- unique(vhfData$indiv)

## Convert all to atm
# Loop through 'em
lsatms <- lapply(uids, function(id) convert_atm(df = vhfData[vhfData$indiv == id,]))

# Name the list for downstream
names(lsatms) <- uids

## Visualize all
pdf('graphics/init_triangulation.pdf')
lapply(seq_along(lsatms), function(x) {
  plot.new()
  mtext(names(lsatms)[[x]])
  lapply(lsatms[[x]]$pid, function(i) visualize_atm(lsatms[[x]], obs_id = i, add_prior = TRUE))
})
dev.off()

# Model all
lsfits <- lapply(lsatms, atm_mcmc, n_mcmc = 10000, n_burn = 1000)
lsfits <- saveRDS(lsfits, "large/lsfits.rds")
lsfits <- readRDS("large/lsfits.rds")

# Visualize all with error
pdf('graphics/final_triangulation.pdf')
lapply(seq_along(lsatms), function(x) {
  print(plot.new())
  mtext(names(lsatms)[[x]])
  lapply(lsatms[[x]]$pid, plyr::failwith(NULL, function(i) {
    visualize_atm(lsatms[[x]], obs_id = i, add_prior = TRUE)
    p_isopleth(df = lsfits[[x]]$mu_ls[[i]]$pdraws, prob_lvls = c(0.95), range_extend = 0,
               kde_n = 50, col_vec = c(4,4))
    points(matrix(lsfits[[x]]$pmode[i, 2:3], ncol = 2), pch = 21, bg = 4)
    legend("topleft", c("Posterior Mode"), pch = 21, pt.bg = 4, bty = "n")
  }
  ))
})

dev.off()

# check convergence ----
pdf('graphics/traceplot.pdf')
lapply(seq_along(lsfits), function(x) {
  plot.new()
  mtext(names(lsfits)[[x]])
  plot_kappa(atm_obj = lsfits[[x]]$kappa_ls, item = "traceplot")
})
dev.off()

pdf('graphics/runmean.pdf')
lapply(seq_along(lsfits), function(x) {
  plot.new()
  mtext(names(lsfits)[[x]])
  plot_kappa(atm_obj = lsfits[[x]]$kappa_ls, item = "run_mean")
})
dev.off()


# --------------------------------------- #
#              Extract Data               #
# --------------------------------------- #
## Extract variances to use in ctmm
lsvars <- lapply(seq_along(lsfits), function(x) {
  rbindlist(lapply(lsfits[[x]]$mu_ls, function(y) {
    xy <- y[['pdraws']]
    data.table(COV.x.y = var(x=xy[,1], y=xy[,2]), pid = y[['pid']], 
               COV.x.x = var(xy[, 1]), COV.y.y = var(xy[, 2]))
  }))[, id := names(lsfits)[[x]]]
})
vars <- rbindlist(lsvars)
# create a dataframe with relocations, IDs, and variances 
# extract relocation data
lsreloc <- lapply(seq_along(lsfits), function(x) {
  rbindlist(lapply(lsfits[[x]]$mu_ls, function(y) {
    xy <- as.matrix(y[['pmode']])
    data.table(utm_x = xy[1,], utm_y=xy[2,], pid = y[['pid']])
  }))[, id := names(lsfits)[[x]]]
})
relocs <- rbindlist(lsreloc)
# join vars and relocs based on date
# add new column of collar, pid to be joined by 
vars <- vars %>% tidyr::unite(cpid, c(id,pid), sep=',', remove=F)
relocs <- relocs %>% tidyr::unite(cpid, c(id,pid), sep=',', remove=F)
df <- full_join(vars, relocs, by='cpid', keep=FALSE)
# now need to get date-time stamp for each reloc for MoveBank 
# get the first date-time for each indiv/obs_id from vhfData
dates <- vhfData %>% 
  dplyr::group_by(indiv, obs_id) %>%
  dplyr::filter(row_number()==1) %>%
  tidyr::unite(cpid, c(indiv, obs_id), sep=',', remove=F) %>%
  dplyr::select(c(-utm_x,-utm_y,-azimuth,-prior_r,-obs_id,-indiv))
df <- left_join(df,dates,by='cpid', keep=FALSE)
# manually went through and removed problematic points from the dataset using Excel 
# input/HomeRanges.txt file describing which points are removed
# remove the two collars that are not in Bloomfield
df <- subset(df, df$indiv != "149.374" &
                 df$indiv != "149.474")
# remove the four individuals that appear in more than one year - keeping the year 
# with the highest number of relocations
df <- subset(df, df$indiv != "149.003" &
                 df$indiv != "149.053" &
                 df$indiv != "149.274" &
                 df$indiv != "149.653" &
                 df$indiv != "149.555")

# reproject into lat/long and save - necessary for MoveBank
coordinates(df) <- c("utm_x", "utm_y")
proj4string(df) <- CRS("+proj=tmerc +lat_0=0 +lon_0=-61.5 +k=0.9999 +x_0=304800 +y_0=0 +ellps=GRS80 +units=m
                                    +no_defs")
df <- spTransform(df, CRS("+init=epsg:4326"))
df <- as.data.frame(df)
df <- df %>%
  rename(lat=utm_y, long=utm_x)

write.csv(df, "output/harestriangulated_razimuth.csv")

# calculate mean time it takes to triangulate the hares 
vhf <- vhfData %>%
  group_by(indiv, GID) %>%
  arrange(datetime, .by_group=TRUE) %>%
  dplyr::summarise(First=first(datetime),
                   Last=last(datetime),
                   difference=difftime(last(datetime), first(datetime), unit="hours"),
                   .groups="keep")

vhf <- drop_na(vhf)
mean(vhf$difference)