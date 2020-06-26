# This script is for estimating the space use using autocorrelated kernels for snowshoe hares 
# over three years (2016-2019). Relocations only taken in summer season.

# Author: Isabella Richmond (code and data shared between Matteo Rizzuto: github.com/matteorizzuto)
# Last Edited: June 25, 2020
# using default autocorrelated Gaussian reference bandwidth with debiased area (Fleming & Calabrese, 2017)
# much of this code was written by Amanda Droghini - cite in paper
# Droghini, A. 2020. Southwest Alaska moose. Git Repository. Available: https://github.com/accs-uaa/southwest-alaska-moose


devtools::install_github("ctmm-initiative/ctmm")
easypackages::packages("tidyverse", "adehabitatHR", "ctmm")

# -------------------------------------#
#           Data Preparation           #
# -------------------------------------#
# upload the data from MoveBank
hares <- read.csv("input/HaresTriangulatedMoveBank.csv")
head(hares)
# convert to a list of telemetry objects
hares.telem <- as.telemetry(hares)

# -------------------------------------#
#           Error Calibration          #
# -------------------------------------#
# upload the sigloc dataset with error ellipses
sighares <- read.csv("output/harestriangulated.csv")
sighares <- sighares %>% rename(location.lat = coords.x2) %>%
  rename(location.long = coords.x1)
# join based on coordinates because hares has one less row than sighares 
hares <- inner_join(hares,sighares, by=c("location.lat", 'location.long'))
# rename ellipses to appropriate column headers for as.telemetry 
# https://groups.google.com/forum/?oldui=1#!searchin/ctmm-user/VHF$20error$20ellipse%7Csort:date/ctmm-user/TrLypiFqG7E/fuiHCi7eFgAJ
hares <- hares %>%
  rename(COV.x.x = Var_X) %>%
  rename(COV.y.y = Var_Y) %>%
  rename(COV.x.y = Cov_XY) # these are the error ellipses values
# remove columns inappropriate for as.telemetry()
hares <- dplyr::select(hares, -c(datetime, Frequency, Date, Time, X, BadPoint, AngleDiff))
hares.telem <- as.telemetry(hares)

names(hares.telem[[1]])
UERE <- uere.fit(hares.telem[1:29]) # only using calibration data
summary(UERE)

uere(hares.telem) <- UERE
names(hares.telem[[3]]) # now the data are calibrated, as VAR is present"

# -------------------------------------#
#           Removing Outliers          #
# -------------------------------------#
# investigate the data 
plot(hares.telem, col=rainbow(length(hares.telem)))
# use ctmm::outlie to investigate outliers 
# blue indicates speed outlier and red indicates distance
ids <- names(hares.telem)
for (i in 1:length(ids)){
  ctmm::outlie(hares.telem[[i]], plot=TRUE, main=ids[i])
  plotName <- paste("outliers",ids[i],sep="")
  filePath <- paste("output/Outliers/",plotName,sep="")
  finalName <- paste(filePath,"png",sep=".")
  dev.copy(png,finalName)
  dev.off()
  rm(plotName,filePath,finalName)
}
# looks like generally high speeds are related to far distances
# less worried about speed because there is a minimum of 15 hours
# between relocs (usually more than that)
# going to look into outliers for each individual 
# using plotOutlier and looking at ctmm::outlie output simultaneously
source("script/function-plotOutliers.R")
# 149.093
subsetOutlier <- subset(hares, tag.local.identifier == "149.093")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # not removing any
# 149.124
subsetOutlier <- subset(hares, tag.local.identifier == "149.124")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # not removing any 
# 149.173
subsetOutlier <- subset(hares, tag.local.identifier == "149.173")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # removing 2 points past -53.980
subsetOutlier <- filter(subsetOutlier, location.long > (-53.980))
plotOutliers(subsetOutlier,1,nrow(subsetOutlier)) # looks better
# 149.213
subsetOutlier <- subset(hares, tag.local.identifier == "149.213")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # removing 2 points past -53.9930
subsetOutlier <- filter(subsetOutlier, location.long > (-53.9930))
plotOutliers(subsetOutlier,1,nrow(subsetOutlier)) # looks better
# 149.233
subsetOutlier <- subset(hares, tag.local.identifier == "149.233")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # removing 2 points past -53.978
subsetOutlier <- filter(subsetOutlier, location.long < (-53.978))
plotOutliers(subsetOutlier,1,nrow(subsetOutlier)) # looks better
# 149.294
subsetOutlier <- subset(hares, tag.local.identifier == "149.294")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # removing 1 point past -53.9800
subsetOutlier <- filter(subsetOutlier, location.long < (-53.9800))
plotOutliers(subsetOutlier,1,nrow(subsetOutlier)) # looks better\
# 149.394
subsetOutlier <- subset(hares, tag.local.identifier == "149.394")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # removing 1 point past 48.354
subsetOutlier <- filter(subsetOutlier, location.lat < (48.354))
plotOutliers(subsetOutlier,1,nrow(subsetOutlier)) # looks better
# 149.423
subsetOutlier <- subset(hares, tag.local.identifier == "149.423")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # removing 1 point past 48.354
subsetOutlier <- filter(subsetOutlier, location.lat < (48.354))
plotOutliers(subsetOutlier,1,nrow(subsetOutlier)) # looks better
# 149.452
subsetOutlier <- subset(hares, tag.local.identifier == "149.452")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # not removing any 
# 149.513
subsetOutlier <- subset(hares, tag.local.identifier == "149.513")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # removing 1 point below 48.350
subsetOutlier <- filter(subsetOutlier, location.lat > (48.350))
plotOutliers(subsetOutlier,1,nrow(subsetOutlier)) # looks better
# 149.535
subsetOutlier <- subset(hares, tag.local.identifier == "149.535")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # removing 1 point past -53.986
subsetOutlier <- filter(subsetOutlier, location.long > (-53.986))
plotOutliers(subsetOutlier,1,nrow(subsetOutlier)) # looks better
# 149.594
subsetOutlier <- subset(hares, tag.local.identifier == "149.594")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # not removing any 
# 149.613
subsetOutlier <- subset(hares, tag.local.identifier == "149.613")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # removing 1 point past -53.9750
subsetOutlier <- filter(subsetOutlier, location.long < (-53.9750))
plotOutliers(subsetOutlier,1,nrow(subsetOutlier)) # looks better
# 149.633
subsetOutlier <- subset(hares, tag.local.identifier == "149.633")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # not remmoving any 
# 149.673
subsetOutlier <- subset(hares, tag.local.identifier == "149.673")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # removing 1 point below -53.9900
subsetOutlier <- filter(subsetOutlier, location.long > (-53.9900))
plotOutliers(subsetOutlier,1,nrow(subsetOutlier)) # looks better
# 150.032
subsetOutlier <- subset(hares, tag.local.identifier == "150.032")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # not removing any 
# 150.052
subsetOutlier <- subset(hares, tag.local.identifier == "150.052")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # removing 1 point past -53.9900 and 1 past 48.349
subsetOutlier <- filter(subsetOutlier, location.long > (-53.9900) & location.lat > (48.349))
plotOutliers(subsetOutlier,1,nrow(subsetOutlier)) # looks better
# 150.072
subsetOutlier <- subset(hares, tag.local.identifier == "150.072")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # removing 1 point past 48.350
subsetOutlier <- filter(subsetOutlier, location.lat > (48.350))
plotOutliers(subsetOutlier,1,nrow(subsetOutlier)) # looks better
# 150.091
subsetOutlier <- subset(hares, tag.local.identifier == "150.091")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # not removing any 
# 150.111
subsetOutlier <- subset(hares, tag.local.identifier == "150.111")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # removing 1 point past -53.976
subsetOutlier <- filter(subsetOutlier, location.long < (-53.976))
plotOutliers(subsetOutlier,1,nrow(subsetOutlier)) # looks better
# 150.154
subsetOutlier <- subset(hares, tag.local.identifier == "150.154")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # not removing any
# 150.173
subsetOutlier <- subset(hares, tag.local.identifier == "150.173")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # not removing any 
# 150.191
subsetOutlier <- subset(hares, tag.local.identifier == "150.191")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # not removing any 
# 150.232
subsetOutlier <- subset(hares, tag.local.identifier == "150.232")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # removing 1 point past 48.353
subsetOutlier <- filter(subsetOutlier, location.lat < (48.353))
plotOutliers(subsetOutlier,1,nrow(subsetOutlier)) # looks better
# 150.273
subsetOutlier <- subset(hares, tag.local.identifier == "150.273")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # not removing any 
# 150.314
subsetOutlier <- subset(hares, tag.local.identifier == "150.314")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # removing 1 point past -53.9775
subsetOutlier <- filter(subsetOutlier, location.long < (-53.9775))
plotOutliers(subsetOutlier,1,nrow(subsetOutlier)) # looks better
# 150.332
subsetOutlier <- subset(hares, tag.local.identifier == "150.332")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # not removing any 
# 150.373
subsetOutlier <- subset(hares, tag.local.identifier == "150.373")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # not removing any 
# 150.392
subsetOutlier <- subset(hares, tag.local.identifier == "150.392")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # not removing any

# create a dataset without outliers to run through the analysis 
haresClean <- hares %>%
  dplyr::filter(!(tag.local.identifier=="149.173" & location.long < (-53.980)|
                  tag.local.identifier=="149.213" & location.long < (-53.9930)|
                  tag.local.identifier=="149.233" & location.long > (-53.978)|
                  tag.local.identifier=="149.294" & location.long > (-53.9800)|
                  tag.local.identifier=="149.394" & location.lat > (48.354)|
                  tag.local.identifier=="149.423" & location.lat > (48.354)|
                  tag.local.identifier=="149.513" & location.lat < (48.350)|
                  tag.local.identifier=="149.535" & location.long < (-53.986)|
                  tag.local.identifier=="149.613" & location.long > (-53.9750)|
                  tag.local.identifier=="149.673" & location.long < (-53.9900)|
                  tag.local.identifier=="150.052" & (location.long < (-53.9900) | location.lat < (48.349))|
                  tag.local.identifier=="150.072" & location.lat < (48.350)|
                  tag.local.identifier=="150.111" & location.long > (-53.976)|
                  tag.local.identifier=="150.232" & location.lat > (48.353) |
                  tag.local.identifier=="150.314" & location.long > (-53.9775)
                  ))
# 19 total outliers removed 
# save this as the cleaned version 
write.csv(haresClean, "output/haresClean.csv")
# convert to telemetry object to use going forward
hares.telem.clean <- as.telemetry(haresClean)


# looking at the variograms to explore space use patterns
# load varioPlot function from function-plotVariograms.R 
# zoom is false because there is only one measurement taken per day (max)
source("script/function-plotVariograms.R")
varioPlot(hares.telem.clean,filePath="output/Variograms/",zoom = FALSE)
# Most of the variograms don't fit well

# ----------------------------------- #
#        Select Model Parameters      #
# ----------------------------------- #
# loop through individuals and guess the parameters for each 
hares.guess.initial <- lapply(hares.telem.clean[1:length(hares.telem.clean)], 
                              function(b) ctmm.guess(b,CTMM=ctmm(error=TRUE),
                                                     interactive=FALSE) )

# then use the guessed parameter values in ctmm.fit for each individual
# Using initial guess parameters and ctmm.select
# ctmm.select will rank models and the top model can be chosen to generate an aKDE
# chose pHREML due to the small sample size (Fleming et al., 2019)
hares.fit <- lapply(1:length(hares.telem.clean), 
                    function(i) ctmm.select(data=hares.telem.clean[[i]],
                                            CTMM=hares.guess.initial[[i]],
                                            verbose=TRUE,trace=TRUE, cores=0,
                                            method = "pHREML") )
# Add seasonal animal ID names to fitModels list
names(hares.fit) <- names(hares.telem)
# The warning "pREML failure: indefinite ML Hessian" is normal if some autocorrelation parameters cannot be well resolved.

# Place model selection parameters for all individuals in dataframe
hares.models.summary <- lapply(hares.fit,function(x) summary(x))
hares.models.summary <- plyr::ldply(hares.models.summary, rbind)

# Place model name in df
modelRows <- lapply(hares.fit,function(x) row.names(summary(x)))
modelRows <- plyr::ldply(modelRows, rbind)
modelRows <- modelRows %>% 
  pivot_longer(cols = -.id,
               values_to="model",names_to="rank",
               values_drop_na = TRUE)

modelSummary <- cbind(modelRows,hares.models.summary)
# Delete duplicate id column. Join doesn't work because .id is not a unique key
modelSummary <- modelSummary[,-4]
# Subset only the highest ranked models
topModels <- distinct(modelSummary,.id, .keep_all=TRUE) 
names(modelSummary) <- enc2utf8(names(modelSummary))
write_csv(modelSummary,"output/haresmodelsummary.csv")
write_csv(topModels,"output/harestopmodels.csv")


# ---------------------------------- #
#        Reassess Variograms         #
# ---------------------------------- #
filePath <- paste("output/Variograms/")

lapply(1:length(hares.telem), 
       function (a) {
         plotName <- paste(names(hares.fit[a]),sep="_")
         plotPath <- paste(filePath,plotName,sep="")
         finalName <- paste(plotPath,"png",sep=".")
         
         plot(ctmm::variogram(hares.telem[[a]],CI="Gauss"),
              CTMM=hares.fit[[a]][1:2],
              col.CTMM=c("red","blue","purple","green"),
              fraction=1,
              level=c(0.5,0.95),
              main=names(hares.fit[a]))
         
         dev.copy(png,finalName)
         dev.off()
         
       }
)

# ---------------------------------- #
#           Create aKDE              #
# ---------------------------------- #
