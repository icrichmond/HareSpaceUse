# This script is for estimating the space use using autocorrelated kernels for snowshoe hares 
# over three years (2016-2019). Relocations only taken in summer season.

# Author: Isabella Richmond (code and data shared between Matteo Rizzuto: github.com/matteorizzuto)
# Last Edited: June 25, 2020
# using default autocorrelated Gaussian reference bandwidth with debiased area (Fleming & Calabrese, 2017)
# much of this code was written by Amanda Droghini - cite in paper
# Droghini, A. 2020. Southwest Alaska moose. Git Repository. Available: https://github.com/accs-uaa/southwest-alaska-moose

devtools::install_github("ctmm-initiative/ctmm")
easypackages::packages("tidyverse", "adehabitatHR", "ctmm", "tidyverse")

# -------------------------------------#
#           Data Preparation           #
# -------------------------------------#
# upload the data from MoveBank
hares <- read.csv("input/HaresTriangulatedMoveBank_Razimuth.csv")
head(hares)
# convert to a list of telemetry objects
hares.telem <- as.telemetry(hares)

# -------------------------------------#
#           Error Calibration          #
# -------------------------------------#
# upload the razimuth dataset with error ellipses
razhares <- read.csv("output/harestriangulated_razimuth.csv")
razhares <- razhares %>% rename(location.lat = lat) %>%
  rename(location.long = long)
# join based on coordinates 
hares <- inner_join(hares,razhares, by=c("location.lat", 'location.long'))
# remove columns inappropriate for as.telemetry()
hares <- dplyr::select(hares, -c(date, id.x, id.y, pid.x, pid.y, indiv, obs_id, cpid))
# recalculate the as.telemetry() object with keep = TRUE so that COV.x.y is included 
hares.telem <- as.telemetry(hares, keep = TRUE)
# can see in the plotting that the ellipses are present
par(mfrow=c(1,1))
plot(hares.telem, error = 1, xlim = c(1000,-1000))
plot(hares.telem, error = 0)

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
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # not removing any 
# 149.294
subsetOutlier <- subset(hares, tag.local.identifier == "149.294")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # removing 1 point past -53.9800
subsetOutlier <- filter(subsetOutlier, location.long < (-53.9800))
plotOutliers(subsetOutlier,1,nrow(subsetOutlier)) # looks better
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
# 149.555
subsetOutlier <- subset(hares, tag.local.identifier == "149.555")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # not removing any 
# 149.594
subsetOutlier <- subset(hares, tag.local.identifier == "149.594")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # not removing any 
# 149.613
subsetOutlier <- subset(hares, tag.local.identifier == "149.613")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # removing 1 point past -53.9750
subsetOutlier <- filter(subsetOutlier, location.long < (-53.9775))
plotOutliers(subsetOutlier,1,nrow(subsetOutlier)) # looks better
# 149.633
subsetOutlier <- subset(hares, tag.local.identifier == "149.633")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # not removing any 
# 149.673
subsetOutlier <- subset(hares, tag.local.identifier == "149.673")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # not removing any
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
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # not removing any
# 150.132
subsetOutlier <- subset(hares, tag.local.identifier == "150.132")
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # not removing any 
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
plotOutliers(subsetOutlier, 1, nrow(subsetOutlier)) # not removing any
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
                    tag.local.identifier=="149.294" & location.long > (-53.9800)|
                    tag.local.identifier=="149.394" & location.lat > (48.354)|
                    tag.local.identifier=="149.423" & location.lat > (48.354)|
                    tag.local.identifier=="149.513" & location.lat < (48.350)|
                    tag.local.identifier=="149.535" & location.long < (-53.986)|
                    tag.local.identifier=="149.613" & location.long > (-53.9775)|
                    tag.local.identifier=="150.052" & (location.long < (-53.9900) | location.lat < (48.349))|
                    tag.local.identifier=="150.072" & location.lat < (48.350)|
                    tag.local.identifier=="150.314" & location.long > (-53.9775)
  ))
# 16 total outliers removed 
# save this as the cleaned version 
write.csv(haresClean, "output/haresClean.csv")

# convert to telemetry object to use going forward
hares.telem.clean <- as.telemetry(haresClean, keep = TRUE)


# looking at the variograms to explore space use patterns
# load varioPlot function from function-plotVariograms.R 
# zoom is false because there is only one measurement taken per day (max)
source("script/function-plotVariograms.R")
varioPlot(hares.telem.clean,filePath="output/Variograms/",zoom = FALSE)
# variograms look ok - could be better

# ----------------------------------- #
#        Select Model Parameters      #
# ----------------------------------- #
# we are selecting model parameters and INCLUDING error ellipses
# including error ellipses is really important for species that aren't highly 
# mobile such as snowshoe hare 
# error ellipses were extracted from the razimuth package
# they are not perfect as the razimuth package does not provide Gaussion-distributed
# ellipses like sigloc, but rather irregular shapes. However, most error shapes for this
# dataset were approximately circular-elliptical, so error variances should be okay 
# sigloc does not calculate errors correctly (Gerber et al., 2018), which is why we switched to razimuth
# keep telemetry error with keep=TRUE for COV.x.x, COV.y.y, and COV.x.y
# loop through individuals and guess the parameters for each 
hares.guess.initial.e <- lapply(hares.telem.clean[1:length(hares.telem.clean)], 
                                function(b) ctmm.guess(b,CTMM=ctmm(error=TRUE),
                                                       interactive=FALSE) )

# then use the guessed parameter values in ctmm.fit for each individual
# Using initial guess parameters and ctmm.select
# ctmm.select will rank models and the top model can be chosen to generate an aKDE
# chose pHREML due to the small sample size (Fleming et al., 2019)
hares.fit.e <- lapply(1:length(hares.telem.clean), 
                      function(i) ctmm.select(data=hares.telem.clean[[i]],
                                              CTMM=hares.guess.initial.e[[i]],
                                              verbose=TRUE,trace=TRUE, cores=0,
                                              method = "pHREML") )
# save large output
saveRDS(hares.fit.e, "large/haresfitraz.rds")
hares.fit.e <- readRDS("large/haresfitraz.rds")

# Add seasonal animal ID names to fitModels list
names(hares.fit.e) <- names(hares.telem)
# The warning "pREML failure: indefinite ML Hessian" is normal if some autocorrelation parameters cannot be well resolved.

# Place model selection parameters for all individuals in dataframe
hares.models.summary.e <- lapply(hares.fit.e,function(x) summary(x))
hares.models.summary.e <- plyr::ldply(hares.models.summary.e, rbind)

# Place model name in df
modelRows.e <- lapply(hares.fit.e,function(x) row.names(summary(x)))
modelRows.e <- plyr::ldply(modelRows.e, rbind)
modelRows.e <- modelRows.e %>% 
  pivot_longer(cols = -.id,
               values_to="model",names_to="rank",
               values_drop_na = TRUE)

modelSummary.e <- cbind(modelRows.e,hares.models.summary.e)
# Delete duplicate id column. Join doesn't work because .id is not a unique key
modelSummary.e <- modelSummary.e[,-4]
# Subset only the highest ranked models
topModels.e <- distinct(modelSummary.e,.id, .keep_all=TRUE) 
names(modelSummary.e) <- enc2utf8(names(modelSummary.e))
write_csv(modelSummary.e,"output/haresmodelsummary.e.csv")
write_csv(topModels.e,"output/harestopmodels.e.csv")
# DOF is large enough (over 4-5) for all individuals except 150.132

# try bootstrapping 150.132 to see if we can compensate for small sample size 
h150132 <- hares.fit.e$'150.132'
# top model for 150.132 was IID anisotropic error model 
h150132 <- h150132$`IID anisotropic error`
h150132.telem <- hares.telem.clean$'150.132'
h150132.boot <- ctmm.boot(h150132.telem,h150132,iterate=TRUE,trace=2)
saveRDS(h150132.boot, "large/boot150132.rds")
summary(h150132.boot$MLE)
# area is still not large enough - remove 150.132 from analysis 
hares.telem.clean <- purrr::list_modify(hares.telem.clean, "150.132"=NULL)
hares.fit.e <- purrr::list_modify(hares.fit.e, "150.132" = NULL)

# save final telemetry and fit objects
saveRDS(hares.telem.clean, "large/harestelemclean_final.rds")
saveRDS(hares.fit.e, "large/haresfite_final.rds")


# ---------------------------------- #
#        Reassess Variograms         #
# ---------------------------------- #
# plot variograms with model fit 
filePath <- paste("output/Variograms/Models/")

lapply(1:length(hares.telem.clean), 
       function (a) {
         plotName <- paste(names(hares.fit.e[a]),sep="_")
         plotPath <- paste(filePath,plotName,sep="")
         finalName <- paste(plotPath,"png",sep=".")
         
         plot(ctmm::variogram(hares.telem.clean[[a]],CI="Gauss"),
              CTMM=hares.fit.e[[a]][1:2],
              col.CTMM=c("red","blue","purple","green"),
              fraction=1,
              level=c(0.5,0.95),
              main=names(hares.fit.e[a]))
         
         dev.copy(png,finalName)
         dev.off()
         
       }
)


# ---------------------------------- #
#       Select Final Models          #
# ---------------------------------- #
# Select only top models from all possible models
finalNames <- names(hares.fit.e)
finalMods <- lapply(1:length(hares.fit.e), 
                    function(i) hares.fit.e[[i]][1][[1]]) 
names(finalMods) <- finalNames
# save final models 
saveRDS(finalMods, "large/finalmodelse.rds")
# clean entire environment, need space for aKDE
remove(list=ls())

# ---------------------------------- #
#           Create aKDE              #
# ---------------------------------- #
# need clean telemetry data and final models to generate home ranges using aKDE
hares.telem.clean <- readRDS("large/harestelemclean_final.rds")
hares.finalmods.e <- readRDS("large/finalmodelse.rds")

#### Specify extent ----
# Get extent for each telemetry set
ee <- lapply(hares.telem.clean,function(x) extent(x))
ee <- data.frame(matrix(unlist(ee), 
                        nrow=length(ee), 
                        byrow=T))
colnames(ee) <- c("min.x","max.x","min.y","max.y")

# Find absolute minimum and maximum
# Pad it to prevent home ranges from getting cut off
eeMatrix <- c(min(ee$min.x)-1000,max(ee$max.x)+1000,min(ee$min.y)-1000,max(ee$min.y)+1000)
eeMatrix<-matrix(data=eeMatrix,nrow=2,ncol=2,dimnames=list(c("min","max")))
colnames(eeMatrix)<-c("x","y")
ee <- as.data.frame(eeMatrix)

rm(eeMatrix)

# Split data into three chunks ----
# Need to do this to avoid memory limit error
# Work laptop has 16 GB of RAM, which is not enough even after increasing memory limit and trying to run on "fresh" (rebooted) computer

# Order calibratedData and finalMods alphabetically by IDs
ids <- names(hares.finalmods.e)
ids <- ids[order(ids)]

hares.telem.clean <- hares.telem.clean[ids]
hares.finalmods.e <- hares.finalmods.e[ids]

# calculate home ranges
# debias = TRUE debiases the distribution for area estimation (AKDEc, Fleming et al., 2019)
homeRanges <- akde(data=hares.telem.clean, debias=TRUE, CTMM=hares.finalmods.e, grid=ee)
plot(homeRanges$'149.124')

# Export homeRanges
saveRDS(homeRanges, "large/akdehomeranges.rds")

# get home range size 
hr_size <- lapply(1:length(homeRanges),
                  function(i)
                    summary(homeRanges[[i]])$CI)
names(hr_size) <- names(homeRanges)
hr_df <- plyr::ldply(hr_size, data.frame)
hr_df

rm(list=ls())
