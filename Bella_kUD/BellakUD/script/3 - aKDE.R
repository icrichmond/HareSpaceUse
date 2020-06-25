# This script is for estimating the space use using autocorrelated kernels for snowshoe hares 
# over three years (2016-2019). Relocations only taken in summer season.

# Author: Isabella Richmond (code and data shared between Matteo Rizzuto: github.com/matteorizzuto)
# Last Edited: June 25, 2020
# using default autocorrelated Gaussian reference bandwidth with debiased area (Fleming & Calabrese, 2017)
# much of this code was written by Amanda Droghini - cite in paper
# Droghini, A. 2020. Southwest Alaska moose. Git Repository. Available: https://github.com/accs-uaa/southwest-alaska-moose


devtools::install_github("ctmm-initiative/ctmm")
easypackages::packages("tidyverse", "adehabitatHR")

# -------------------------------------#
#           Data Preparation           #
# -------------------------------------#
# upload the data from MoveBank
hares <- read.csv("input/HaresTriangulatedMoveBank.csv")
head(hares)
# convert to a list of telemetry objects
hares.telem <- as.telemetry(hares)

#### Data Visualization & Variograms ####
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
# going to look into outliers for each individual 


# looking at the variograms to explore space use patterns
# load varioPlot function from function-plotVariograms.R 
# zoom is false because there is only one measurement taken per day (max)
varioPlot(hares.telem,filePath="output/Variograms/",zoom = FALSE)
# Most of the variograms don't fit well

# ----------------------------------- #
#        Select Model Parameters      #
# ----------------------------------- #
# loop through individuals and guess the parameters for each 
hares.guess.initial <- lapply(hares.telem[1:length(hares.telem)], 
                              function(b) ctmm.guess(b,CTMM=ctmm(error=TRUE),
                                                     interactive=FALSE) )

# then use the guessed parameter values in ctmm.fit for each individual
# Using initial guess parameters and ctmm.select
# ctmm.select will rank models and the top model can be chosen to generate an aKDE
hares.fit <- lapply(1:length(hares.telem), 
                    function(i) ctmm.select(data=hares.telem[[i]],
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

