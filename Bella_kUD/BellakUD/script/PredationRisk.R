# Author: Isabella Richmond (code and data shared between Matteo Rizzuto: github.com/matteorizzuto)
# Last Edited: April 22, 2020

# This script is for estimating the predation risk of snowshoe hare habitat using structural
# complexity of the environment



# load required packages 
easypackages::packages("tidyr")


# --------------------------------------- #
#             Data Preparation            #
# --------------------------------------- #

# import telemetry datasetss
VHF2019 <- read.csv("input/TelemetryPoints_VHF_2019.csv")
View(VHF2019)