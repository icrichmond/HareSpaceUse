# This script is for estimating the space use of snowshoe hare collared in Bloomfield, NL 
# over three years (2016-2019). Relocations only taken in summer season.
# Author: Isabella Richmond
# Last Edited: February 25, 2020

# load required packages 
easypackages::packages("sp", "maptools", "tmap", "tmaptools", 
              "adehabitatHR", "adehabitatHS", "adehabitatLT", "ellipse", 
              "nleqslv", "adehabitatMA", "adehabitatHR","dplyr", "gdtools", "ggmap", "ggplot2", 
              "ggrepel", "ggsci", "ggthemes", "maps", "raster", "spatial", 
              "tidyr", "readr","rgdal", "rgeos", "reshape2", "dismo")

### NOTE: one package (sigloc) is no longer maintained - need to manaully load package ###

# import telemetry dataset
VHF <- read_csv("input/TelemetryPoints_VHF_2019.csv")
View(VHF)

# following section developed from Matteo Rizzuto's code
# process the dataset so it is ready for home range estimation.