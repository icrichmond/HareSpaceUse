# This script is for creating an index of habitat complexity/predation risk for snowshoe hares
# in Newfoundland, Canada. All habitat data was collected in July 2019 (after green-up)
# Refer to KernelEstimation.R script for the calculation of kUD and home range areas 
# that are used in this script

# Author: Isabella Richmond
# Last Edited: June 10, 2020

# load required packages 
library(easypackages)
devtools::install_github("vqv/ggbiplot")
easypackages::packages("tidyverse", "vegan", "sf", "ggbiplot")

# load required data
# home range areas 
homeranges <- read_csv("output/homerangeareas.csv")
kernel95 <- read_sf("output/hares.kudhr.95.shp")
kernel50 <- read_sf("output/hares.kudhr.50.shp")
# habitat complexity/predation risk data
predrisk <- read_csv("input/HC_CleanData_2019.csv")
# visualize home ranges 
ggplot(kernel95) +
  geom_sf(fill = "dark grey") +
  geom_sf(data = kernel50, aes(fill = "light grey"))
  xlab("Longitude") + ylab("Latitude")

# prepare home range data - we are using range ratio with 50:95 home range areas 
rangeuse <- pivot_longer(homeranges, cols = 2:ncol(homeranges), names_to = "CollarFrequency", values_to = "HomeRangeArea")
rangeuse <- pivot_wider(rangeuse, names_from = "Kernel", values_from = "HomeRangeArea")
# calculating a 50:95 range use ratio (Webber et al., 2020)
rangeuse <- add_column(rangeuse, ratio = rangeuse$`50`/rangeuse$`95`)
write_csv(rangeuse, "output/rangeuseratio.csv")

# perform ordination on predation risk data to see which variables are most important
# subset data to only numeric variables
# drop any NA values

ord = prcomp(drop_na(predrisk[,c(4:20)]), scale = TRUE, center = TRUE)
summary(ord)
# visualize PCA
ggbiplot(ord, choices=c(1,2))
ggsave("graphics/predriskpca.jpg", height = 200, width = 200, units = "mm")
# add PC1 and PC2 to predrisk dataset
# drop NA from dataset
predriskpca <- drop_na(predrisk)
predriskpca <- predriskpca %>% add_column(PC1 = ord$x[,1]) %>%
  add_column(PC2 = ord$x[,2])
write_csv(predriskpca, "output/predationriskpca.csv")
