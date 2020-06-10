# This script is for creating an index of habitat complexity/predation risk for snowshoe hares
# in Newfoundland, Canada. All habitat data was collected in July 2019 (after green-up)
# Refer to KernelEstimation.R script for the calculation of kUD and home range areas 
# that are used in this script

# Author: Isabella Richmond
# Last Edited: June 10, 2020


#------------------------------#
#      Horizontal Complexity   #
#------------------------------#
# NOTE: before this analysis, the horizontal complexity data was calculated as per 
# code below and populated predrisk.csv
# horizontal complexity data was measured using a Nudds board (Nudds, 1977)
# to choose what distance we use to measure horizontal complexity, we need to determine
# which distance has the greatest variation in scores 
# hc refers to raw dataset with only horizontal complexity
hcsum <- hc %>%
  group_by(Distance, Score) %>%
  tally()

ggplot(hcsum, aes(x=Distance,y=n, fill=Score))+
  geom_col()+
  scale_fill_brewer(palette = 1)+
  theme(panel.background = element_rect(fill = "darkgrey"))+
  labs(fill="Score")
# inspection of the graph shows that a distance of 10 metres experiences the most
# variation in horizontal complexity scores
# therefore, use values from 10 m only going forward 
hc <- hc %>% filter(Distance == 10)
# terrestrial predators are all under 1 m tall (coyote & lynx), subset these heights
# from the dataset and average for each plot 
hc <- hc %>% filter(Height %in% c(0.5,1))
# calculate an average horizontal complexity score for each plot
# to do this we have to change the Score variable to an integer
hc <- hc %>% mutate(score = as.numeric(Score))
hcmean <- hc %>%
  group_by(Plot) %>%
  dplyr::summarise(meanhc = mean(score))
cchc <- inner_join(cc, hcmean, by = "Plot")

#------------------------------#
#        Data Preparation      #
#------------------------------#

# load required packages 
library(easypackages)
devtools::install_github("vqv/ggbiplot")
easypackages::packages("tidyverse", "vegan", "sf", "ggbiplot")

# load required data
# habitat complexity/predation risk data
predrisk <- read_csv("input/HC_CleanData_2019.csv")


#------------------------------#
#          Ordination          #
#------------------------------#

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
