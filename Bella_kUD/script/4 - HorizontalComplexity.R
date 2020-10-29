# This script is used to calculate the mean horizontal complexity at each sampling point 
# as per Nudds, 1977. Mean horizontal complexity will then be used to assess the structural 
# complexity of different habitats, in conjunction with other measured variables

# Author: Isabella Richmond
# Last Edited: August 15, 2020

# load required packages 
library(easypackages)
easypackages::packages("tidyverse")

#------------------------------#
#      Horizontal Complexity   #
#------------------------------#
# NOTE: before this analysis, the horizontal complexity data was calculated as per 
# code below and populated predrisk.csv
# horizontal complexity data was measured using a Nudds board (Nudds, 1977)
# to choose what distance we use to measure horizontal complexity, we need to determine
# which distance has the greatest variation in scores 
# hc refers to raw dataset with only horizontal complexity
hc <- read_csv("input/CSRawData_HorizontalComplexity.csv")
hcsum <- hc %>%
  group_by(Distance, Score) %>%
  tally()
hcsum <- hcsum %>% 
  mutate(Score = as.factor(Score))
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

write_csv(hcmean, "output/meanhorizontalcomplexity.csv")
# this data was then taken and used to populate horizontal complexity .csv spreadsheet
# import fully cleaned spreadsheet from input folder as HC_CleanData_2019.csv