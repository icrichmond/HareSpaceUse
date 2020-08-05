# This script is for creating an index of habitat complexity/predation risk for snowshoe hares
# in Newfoundland, Canada. All habitat data was collected in July 2019 (after green-up)
# Refer to KernelEstimation.R script for the calculation of kUD and home range areas 
# that are used in this script

# Author: Isabella Richmond
# Last Edited: June 16, 2020

#------------------------------#
#        Data Preparation      #
#------------------------------#

# load required packages 
library(easypackages)
devtools::install_github("vqv/ggbiplot")
easypackages::packages("tidyverse", "vegan", "sf", "ggbiplot", "data.table")

# load required data
# habitat complexity/predation risk data
predrisk <- read_csv("input/HC_CleanData_2019.csv")
head(predrisk)
# drop na's from dataset so that subsets are same number of rows and there are no 
# na's present for PCA
predrisk <- drop_na(predrisk)
# separate the data into two datasets - one for understorey variables and one for 
# overstorey variables 
under <- subset(predrisk, select = c(Plot, Date, Observers, ShrubCountTotal, FallenLogsTotal,
                                     HorizComplex, LeafLittDepth, AvgShrubHeight, AvgShrubDiam,
                                     StumpCount, AvgFallenDist))
over <- subset(predrisk, select = c(Plot, Date, Observers, CanopyIntOver, CanopyIntUnder,
                                    CanopyIntTotal, CanClos, DeadCount, AvgOverDist,
                                    AvgOverDBH,AvgUnderDist,AvgUnderDBH))
#------------------------------#
#          Ordination          #
#------------------------------#
# perform ordination on all variables to see how informative the axes are 
predord <- prcomp(predrisk[,c(4:20)], scale=TRUE, center=TRUE)
summary(predord)
ggbiplot(predord, choices=c(1,2))
# remove highly correlated variables, CanopyIntOver, AvgOverDBH, AvgUnderDBH, AvgOverDist
newpredrisk <- select(predrisk, c(-AvgOverDBH, -CanopyIntOver, -AvgUnderDBH, -AvgOverDist,-FallenLogsTotal))
newpredord <- prcomp(newpredrisk[,c(4:15)],scale=TRUE,center=TRUE)
summary(newpredord)
ggbiplot(newpredord)
# removing the highly correlated variable does not improve performance, split variables
# into understorey and overstorey

# perform ordination on understorey predation data to see which variables are most 
# important
# subset data to only numeric variables
# drop any NA values
underord = prcomp(under[,c(4:11)], scale = TRUE, center = TRUE)
summary(underord)
# visualize understorey PCA
ggbiplot(underord, choices=c(1,2))
ggsave("graphics/understoreypredriskpca.jpg", height = 200, width = 200, units = "mm")

# perform ordination on overstorey predation data to see which variables are most 
# important 
# subset data to only numeric variables 
# drop any NA values 
overord = prcomp(over[,c(4:12)], scale = TRUE, center = TRUE)
summary(overord)
ggbiplot(overord, choices=c(1,2))
# remove highly correlated variables (same magnitude and direction of vectors)
# Remove AvgOverDBH, CanopyIntUnder, CanopyIntOVer
newover <- select(over, c(-AvgOverDBH, -CanopyIntUnder, -CanopyIntOver))
newoverord = prcomp(newover[,c(4:9)], scale = TRUE, center = TRUE)
summary(newoverord)
ggbiplot(newoverord, choices=c(1,2))
ggsave("graphics/overstoreypredriskpca.jpg", height = 200, width = 200, units = "mm")


# extract PC1 from understorey PCA and overstorey PCA for modelling 
predriskpca <- predrisk %>% add_column(underPCA = underord$x[,1]) %>%
  add_column(overPCA = newoverord$x[,1])
write_csv(predriskpca, "output/predationriskpca.csv")

# regress PCA axes against variables to determine their relationship 
underpca <- under %>% add_column(underPCA = underord$x[,1])
undermelt <- melt(underpca,id.vars = c("Plot","Date","Observers", "underPCA") )
ggplot(data=undermelt, aes(x=underPCA, y = value))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~variable, scales="free_y")
ggsave("graphics/underPCAvariables.png")
# understorey PCA axis has a positive relationship with ShrubCountTotal, HorizComplex,
# AvgFallenDist, StumpCount, LeafLittDepth, and ShrubHeight
# understorey PCA axis has a negative relationship with FallenLogsTotal and AvgShrubDiam
# therefore, understorey PCA axis overall has a positive relationship with predation risk
# and/or habitat complexity, which is what we want 

overpca <- newover %>% add_column(overPCA = newoverord$x[,1])
overmelt <- melt(overpca,id.vars = c("Plot","Date","Observers", "overPCA") )
ggplot(data=overmelt, aes(x=overPCA, y = value))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~variable, scales="free_y")
ggsave("graphics/overPCAvariables.png")
# overstorey PCA axis has a negative relationship with CanopyIntTotal, CanClos,
# AvgUnderDist, and DeadCount
# overstorey PCA axis has a positive relationship with AvgOverDist and AvgUnderDBH
# therefore, overstorey PCA axis overall has a negative relationship with predation risk
# and/or habitat complexity, which means we will multiply it be -1 so that the
# relationship is inverted and it is easier to interpret
predriskpca <- predriskpca %>% mutate(overPCA = (overPCA*-1))
write_csv(predriskpca, "output/predationriskpca.csv")




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