# Author: Isabella Richmond (code and data shared between Matteo Rizzuto: github.com/matteorizzuto)
# Last Edited: June 17, 2020

# This script is for the analysis of...


# load required packages 
easypackages::packages("tidyverse", "lme4", "glmmTMB")

# --------------------------------------- #
#             Data Preparation            #
# --------------------------------------- #
ranges <- read_csv("output/RangeStoichRisk.csv")
head(ranges)
# melt data so that it can be compared in t-test form 
rangesCN <- subset(ranges, select = c(CollarFrequency, MeanCN_Core, MeanCN_Home))
rangesCN <- reshape2::melt(rangesCN, id.vars=c("CollarFrequency"))
ggplot(rangesCN, aes(variable, value))+
  geom_boxplot()+
  stat_compare_means()

rangesCP <- subset(ranges, select = c(CollarFrequency, MeanCP_Core, MeanCP_Home))
rangesCP <- reshape2::melt(rangesCP, id.vars=c("CollarFrequency"))
ggplot(rangesCP, aes(variable, value))+
  geom_boxplot()+
  stat_compare_means()

rangesover <- subset(ranges, select = c(CollarFrequency, MeanOverPCA_Core, MeanOverPCA_Home))
rangesover <- reshape2::melt(rangesover, id.vars=c("CollarFrequency"))
ggplot(rangesover, aes(variable, value))+
  geom_boxplot()+
  stat_compare_means()

rangesunder <- subset(ranges, select = c(CollarFrequency, MeanUnderPCA_Core, MeanUnderPCA_Home))
rangesunder <- reshape2::melt(rangesunder, id.vars=c("CollarFrequency"))
ggplot(rangesunder, aes(variable, value))+
  geom_boxplot()+
  stat_compare_means()



# standardize the explanatory variables
cscchc_stack <- cscchc_stack %>%
  add_column(kUD_s = scale(cscchc_stack$kUD, center = TRUE, scale = TRUE)) %>%
  add_column(VAAN_CN_s = scale(cscchc_stack$VAAN_CN, center = TRUE, scale = TRUE)) %>%
  add_column(CoverValue_s = scale(cscchc_stack$CoverValue, center = TRUE, scale = TRUE)) %>%
  add_column(meanhc_s = scale(cscchc_stack$meanhc, center = TRUE, scale = TRUE))

# --------------------------------------- #
#        Linear Mixed Effects Models      #
# --------------------------------------- #
# general linear model with no random effects
mod1 <- glm(ratio ~ MeanPC1 + MeanPC2 + MeanCN + MeanCP, data = rangeuse)
plot(mod1)
qqnorm(residuals(mod1))
qqline(residuals(mod1))
summary(mod1)
# model diagnostics look okay 
# want to model with individual as a random effect 

# linear mixed effect model with random variable
mod2 <- lmer(ratio ~ MeanPC1 + MeanPC2 + MeanCN + MeanCP + (1|CollarFrequency), data = rangeuse)
# problem is that you can't include individual as a random effect when there is one measurement per individual 


# results do not change without outliers. Biologically significant so leaving them in.
# do an AIC to compare canopy cover to horizontal complexity
overlap1 <- glm(overlap ~ CoverValue + meanhc, family = Gamma, data = overlapdata)
overlap2 <- glm(overlap ~ CoverValue, family= Gamma, data = overlapdata)
overlap3 <- glm(overlap ~ meanhc, family = Gamma, data = overlapdata)
overlap4 <- glm(overlap ~ 1, family = Gamma, data = overlapdata)
overlapmodels <- list(overlap1,overlap2,overlap3,overlap4)
overlapmodels <- aictab(cand.set = overlapmodels)
print(overlapmodels)
# intercept is top model in this AIC analysis 

