# Author: Isabella Richmond
# Last Edited: July 21, 2020

# This script is for the analysis of...


# load required packages 
easypackages::packages("tidyverse", "lme4", "glmmTMB", "data.table")

# --------------------------------------- #
#             Data Preparation            #
# --------------------------------------- #
# import stoich and kUD values at each complexity sampling point 
stoichkud <- fread("output/cs_stoich_kud.csv")
# load predation risk data
predrisk <- fread("output/predationriskpca.csv")

# combine all variables so there are values for each sampling plot 
full <- stoichkud %>%
  dplyr::select(-V1) %>%
  tibble::add_column(overPCA = predrisk$overPCA)%>%
  tibble::add_column(underPCA = predrisk$underPCA)

# get median value of all collars for each sampling plot 
full <- full %>%
  mutate(kUDmean = rowMeans(dplyr::select(., starts_with("X")))) %>%
  rowwise() %>% 
  mutate(kUDmedian =  median(c(!!! rlang::syms(grep('X', names(.), value=TRUE)))))

# start by standardizing the explanatory variables
full <- full %>%
  add_column(kUDmedian_s = scale(full$kUDmedian, center = TRUE, scale = TRUE)) %>%
  add_column(kUDmean_s= scale(full$kUDmean, center = TRUE, scale = TRUE)) %>%
  add_column(VAAN_CN_s = scale(full$VAAN_CN, center = TRUE, scale = TRUE)) %>%
  add_column(VAAN_CP_s = scale(full$VAAN_CP, center = TRUE, scale = TRUE)) %>%
  add_column(overPCA_s = scale(full$overPCA, center = TRUE, scale = TRUE)) %>%
  add_column(underPCA_s = scale(full$underPCA, center = TRUE, scale = TRUE))

# if analyzing with individual as a random variable, stack the data so that it is 
# ready for analysis. If using median kUD value, keep the dataframe as is and
# skip the next two lines.
full_stack <- pivot_longer(full, cols = starts_with("X"), names_to = "CollarID", names_prefix = "X", values_to = "kUD")
# fix dataset so appropriate columns are factors
full_stack <- full_stack %>% mutate(Plot = as.factor(Plot)) %>%
  mutate(CollarID = as.factor(CollarID))
head(full_stack)

# --------------------------------------- #
#   Generalized Linear Models with kUD    #
# --------------------------------------- #
# now we have plots with PCA values for complexity, 
# kernel utilization values for all 30 hares, 
# median kUD values for every plot and stoich values for 
# lowland blueberry C:N and C:P

# first going to test if there is a relationship between median values and 
# predation risk/food quality
medianmodel <- glm(kUDmedian_s ~ overPCA_s + underPCA_s + VAAN_CN_s + VAAN_CP_s, data = full)
plot(medianmodel)
qqnorm(residuals(medianmodel))
qqline(residuals(medianmodel))
# residuals are extremely non-normal due to skew in the data

meanmodel <- glm(kUDmean_s ~ overPCA_s + underPCA_s + VAAN_CN_s + VAAN_CP_s, data = full)
plot(meanmodel)
qqnorm(residuals(medianmodel))
qqline(residuals(medianmodel))
summary(meanmodel)


# in order to make the model work, nAGQ is set to zero. This 
# indicates less precision with respect to the effects of the 
# random variables
model1 <- lmer(kUD ~ overPCA_s + underPCA_s + VAAN_CN_s + VAAN_CP_s + (1|CollarID) + (1|Plot), data=full_stack)
plot(model1)
qqnorm(residuals(model1))
qqline(residuals(model1))
summary(model1)
# residuals are extremely non-normal