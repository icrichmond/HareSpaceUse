# Author: Isabella Richmond
# Last Edited: July 24, 2020

# This script is for the analysis of the effects of habitat complexity and food 
# quality on space use by snowshoe hare 


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

# stack the data so that it is ready for analysis - will include individual as 
# a random effect to control for pseudoreplication
full_stack <- pivot_longer(full, cols = starts_with("X"), names_to = "CollarID", names_prefix = "X", values_to = "kUD")
# fix dataset so appropriate columns are factors
full_stack <- full_stack %>% mutate(Plot = as.factor(Plot)) %>%
  mutate(CollarID = as.factor(CollarID))
head(full_stack)

# set zeroes in the dataset to NA 
full_stack[full_stack == 0] <- NA

# standardize the explanatory variables
full_stack_s <- full_stack %>%
  add_column(VAAN_CN_s = scale(full_stack$VAAN_CN, center = TRUE, scale = TRUE)) %>%
  add_column(VAAN_CP_s = scale(full_stack$VAAN_CP, center = TRUE, scale = TRUE)) %>%
  add_column(overPCA_s = scale(full_stack$overPCA, center = TRUE, scale = TRUE)) %>%
  add_column(underPCA_s = scale(full_stack$underPCA, center = TRUE, scale = TRUE))


# --------------------------------------- #
#   Linear Mixed Effect Models with kUD   #
# --------------------------------------- #
# now we have plots with PCA values for complexity, 
# kernel utilization values for all 30 hares, 
# and stoich values for lowland blueberry C:N and C:P

# linear mixed effect model with individual (CollarID) and plot as random effects
# global model
model1 <- lmer(kUD ~ overPCA_s + underPCA_s + VAAN_CN_s + VAAN_CP_s + (1|CollarID) + (1|Plot), data=full_stack_s)
plot(model1)
qqnorm(residuals(model1))
qqline(residuals(model1))
summary(model1)
# residuals are extremely non-normal

# habitat complexity model
model2 <- lmer(kUD ~ overPCA_s + underPCA_s + (1|CollarID) + (1|Plot), data=full_stack_s)
plot(model1)
qqnorm(residuals(model2))
qqline(residuals(model2))
summary(model2)
# residuals are extremely non-normal

# stoichiometry model 
model3 <- lmer(kUD ~ VAAN_CN_s + VAAN_CP_s + (1|CollarID) + (1|Plot), data=full_stack_s)
plot(model3)
qqnorm(residuals(model3))
qqline(residuals(model3))
summary(model3)
# residuals are extremely non-normal

# null model
