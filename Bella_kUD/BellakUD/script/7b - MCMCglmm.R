# Author: Isabella Richmond
# Last Edited: October 13, 2020

# This script is for the analysis of the effects of habitat complexity and food 
# quality on space use by snowshoe hare using a Bayesian approach instead of 
# Frequentist (see 7 - Modelling.R for Frequentist approach)


# load required packages 
easypackages::packages("tidyverse", "glmmTMB", "data.table", "AICcmodavg", "ggeffects",
                       "broom.mixed", "ggpubr", "patchwork", "MCMCglmm", "lme4", "nadiv")

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
#          Univariate MCMCglmm            #
# --------------------------------------- #
# we are using MCMCglmm to build our linear mixed effects models 
# MCMCglmm uses a Monte Carlo Markov Chain and is a Bayesian approach 
# see Houslay & Wilson, 2017 for why we switched from Frequentist to Bayesian

# need to set prior 
prior0 <- list(R=list(V=diag(12), nu=0.002),
               G=list(G1=list(V=diag(2), nu=2,
                              alpha.mu = rep(0,2),
                              alpha.V=diag(25^2,2,2))))



