# Author: Isabella Richmond
# Last Edited: July 24, 2020

# This script is for the analysis of the effects of habitat complexity and food 
# quality on space use by snowshoe hare 


# load required packages 
easypackages::packages("tidyverse", "lme4","glmmTMB", "data.table", "AICcmodavg")

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
global <- glmmTMB(kUD ~ overPCA_s + underPCA_s + VAAN_CN_s + VAAN_CP_s + (1|CollarID) + (1|Plot), family = beta_family, data=full_stack_s)
qqnorm(residuals(global))
qqline(residuals(global))
hist(residuals(global))
summary(global)
# residuals are extremely non-normal

# habitat complexity model
pred <- glmmTMB(kUD ~ overPCA_s + underPCA_s + (1|CollarID) + (1|Plot), family = beta_family, data=full_stack_s)
plot(pred)
qqnorm(residuals(pred))
qqline(residuals(pred))
summary(pred)
# residuals are extremely non-normal

# stoichiometry model 
stoich <- glmmTMB(kUD ~ VAAN_CN_s + VAAN_CP_s + (1|CollarID) + (1|Plot), family = beta_family, data=full_stack_s)
plot(stoich)
qqnorm(residuals(stoich))
qqline(residuals(stoich))
summary(stoich)
# residuals are extremely non-normal

# null model
null <- glmmTMB(kUD ~ 1, data=full_stack_s)
plot(null)

# --------------------------------------- #
#           Data Transformation           #
# --------------------------------------- #
# ln and log10 transform the data to fix the wacky residuals 
# not ideal but no error distribution fits the data 
# log function in r is by default natural logarithm (ln)
full_stack_s <- full_stack_s %>%
  mutate(lnKUD = log(kUD)) %>%
  mutate(logKUD = log10(kUD))

globalt <- glmmTMB(logKUD ~ overPCA_s + underPCA_s + VAAN_CN_s + VAAN_CP_s + (1|CollarID) + (1|Plot), data=full_stack_s)
qqnorm(residuals(globalt))
qqline(residuals(globalt))
hist(residuals(globalt))

globalt <- glmmTMB(lnKUD ~ overPCA_s + underPCA_s + VAAN_CN_s + VAAN_CP_s + (1|CollarID) + (1|Plot), data=full_stack_s)
qqnorm(residuals(globalt))
qqline(residuals(globalt))
hist(residuals(globalt))

# proceeding with log10 transformation - tail is smaller and residuals
# are not perfect but have improved drastically. SJL said GLM is very
# robust to imperfect residuals. Could do hurdle model but a threshold 
# based hurdle model is more controversial.

# --------------------------------------- #
# Transformed Linear Mixed Effects Models #
# --------------------------------------- #
# linear mixed effect model with individual (CollarID) and plot as random effects
# global model
global_log <- lmer(logKUD ~ overPCA_s + underPCA_s + VAAN_CN_s + VAAN_CP_s + (1|CollarID) + (1|Plot), data=full_stack_s)
pred_log <- lmer(logKUD ~ overPCA_s + underPCA_s + (1|CollarID) + (1|Plot), data=full_stack_s)
stoich_log <- lmer(logKUD ~ VAAN_CN_s + VAAN_CP_s + (1|CollarID) + (1|Plot), data=full_stack_s)
null_log <- lmer(logKUD ~ (1|CollarID) + (1|Plot), data=full_stack_s)

# Use AICc to evaluate the competing hypotheses 
# create list of models 
models <- list(global_log, pred_log, stoich_log, null_log)
# use imap to loop through list of models using function at start of script and 
# create diagnostic figures 
#source("script/function-residPlots.R")
#models.residplots <- imap(models, resid_plots) 
# save all diagnostic plots to a pdf 
#pdf("graphics/StoichModels_2Step/ModelDiagnostics_GzLM/ABBA_Qty_P_gamma.pdf")
#ABBA.Qty_P.residplots
#dev.off()

# models do not meet assumptions - really bad fit. Will have to try a different error structure.
# create an AICc table to show the "best model" to use as a prediction of spatial stoichiometry
modelsnames <- list("Mod 1 = Global" = global_log, "Mod 2 = Habitat Complexity" = pred_log, "Mod 3 = Food Quality" = stoich_log, "Mod 4 = Null" = null_log)
models.aic <- aictab(cand.set = models)
print(models.aic)
write.csv(models.aic, "output/AIC_2Step/ABBA_Qty_P_gamma.csv")
# save the summary tables of the models 
summary.models <-map_df(models, broom::tidy, .id="model")
write_csv(summary.models, path = "output/Summary_2Step/summary.ABBA.Qty_P_gamma.csv")
# calculate pseudo R^2 - just another check of significance determination
performance::r2_nakagawa(pred)