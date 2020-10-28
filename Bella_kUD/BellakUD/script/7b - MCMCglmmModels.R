# Author: Isabella Richmond
# Last Edited: October 28, 2020

# This script is for the analysis of the effects of habitat complexity and food 
# quality on space use by snowshoe hare using Bayesian statistics


# load required packages
easypackages::libraries("MCMCglmm", "data.table", "tidyverse", "coda", "parallel")

##################  DATA PREPARATION ##################

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

# omit the NAs so that MCMCglmm will work
full_stack_s2 <- na.omit(full_stack_s)

##################  MODELS ##################
# now we have plots with PCA values for complexity, 
# kernel utilization values for all 30 hares, 
# and stoich values for lowland blueberry C:N and C:P

### Global model
# uninformative prior. R = residual setup and G = random effect setup 
prior0 <- list(R = list(V = diag(1), nu = 6),
               G = list(G1 = list(V = diag(2), nu = 6),
                        G2 = list(V = diag(2), nu = 6),
                        G3 = list(V = diag(2), nu = 6),
                        G4 = list(V = diag(2), nu = 6)
                       ))
# us in random effect represents "unstructured" covariance matrix for CollarID - meaning we want to calculate the variance due to differences among individuals
# and the covariance between the variances
# us(1+var):individual fits a random intercept-slope model with a covariance term
# nitt = total number of iterations
# burin = discarded iterations
# thin = number of iterations to discard in between successive stored samples (reduces autocorrelation)
global <- MCMCglmm(kUD ~ overPCA_s + underPCA_s + VAAN_CN_s + VAAN_CP_s + 
                     overPCA_s*VAAN_CN_s + overPCA_s*VAAN_CP_s + underPCA_s*VAAN_CN_s + underPCA_s*VAAN_CP_s, 
                   random =~ us(1 + underPCA_s):CollarID + us(1 + VAAN_CN_s):CollarID + us(1 + overPCA_s):CollarID + us(1 + VAAN_CP_s):CollarID, 
                   #family = "gaussian",
                   prior = prior0,
                   nitt=420000,
                   burnin=20000,
                   thin=100,
                   verbose = TRUE,
                   data = full_stack_s2,
                   pr=T, saveX = TRUE,saveZ = TRUE)
saveRDS(global, "large/globalMCMC.RDS")

# show diagnostic plots for random variables and fixed effects to check autocorrelation
plot(global$VCV)
plot(global$Sol)
# look at Geweke plots to check convergence 
# want Z score to be within the confidence intervals (dashed lines on plots)
# if most points are within dashed lines, no evidence against convergence
geweke.diag(global$Sol)
geweke.plot(global$Sol)
# look at Gelman plots
# compares two parallel chains to test convergence. If both quantiles are appox 1.0, 
# effective convergence can be diagnosed
# create parallel chains 
globalchains <- mclapply(1:4, function(i){
  MCMCglmm(kUD ~ overPCA_s + underPCA_s + VAAN_CN_s + VAAN_CP_s + 
             overPCA_s*VAAN_CN_s + overPCA_s*VAAN_CP_s + underPCA_s*VAAN_CN_s + underPCA_s*VAAN_CP_s, 
           random =~ us(1 + underPCA_s):CollarID + us(1 + VAAN_CN_s):CollarID + us(1 + overPCA_s):CollarID + us(1 + VAAN_CP_s):CollarID, 
           family = "gaussian",
           prior = prior0,
           nitt=420000,
           burnin=20000,
           thin=100,
           verbose = TRUE,
           data = full_stack_s2,
           pr=T, saveX = TRUE,saveZ = TRUE)
  
}, mc.cores=1)
globalchains <- lapply(globalchains, function(m) m$Sol)
globalchains <- do.call(mcmc.list, globalchains)
saveRDS(globalchains, "large/globalchainsMCMC.RDS")

gelman.diag(globalchains)
gelman.plot(globalchains)
# all diagnostics are good

# summary of model
summary(global)


### Stoich model
prior1 <- list(R = list(V = diag(1), nu = 6),
               G = list(G1 = list(V = diag(2), nu = 6),
                        G2 = list(V = diag(2), nu = 6)
               ))

stoich <- MCMCglmm(kUD ~ VAAN_CN_s + VAAN_CP_s + VAAN_CN_s*VAAN_CP_s,
                   random =~ us(1 + VAAN_CN_s):CollarID + us(1 + VAAN_CP_s):CollarID, 
                   #family = "gaussian",
                   prior = prior1,
                   nitt=420000,
                   burnin=20000,
                   thin=100,
                   verbose = TRUE,
                   data = full_stack_s2,
                   pr=T, saveX = TRUE,saveZ = TRUE)
saveRDS(stoich, "large/stoichMCMC.RDS")

### Risk model
prior2 <- list(R = list(V = diag(1), nu = 6),
               G = list(G1 = list(V = diag(2), nu = 6),
                        G2 = list(V = diag(2), nu = 6)
               ))

pred <- MCMCglmm(kUD ~ overPCA_s + underPCA_s + overPCA_s*underPCA_s, 
                 random =~ us(1 + underPCA_s):CollarID + us(1 + overPCA_s):CollarID, 
                   #family = "gaussian",
                   prior = prior2,
                   nitt=420000,
                   burnin=20000,
                   thin=100,
                   verbose = TRUE,
                   data = full_stack_s2,
                   pr=T, saveX = TRUE,saveZ = TRUE)
saveRDS(pred, "large/predMCMC.RDS")

### Intercept model
prior3  <- list(R = list(V = diag(1), nu = 6),
                         G = list(G1 = list(V = diag(1), nu = 6)
                         ))

intercept <- MCMCglmm(kUD ~ 1, 
                      random =~ CollarID, 
                      #family = "gaussian",
                      prior = prior3,
                      nitt=420000,
                      burnin=20000,
                      thin=100,
                      verbose = TRUE,
                      data = full_stack_s2,
                      pr=T, saveX = TRUE,saveZ = TRUE)
saveRDS(intercept, "large/interceptMCMC.RDS")

##################  DEVIANCE INFORMATION CRITERIA (DIC) ##################
# perform DIC and rank all four models
# because DIC values are negative, we want to subtract the max value to get 
# the smaller delta DIC (not the one with the smaller absolute value)
# smaller DIC is preferred to larger DIC 
DIC <- data.table(DIC = c(global$DIC, stoich$DIC, pred$DIC, intercept$DIC))
DIC$deltaDIC <- DIC$DIC - min(DIC$DIC)
DIC
# intercept model is top ranked model


