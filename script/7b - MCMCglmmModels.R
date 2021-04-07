# Author: Isabella Richmond
# Last Edited: October 29, 2020

# This script is for the analysis of the effects of habitat complexity and food 
# quality on space use by snowshoe hare using Bayesian statistics

# Quinn Webber helped with most of this code - https://qwebber.weebly.com/

# load required packages
easypackages::libraries("MCMCglmm", "data.table", "tidyverse", "coda", "parallel", "MuMIn")

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

# add sex of each individual to dataset 
sex <- read.csv("input/Trapping.csv")
h2017 <- read.csv("input/VHF_CleanData_2017.csv")
h2018 <- read.csv("input/VHF_CleanData_2018.csv")
h2019 <- read.csv("input/VHF_CleanData_2019.csv")
sex <- sex %>% 
  dplyr::rename(EarTag = Ear_tag) %>%
  dplyr::group_by(EarTag, Sex) %>%
  dplyr::summarise() %>%
  drop_na()
sex <- sex[-c(1:2,129:131),] # remove rows without identifiers 
sex1 <- sex[31:126,]
sex2 <- transform(sex[1:30,], EarTag = sprintf('A%s',EarTag)) # put A in front of EarTag # for rows 1:32
sex <- rbind(sex1,sex2)
sex$Sex[sex$Sex == ""] <- NA # set blank values as NA 
sex$Sex[sex$Sex == "Escaped"] <- NA # set values with escaped designation as NA 
sex <- drop_na(sex)
sex <- sex[!duplicated(sex), ] #remove duplicate rows
h2017 <- h2017 %>% 
  drop_na() %>%
  dplyr::select(Frequency, EarTag) %>%
  dplyr::group_by(Frequency, EarTag) %>%
  dplyr::summarise()
h2018 <- h2018 %>% 
  drop_na() %>%
  dplyr::select(Frequency, EarTag) %>%
  dplyr::group_by(Frequency, EarTag) %>%
  dplyr::summarise()
h2019 <- h2019 %>% 
  drop_na() %>%
  dplyr::select(Frequency, Eartag) %>%
  dplyr::group_by(Frequency, Eartag) %>%
  dplyr::summarise() %>%
  dplyr::rename(EarTag = Eartag)
h <- rbind(h2017, h2018, h2019)
hsex <- inner_join(h, sex, by = "EarTag")
hsex <- hsex %>%
  dplyr::group_by(EarTag, Frequency, Sex) %>%
  dplyr::summarise() %>%
  drop_na() %>%
  dplyr::rename(CollarID = Frequency)
hsex$CollarID <- as.factor(hsex$CollarID)
full_stack_s <- inner_join(full_stack, hsex, by = "CollarID")

# standardize the explanatory variables
full_stack_s <- full_stack_s %>%
  add_column(VAAN_CN_s = scale(full_stack_s$VAAN_CN, center = TRUE, scale = TRUE)) %>%
  add_column(VAAN_CP_s = scale(full_stack_s$VAAN_CP, center = TRUE, scale = TRUE)) %>%
  add_column(overPCA_s = scale(full_stack_s$overPCA, center = TRUE, scale = TRUE)) %>%
  add_column(underPCA_s = scale(full_stack_s$underPCA, center = TRUE, scale = TRUE))

# omit the NAs so that MCMCglmm will work
full_stack_s2 <- na.omit(full_stack_s)
# save final dataframe 
saveRDS(full_stack_s2, "large/full_stack_s2.rds")
rm(list=ls())
full_stack_s2 <- readRDS("large/full_stack_s2.rds")
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
saveRDS(global, "large/globalMCMC.rds")

# show diagnostic plots for random variables and fixed effects to check autocorrelation
pdf("graphics/traceplots_global.pdf")
plot(global$VCV)
plot(global$Sol)
dev.off()
# look at Geweke plots to check convergence 
# want Z score to be within the confidence intervals (dashed lines on plots)
# if most points are within dashed lines, no evidence against convergence
geweke.diag(global$Sol)
pdf("graphics/gewekeplots_global.pdf")
geweke.plot(global$Sol)
dev.off()
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
pdf("graphics/gelmanplots_global.pdf")
gelman.plot(globalchains)
dev.off()
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
# show diagnostic plots for random variables and fixed effects to check autocorrelation
pdf("graphics/traceplots_stoich.pdf")
plot(stoich$VCV)
plot(stoich$Sol)
dev.off()
# look at Geweke plots to check convergence 
# want Z score to be within the confidence intervals (dashed lines on plots)
# if most points are within dashed lines, no evidence against convergence
geweke.diag(stoich$Sol)
pdf("graphics/gewekeplots_stoich.pdf")
geweke.plot(stoich$Sol)
dev.off()

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
# show diagnostic plots for random variables and fixed effects to check autocorrelation
pdf("graphics/traceplots_pred.pdf")
plot(pred$VCV)
plot(pred$Sol)
dev.off()
# look at Geweke plots to check convergence 
# want Z score to be within the confidence intervals (dashed lines on plots)
# if most points are within dashed lines, no evidence against convergence
geweke.diag(pred$Sol)
pdf("graphics/gewekeplots_pred.pdf")
geweke.plot(pred$Sol)
dev.off()


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
# get full DIC table from MuMIn 
DICtable <- model.sel(global,stoich,pred,intercept,rank="DIC")
# intercept model is top ranked model