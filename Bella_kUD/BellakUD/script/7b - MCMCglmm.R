# Author: Isabella Richmond
# Last Edited: October 27, 2020

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
                   family = "gaussian",
                   prior = prior0,
                   nitt=420000,
                   burnin=20000,
                   thin=100,
                   verbose = TRUE,
                   data = full_stack_s2,
                   pr=T, saveX = TRUE,saveZ = TRUE)

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
                   family = "gaussian",
                   prior = prior1,
                   nitt=420000,
                   burnin=20000,
                   thin=100,
                   verbose = TRUE,
                   data = full_stack_s2,
                   pr=T, saveX = TRUE,saveZ = TRUE)

### Risk model
prior2 <- list(R = list(V = diag(1), nu = 6),
               G = list(G1 = list(V = diag(2), nu = 6),
                        G2 = list(V = diag(2), nu = 6)
               ))

pred <- MCMCglmm(kUD ~ overPCA_s + underPCA_s + overPCA_s*underPCA_s, 
                 random =~ us(1 + underPCA_s):CollarID + us(1 + overPCA_s):CollarID, 
                   family = "gaussian",
                   prior = prior2,
                   nitt=420000,
                   burnin=20000,
                   thin=100,
                   verbose = TRUE,
                   data = full_stack_s2,
                   pr=T, saveX = TRUE,saveZ = TRUE)

### Intercept model
prior3 <- prior2 <- list(R = list(V = diag(1), nu = 6),
                         G = list(G1 = list(V = diag(1), nu = 6)
                         ))

intercept <- MCMCglmm(kUD ~ 1, 
                      random =~ CollarID, 
                      family = "gaussian",
                      prior = prior3,
                      nitt=420000,
                      burnin=20000,
                      thin=100,
                      verbose = TRUE,
                      data = full_stack_s2,
                      pr=T, saveX = TRUE,saveZ = TRUE)

##################  DEVIANCE INFORMATION CRITERIA (DIC) ##################
# perform DIC and rank all four models
# because DIC values are negative, we want to subtract the max value to get 
# the smaller delta DIC (not the one with the smaller absolute value)
# smaller DIC is preferred to larger DIC 
DIC <- data.table(DIC = c(global$DIC, stoich$DIC, pred$DIC, intercept$DIC))
DIC$deltaDIC <- DIC$DIC - min(DIC$DIC)
DIC
# intercept model is top ranked model

##################  CORRELATIONS ##################
# extract data from global model
dfglobal <- data.table(Trait = attr(colMeans(global$Sol), "names"),
                 Value = colMeans(global$Sol)) 
dfglobal$Trait <- gsub(pattern = "s.",replacement = "s-", as.factor(dfglobal$Trait)) 
dfglobal$Trait <- gsub(pattern = "D.",replacement = "D-", as.factor(dfglobal$Trait)) 
dfglobal$Trait <- gsub(pattern = ").",replacement = ")-", as.factor(dfglobal$Trait)) 
# extract rows of interest
print(dfglobal, max = 257)
dfglobalunder <-  dfglobal[41:71,]
dfglobalcn <- dfglobal[103:133,]
dfglobalover <- dfglobal[165:195,]
dfglobalcp <- dfglobal[227:257,]
dfglobalunder[, c('Trait', 'collar' ,'ID') := tstrsplit(Trait, '-', fixed = TRUE)][, c("collar") := NULL]
dfglobalcn[, c('Trait', 'collar' ,'ID') := tstrsplit(Trait, '-', fixed = TRUE)][, c("collar") := NULL]
dfglobalover[, c('Trait', 'collar' ,'ID') := tstrsplit(Trait, '-', fixed = TRUE)][, c("collar") := NULL]
dfglobalcp[, c('Trait', 'collar' ,'ID') := tstrsplit(Trait, '-', fixed = TRUE)][, c("collar") := NULL]

cnglobal <- dfglobalcn[Trait=="VAAN_CN_s"][,c("Trait") := NULL]
setnames(cnglobal, "Value", "VAAN_CN_s")
cpglobal <- dfglobalcp[Trait=="VAAN_CP_s"][,c("Trait") := NULL]
setnames(cpglobal, "Value", "VAAN_CP_s")
underglobal <- dfglobalunder[Trait == "underPCA_s"][,c("Trait") := NULL]
setnames(underglobal, "Value", "underPCA_s")
overglobal <- dfglobalover[Trait == "overPCA_s"][,c("Trait") := NULL]
setnames(overglobal, "Value", "overPCA_s")
g1 <- merge(cnglobal,cpglobal, by = "ID")
g2 <- merge(underglobal,overglobal, by="ID")
g3 <- merge(g1, g2, by="ID")


# extract data from stoich model
dfstoich <- data.table(Trait = attr(colMeans(stoich$Sol), "names"),
                       Value = colMeans(stoich$Sol)) 
dfstoich$Trait <- gsub(pattern = "s.",replacement = "s-", as.factor(dfstoich$Trait)) 
dfstoich$Trait <- gsub(pattern = "D.",replacement = "D-", as.factor(dfstoich$Trait)) 
dfstoich$Trait <- gsub(pattern = ").",replacement = ")-", as.factor(dfstoich$Trait)) 
# extract rows of interest - C:N slopes
dfstoich[1:120]
dfCN <-  dfstoich[36:66,]
dfCN[, c('Trait', 'collar' ,'ID') := tstrsplit(Trait, '-', fixed = TRUE)][, c("collar") := NULL]
# extract rows of interest - C:P slopes
dfCP <-  dfstoich[98:128,]
dfCP[, c('Trait', 'collar' ,'ID') := tstrsplit(Trait, '-', fixed = TRUE)][, c("collar") := NULL]


# extract data from pred model
dfpred <- data.table(Trait = attr(colMeans(pred$Sol), "names"),
                       Value = colMeans(pred$Sol)) 
dfpred$Trait <- gsub(pattern = "s.",replacement = "s-", as.factor(dfpred$Trait)) 
dfpred$Trait <- gsub(pattern = "D.",replacement = "D-", as.factor(dfpred$Trait)) 
dfpred$Trait <- gsub(pattern = ").",replacement = ")-", as.factor(dfpred$Trait)) 
# extract rows of interest - understory
dfpred[1:100]
dfunder <-  dfpred[36:66,]
dfunder[, c('Trait', 'collar' ,'ID') := tstrsplit(Trait, '-', fixed = TRUE)][, c("collar") := NULL]
# extract rows of interest - overstory
dfover <-  dfpred[98:128,]
dfover[, c('Trait', 'collar' ,'ID') := tstrsplit(Trait, '-', fixed = TRUE)][, c("collar") := NULL]


# make dataframes and merge for plotting
cn <- dfCN[Trait=="VAAN_CN_s"][,c("Trait") := NULL]
setnames(cn, "Value", "VAAN_CN_s")
cp <- dfCP[Trait=="VAAN_CP_s"][,c("Trait") := NULL]
setnames(cp, "Value", "VAAN_CP_s")
under <- dfunder[Trait == "underPCA_s"][,c("Trait") := NULL]
setnames(under, "Value", "underPCA_s")
over <- dfover[Trait == "overPCA_s"][,c("Trait") := NULL]
setnames(over, "Value", "overPCA_s")
aa <- merge(cn,cp, by = "ID")
ab <- merge(under,over, by="ID")
ac <- merge(aa, ab, by="ID")

# plot correlation
ggplot(g3, aes(VAAN_CN_s, underPCA_s)) +
  geom_point(size = 2, alpha = 0.65) +
  geom_smooth(method = "lm", se = F) +
  ylab("Understory Complexity") +
  xlab("Lowbush Blueberry C:N") +
  theme(legend.position = 'none',
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key = element_blank(),
        axis.text=element_text(size=12, color = "black"),
        plot.title = element_text(size = 20),
        axis.title=element_text(size=20),
        strip.text = element_text(size=12,face = "bold"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1))

cor.test(g3$VAAN_CN_s, g3$underPCA_s)
