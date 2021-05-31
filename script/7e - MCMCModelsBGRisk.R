# Author: Isabella Richmond
# Last Edited: May 26, 2021

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
h2017 <- add_column(h2017, Year = 2017)
h2018 <- h2018 %>% 
  drop_na() %>%
  dplyr::select(Frequency, EarTag) %>%
  dplyr::group_by(Frequency, EarTag) %>%
  dplyr::summarise()
h2018 <- add_column(h2018, Year = 2018)
h2019 <- h2019 %>% 
  drop_na() %>%
  dplyr::select(Frequency, Eartag) %>%
  dplyr::group_by(Frequency, Eartag) %>%
  dplyr::summarise() %>%
  dplyr::rename(EarTag = Eartag)
h2019 <- add_column(h2019, Year = 2019)
h <- rbind(h2017, h2018, h2019)
hsex <- inner_join(h, sex, by = "EarTag")
hsex <- hsex %>%
  dplyr::group_by(EarTag, Frequency, Sex, Year) %>%
  dplyr::summarise() %>%
  drop_na() %>%
  dplyr::rename(CollarID = Frequency)
hsex$CollarID <- as.factor(hsex$CollarID)
full_stack_s <- inner_join(full_stack, hsex, by = "CollarID")
full_stack_s$Sex <- as.factor(full_stack_s$Sex)

# standardize the explanatory variables
full_stack_s <- full_stack_s %>%
  add_column(VAAN_CN_s = scale(full_stack_s$VAAN_CN, center = TRUE, scale = TRUE)) %>%
  add_column(VAAN_CP_s = scale(full_stack_s$VAAN_CP, center = TRUE, scale = TRUE)) %>%
  add_column(overPCA_s = scale(full_stack_s$overPCA, center = TRUE, scale = TRUE)) %>%
  add_column(underPCA_s = scale(full_stack_s$underPCA, center = TRUE, scale = TRUE))

# omit the NAs so that MCMCglmm will work
full_stack_s2 <- na.omit(full_stack_s)

##################  MODEL ##################
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
                     overPCA_s*VAAN_CN_s + overPCA_s*VAAN_CP_s + underPCA_s*VAAN_CN_s + 
                     underPCA_s*VAAN_CP_s + Sex + overPCA_s*Year + underPCA_s*Year, 
                   random =~ us(1 + underPCA_s):CollarID + us(1 + VAAN_CN_s):CollarID + 
                     us(1 + overPCA_s):CollarID + us(1 + VAAN_CP_s):CollarID, 
                   #family = "gaussian",
                   prior = prior0,
                   nitt=420000,
                   burnin=20000,
                   thin=100,
                   verbose = TRUE,
                   data = full_stack_s2,
                   pr=T, saveX = TRUE,saveZ = TRUE)
saveRDS(global, "large/globalMCMC_BGrisk.rds")

# summary of model
summary(global)

# get posterior means, 95% CI, and quantiles
m <- global$Sol
# select the variables we care about (risk * year)
m <- m[,13:14]
summary(m)

plot.estimates <- function(x) {
  if (class(x) != "summary.mcmc")
    x <- summary(x)
  n <- dim(x$statistics)[1]
  par(mar=c(2, 7, 4, 1))
  plot(x$statistics[,1], n:1,
       yaxt="n", ylab="",
       xlim=range(x$quantiles)*1.2,
       pch=19,
       main="Posterior means and 95% credible intervals")
  grid()
  axis(2, at=n:1, rownames(x$statistics), las=2)
  arrows(x$quantiles[,1], n:1, x$quantiles[,5], n:1, code=0)
  abline(v=0, lty=2)
}

plot.estimates(m)
