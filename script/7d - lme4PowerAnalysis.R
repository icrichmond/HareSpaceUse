# Author: Isabella Richmond 
# This script is for testing the power of our dataset to detect individual-level
# reaction norms. The power analysis is based off Martin et al's (2011) paper. 
# Although it is a frequentist test and we use Bayesian modelling, we believe it 
# is useful to provide confidence in our sampling regime and analyses

#### Load Packages ####
easypackages::packages("lmerTest", "pamm")

#### Load Data ####
full_stack_s2 <- readRDS("large/full_stack_s2.rds")

#### lmer Models ####
global_pa <- lmerTest::lmer(kUD ~ overPCA_s + underPCA_s + VAAN_CN_s + VAAN_CP_s + 
                              overPCA_s*VAAN_CN_s + overPCA_s*VAAN_CP_s + 
                              underPCA_s*VAAN_CN_s + underPCA_s*VAAN_CP_s + 
                              (1 + VAAN_CN_s | CollarID)  + (1 + VAAN_CP_s | CollarID) + 
                              (1 + overPCA_s | CollarID) + (1 + underPCA_s | CollarID), 
                              data=full_stack_s2)
# model doesn't converge but that is okay because we are just using it for power analysis

#### Power Analysis ####
# use EAMM function from pamm package, from Martin et al (2011)
# test each environmental characteristic and individual collars 
pa_g_over <- EAMM(numsim = 10, mer.model = list(global_pa, env = "overPCA_s", random = "CollarID"))
pa_g_under <- EAMM(numsim = 10, mer.model = list(global_pa, env = "underPCA_s", random = "CollarID"))
pa_g_cn <- EAMM(numsim = 10, mer.model = list(global_pa, env = "VAAN_CN_s", random = "CollarID"))
pa_g_cp <- EAMM(numsim = 10, mer.model = list(global_pa, env = "VAAN_CP_s", random = "CollarID"))
# power is good - we can assess individual level reaction norms 
saveRDS(pa_g_over, "large/powerOverstory.rds")
saveRDS(pa_g_under, "large/powerUnderstory.rds")
saveRDS(pa_g_cn, "large/powerCN.rds")
saveRDS(pa_g_cp, "large/powerCP.rds")
