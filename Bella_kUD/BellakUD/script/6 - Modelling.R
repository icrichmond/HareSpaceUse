# Author: Isabella Richmond
# Last Edited: July 29, 2020

# This script is for the analysis of the effects of habitat complexity and food 
# quality on space use by snowshoe hare 


# load required packages 
easypackages::packages("tidyverse", "lme4","glmmTMB", "data.table", "AICcmodavg", "ggeffects",
                       "broom.mixed")

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
global_log <- lmerTest::lmer(logKUD ~ overPCA_s + underPCA_s + VAAN_CN_s + VAAN_CP_s + 
                     overPCA_s*underPCA_s + VAAN_CN_s*VAAN_CP_s + 
                     overPCA_s*VAAN_CN_s + overPCA_s*VAAN_CP_s + 
                     underPCA_s*VAAN_CN_s + underPCA_s*VAAN_CP_s + 
                     (1|CollarID) + (1|Plot), data=full_stack_s)
pred_log <- lmerTest::lmer(logKUD ~ overPCA_s + underPCA_s + overPCA_s*underPCA_s + (1|CollarID) + (1|Plot), data=full_stack_s)
stoich_log <- lmerTest::lmer(logKUD ~ VAAN_CN_s + VAAN_CP_s + VAAN_CN_s*VAAN_CP_s + (1|CollarID) + (1|Plot), data=full_stack_s)
null_log <- lmerTest::lmer(logKUD ~ 1 + (1|CollarID) + (1|Plot), data=full_stack_s)

# Use AICc to evaluate the competing hypotheses 
# create list of models 
models <- list(global_log, pred_log, stoich_log, null_log)
# use imap to loop through list of models using function at start of script and 
# create diagnostic figures 
source("script/function-residPlots.R")
models.residplots <- imap(models, resid_plots) 
# save all diagnostic plots to a pdf 
pdf("graphics/lmem_log_diagnostics.pdf")
models.residplots
dev.off()
# models are not perfect but they are good enough (or as good as they will be)
# create an AICc table to show the "best model"
modelsnames <- list("Mod 1 = Global" = global_log, "Mod 2 = Habitat Complexity" = pred_log, "Mod 3 = Food Quality" = stoich_log, "Mod 4 = Null" = null_log)
models.aic <- aictab(cand.set = modelsnames)
print(models.aic)
write.csv(models.aic, "output/lmem_log_aic.csv")
# Mod3 - stoich model, is top ranking model 
# save the summary tables of the models 
summary(stoich_log)
summary.models <-map_df(models, broom.mixed::tidy, .id="model")
write_csv(summary.models, path = "output/lmem_log_summary.csv")
# calculate pseudo R^2 of top model- just another check of significance determination
performance::r2_nakagawa(stoich_log)

# --------------------------------------- #
#             Random Slopes               #
# --------------------------------------- #
# investigating relationships between explanatory variables and logKUD 
# to see if random variables should be allowed to have random slopes as
# well as intercepts. Only doing this for CollarID because individuals 
# should have different reactions due to behavioural differences 
# but Plot would be expected to have similar slopes/no individuality. 
# Just want to control for random variation across plots, so stick to 
# random intercepts.

# looking at individual collar relationships with C:N and logKUD
ggplot(data = full_stack_s, aes(x = VAAN_CN_s, y = logKUD, col = CollarID, group = CollarID))+ 
  geom_point(size     = 1.2, alpha = .8, position = "jitter")+ 
  theme_minimal()+
  scale_color_manual(values = cols)+
  geom_smooth(method = lm, se = FALSE, size = .5, alpha = .8)
ggsave("graphics/VAANCN_kUD_collar_vis.png")
# pretty different across individuals - will change model so that there are random slopes 
# for CollarID with VAAN C:N

# now look at individual collar relationships with C:P and logKUD
ggplot(data = full_stack_s, aes(x = VAAN_CP_s, y = logKUD, col = CollarID, group = CollarID))+ 
  geom_point(size     = 1.2, alpha = .8, position = "jitter")+ 
  theme_minimal()+
  scale_color_manual(values = cols)+
  geom_smooth(method = lm, se = FALSE, size = .5, alpha = .8)
# the slopes of these lines look a lot more similar across individuals than with C:N
# only use random intercepts 

# now look at individual collar relationships with overstory complexity and logKUD
ggplot(data = full_stack_s, aes(x = overPCA_s, y = logKUD, col = CollarID, group = CollarID))+ 
  geom_point(size     = 1.2, alpha = .8, position = "jitter")+ 
  theme_minimal()+
  scale_color_manual(values = cols)+
  geom_smooth(method = lm, se = FALSE, size = .5, alpha = .8)
# the slopes of these lines look a lot more similar across individuals than with C:N
# only use random intercepts

# now look at individual collar relationships with overstory complexity and logKUD
ggplot(data = full_stack_s, aes(x = underPCA_s, y = logKUD, col = CollarID, group = CollarID))+ 
  geom_point(size     = 1.2, alpha = .8, position = "jitter")+ 
  theme_minimal()+
  scale_color_manual(values = cols)+
  geom_smooth(method = lm, se = FALSE, size = .5, alpha = .8)
# the slopes of these lines look very different. Include random slopes with CollarID

# --------------------------------------- #
#        LMEM with Random Slopes          #
# --------------------------------------- #
# linear mixed effect model with individual (CollarID) and plot as random effects
# assigning CollarID random intercepts AND slopes based on plots above 

# global model
global_log_slope <- lmerTest::lmer(logKUD ~ overPCA_s + underPCA_s + VAAN_CN_s + VAAN_CP_s + 
                               overPCA_s*underPCA_s + VAAN_CN_s*VAAN_CP_s + 
                               overPCA_s*VAAN_CN_s + overPCA_s*VAAN_CP_s + 
                               underPCA_s*VAAN_CN_s + underPCA_s*VAAN_CP_s + 
                               (1 + VAAN_CN_s | CollarID)  + (1 + underPCA_s|CollarID) + (1|Plot), data=full_stack_s)
# model is near singular - potentially overfitted 
pred_log_slope <- lmerTest::lmer(logKUD ~ overPCA_s + underPCA_s + overPCA_s*underPCA_s + (1 +  underPCA_s|CollarID) + (1|Plot), data=full_stack_s)
# model is near singular - potentially overfitted
# model does not converge if overPCA_s is included in random effect - less variable, leave out
stoich_log_slope <- lmerTest::lmer(logKUD ~ VAAN_CN_s + VAAN_CP_s + VAAN_CN_s*VAAN_CP_s + (1 + VAAN_CN_s|CollarID) + (1|Plot), data=full_stack_s)
# model does not converge if VAAN_CP_s is included in random effect - less variable, leave out
null_log_slope <- lmerTest::lmer(logKUD ~ 1 + (1|CollarID) + (1|Plot), data=full_stack_s)
# null model with random intercepts does not significantly differ from null model with 
# random slopes and intercepts and converges, move forward with this 

models_slope <- list(global_log_slope, pred_log_slope, stoich_log_slope, null_log_slope)
# create an AICc table to show the "best model"
modelsnames_slope <- list("Mod 1 = Global" = global_log_slope, "Mod 2 = Habitat Complexity" = pred_log_slope, "Mod 3 = Food Quality" = stoich_log_slope, "Mod 4 = Null" = null_log_slope)
models.aic_slope <- aictab(cand.set = modelsnames_slope)
print(models.aic_slope)
write.csv(models.aic_slope, "output/lmem_log_aic_slope.csv")
# Mod3 - stoich model, is top ranking model 
# save the summary tables of the models 
summary(stoich_log_slope_slope)
summary.models_slope <-map_df(models_slope, broom.mixed::tidy, .id="model")
write_csv(summary.models_slope, path = "output/lmem_log_slope_summary.csv")
# calculate pseudo R^2 of top model- just another check of significance determination
performance::r2_nakagawa(stoich_log_slope)

# compare models - one with random intercepts and one with random intercepts and slopes
anova(stoich_log, stoich_log_slope, refit=FALSE)
# stoich model with random slopes is a significantly better model fit than without
# present models with random slopes 

# --------------------------------------- #
#                Plotting                 #
# --------------------------------------- #
# plot the relationship between kUD and C:N for each individual 
# stoich was top model and C:N was the only variable where std. error wasn't larger than 
# the estimate - although it is still insignificant 
cols <- palette(rainbow(31))
# extract the prediction dataframe 
predicts <- ggpredict(stoich_log_slope, terms = c("VAAN_CN_s"))
ggplot(predicts)+
  geom_line(aes(x=x, y=predicted))+
  geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
              fill = "lightgrey", alpha = 0.5) +  # error band
  geom_point(data = full_stack_s,                      # adding the raw data (scaled values)
             aes(x = VAAN_CN_s, y = logKUD, colour = CollarID)) + 
  labs(x = "Lowbush Blueberry C:N (scaled)", y = "Kernel Utilization Distribution (log)", 
       title = NULL) + 
  theme_minimal()
ggsave("graphics/VAANCN_kUD_vis.png")

# look at the intercepts and slopes of the random effect of CollarID (individual)
ggpredict(stoich_log_slope, terms=c("VAAN_CP_s", "CollarID"), type = "re") %>%
  plot() + 
  labs(x = "Lowbush Blueberry C:N", y = "Kernel Utilization Distribution (log)", title=NULL)+
  scale_fill_manual(values = cols)
ggsave("graphics/CollarSlopes_CP_vis.png")