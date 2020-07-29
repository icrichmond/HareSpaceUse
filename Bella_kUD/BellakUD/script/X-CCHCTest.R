# this is to compare models with just hc and cc as opposed to the PCA axes
hc <- read_csv("input/CSRawData_HorizontalComplexity.csv")
hcsum <- hc %>%
  group_by(Distance, Score) %>%
  tally()
hcsum <- hcsum %>% 
  mutate(Score = as.factor(Score))
ggplot(hcsum, aes(x=Distance,y=n, fill=Score))+
  geom_col()+
  scale_fill_brewer(palette = 1)+
  theme(panel.background = element_rect(fill = "darkgrey"))+
  labs(fill="Score")
# inspection of the graph shows that a distance of 10 metres experiences the most
# variation in horizontal complexity scores
# therefore, use values from 10 m only going forward 
hc <- hc %>% filter(Distance == 10)
# terrestrial predators are all under 1 m tall (coyote & lynx), subset these heights
# from the dataset and average for each plot 
hc <- hc %>% filter(Height %in% c(0.5,1))
# calculate an average horizontal complexity score for each plot
# to do this we have to change the Score variable to an integer
hc <- hc %>% mutate(score = as.numeric(Score))
hcmean <- hc %>%
  group_by(Plot) %>%
  dplyr::summarise(meanhc = mean(score))

# canopy closure
cc <- read_csv('input/CSRawData_CanopyClosure.csv')
# join canopy closure and horiz complexity
cchc <- inner_join(cc, hcmean, by="Plot")
# read in stoich and KUD values 
stoichkud <- read_csv("output/cs_stoich_kud.csv")
# join all data 
full <- inner_join(cchc, stoichkud, by="Plot")
# remove X1 column 
full <- full %>% select(-X1)

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
  add_column(CoverValue_s = scale(full_stack$CoverValue, center = TRUE, scale = TRUE)) %>%
  add_column(meanhc_s = scale(full_stack$meanhc, center = TRUE, scale = TRUE))
# transform response variable so it meets assumptions
full_stack_s <- full_stack_s %>%
  mutate(lnKUD = log(kUD)) %>%
  mutate(logKUD = log10(kUD))

# --------------------------------------- #
#        LMEM with Random Slopes          #
# --------------------------------------- #
# linear mixed effect model with individual (CollarID) and plot as random effects
# assigning CollarID random intercepts AND slopes based on plots above 

# global model
global_log_slope <- lmerTest::lmer(logKUD ~ meanhc_s + CoverValue_s + VAAN_CN_s + VAAN_CP_s + 
                                     meanhc_s*CoverValue_s + VAAN_CN_s*VAAN_CP_s + 
                                     meanhc_s*VAAN_CN_s + meanhc_s*VAAN_CP_s + 
                                     CoverValue_s*VAAN_CN_s + CoverValue_s*VAAN_CP_s + 
                                     (1 + VAAN_CN_s | CollarID)  + (1 + CoverValue_s|CollarID) + (1|Plot), data=full_stack_s)
# model is near singular - potentially overfitted 
pred_log_slope <- lmerTest::lmer(logKUD ~ meanhc_s + CoverValue_s + meanhc_s*CoverValue_s + (1 +  CoverValue_s|CollarID) + (1|Plot), data=full_stack_s)
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
summary(stoich_log_slope)
summary.models_slope <-map_df(models_slope, broom.mixed::tidy, .id="model")
write_csv(summary.models_slope, path = "output/lmem_log_slope_summary.csv")
# calculate pseudo R^2 of top model- just another check of significance determination
performance::r2_nakagawa(stoich_log_slope)


