# Author: Isabella Richmond
# Last Edited: October 28, 2020

# this script is to plot the varying slopes and correlations between slopes from the 
# MCMCglmm models 
# a lot of the code was adapted from Quinn Webber (https://qwebber.weebly.com/)


# load required packages
easypackages::libraries("data.table", "tidyverse", "patchwork", "ggpubr", "MCMCglmm")

##################  DATA PREPARATION ##################
# dataset
full_stack_s2 <- readRDS("large/full_stack_s2.rds")
# models 
global <- readRDS("large/globalMCMC.rds")
stoich <- readRDS("large/stoichMCMC.RDS")
pred <- readRDS("large/predMCMC.RDS")
intercept <- readRDS("large/interceptMCMC.RDS")

##################  SLOPES ##################
# want to plot the intercepts and slopes of each random effect 

## CN:CollarID
df_cn <- cbind(full_stack_s2,fit = predict(stoich, marginal = NULL)) %>%
  dplyr::group_by(CollarID, VAAN_CN_s) %>%
  dplyr::summarise(fit = mean(fit),
            kUD = mean(kUD)) %>%
  tidyr::gather(Type, Value,
         fit:kUD)

df_fit_cn = setDT(df_cn)[Type == "fit"]
# extract slopes and intercepts for individuals
intCN = lm(Value ~ VAAN_CN_s, data = df_fit_cn)$coefficients[1]
slopeCN = lm(Value ~ VAAN_CN_s, data = df_fit_cn)$coefficients[2]
# plot
cols <- palette(rainbow(31))
n <- ggplot(df_fit_cn, aes(x = VAAN_CN_s, y = Value, group = factor(CollarID))) +
  geom_smooth(
    aes(VAAN_CN_s, Value, colour=CollarID),
    size = 0.25,
    method = lm,
    se = FALSE) +
  labs(x = "Lowbush Blueberry C:N", y = "Kernel Utilization Distribution", color="Individual") +
  scale_colour_grey(name="CollarID")+
  ggtitle('A')+
  theme(
    plot.title = element_text(size = 9),
    axis.text = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 9),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour="black"),
    axis.ticks = element_line(colour="black")
  )

#CP:CollarID
df_cp <- cbind(full_stack_s2,fit = predict(stoich, marginal = NULL)) %>%
  dplyr::group_by(CollarID, VAAN_CP_s) %>%
  dplyr::summarise(fit = mean(fit),
            kUD = mean(kUD)) %>%
  tidyr::gather(Type, Value,
         fit:kUD)

df_fit_cp = setDT(df_cp)[Type == "fit"]
# extract slopes and intercepts for individuals
intCP = lm(Value ~ VAAN_CP_s, data = df_fit_cp)$coefficients[1]
slopeCP = lm(Value ~ VAAN_CP_s, data = df_fit_cp)$coefficients[2]
# plot
p <- ggplot(df_fit_cp, aes(x = VAAN_CP_s, y = Value, group = factor(CollarID))) +
  geom_smooth(
    aes(VAAN_CP_s, Value, colour=CollarID),
    size = 0.25,
    method = lm,
    se = FALSE) +
  labs(x = "Lowbush Blueberry C:P", y = "Kernel Utilization Distribution", color="Individual") +
  scale_colour_grey(name="CollarID")+
  ggtitle('B')+
  theme(
    plot.title = element_text(size = 9),
    axis.text = element_text(size = 8, color = "black"),
    axis.title.y = element_blank(),
    axis.title = element_text(size = 9),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour="black"),
    axis.ticks = element_line(colour="black")
  )

# understory:CollarID
df_under <- cbind(full_stack_s2,fit = predict(pred, marginal = NULL)) %>%
  dplyr::group_by(CollarID, underPCA_s) %>%
  dplyr::summarise(fit = mean(fit),
            kUD = mean(kUD)) %>%
  tidyr::gather(Type, Value,
         fit:kUD)

df_fit_under = setDT(df_under)[Type == "fit"]
# extract slopes and intercepts for individuals
intunder = lm(Value ~ underPCA_s, data = df_fit_under)$coefficients[1]
slopeunder = lm(Value ~ underPCA_s, data = df_fit_under)$coefficients[2]
# plot
u <- ggplot(df_fit_under, aes(x = underPCA_s, y = Value, group = factor(CollarID))) +
  geom_smooth(
    aes(underPCA_s, Value, colour=CollarID),
    size = 0.25,
    method = lm,
    se = FALSE) +
  labs(x = "Understory Complexity", y = "Kernel Utilization Distribution", colour="Individual") +
  scale_colour_grey(name="CollarID")+
  ggtitle('C')+
  theme(
    plot.title = element_text(size = 9),
    axis.text = element_text(size = 8, color = "black"),
    axis.title = element_text(size = 9),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour="black"),
    axis.ticks = element_line(colour="black")
  )

#Overstory:CollarID
df_over <- cbind(full_stack_s2,fit = predict(pred, marginal = NULL)) %>%
  dplyr::group_by(CollarID, overPCA_s) %>%
  dplyr::summarise(fit = mean(fit),
            kUD = mean(kUD)) %>%
  tidyr::gather(Type, Value,
         fit:kUD)

df_fit_over = setDT(df_over)[Type == "fit"]
# extract slopes and intercepts for individuals
intover = lm(Value ~ overPCA_s, data = df_fit_over)$coefficients[1]
slopeover = lm(Value ~ overPCA_s, data = df_fit_over)$coefficients[2]
# plot
o <- ggplot(df_fit_over, aes(x = overPCA_s, y = Value, group = factor(CollarID))) +
  geom_smooth(
    aes(overPCA_s, Value, colour=CollarID),
    size = 0.25,
    method = lm,
    se = FALSE) +
  labs(x = "Overstory Complexity", y = "Kernel Utilization Distribution", color="Individual") +
  scale_colour_grey(name="CollarID")+
  ggtitle('D')+
  theme(
    plot.title = element_text(size = 9),
    axis.text = element_text(size = 8, color = "black"),
    axis.title.y = element_blank(),
    axis.title = element_text(size = 9),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour="black"),
    axis.ticks = element_line(colour="black")
  )

# plot all random slopes together using patchwork 
(n|p)/(u|o) + plot_layout(guides = 'collect')
ggsave("graphics/varyingslopesMCMC_grayscale.tiff", dpi = 600)

##################  CORRELATIONS ##################
# extract data from global model
dfglobal <- data.table(Trait = attr(colMeans(global$Sol), "names"),
                       Value = colMeans(global$Sol)) 
dfglobal$Trait <- gsub(pattern = "s.",replacement = "s-", as.factor(dfglobal$Trait)) 
dfglobal$Trait <- gsub(pattern = "D.",replacement = "D-", as.factor(dfglobal$Trait)) 
dfglobal$Trait <- gsub(pattern = ").",replacement = ")-", as.factor(dfglobal$Trait)) 
# extract rows of interest
print(dfglobal, max = 243)
dfglobalunder <-  dfglobal[41:69,]
dfglobalcn <- dfglobal[99:127,]
dfglobalover <- dfglobal[157:185,]
dfglobalcp <- dfglobal[215:243,]
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
dfCN <-  dfstoich[36:64,]
dfCN[, c('Trait', 'collar' ,'ID') := tstrsplit(Trait, '-', fixed = TRUE)][, c("collar") := NULL]
# extract rows of interest - C:P slopes
dfCP <-  dfstoich[94:122,]
dfCP[, c('Trait', 'collar' ,'ID') := tstrsplit(Trait, '-', fixed = TRUE)][, c("collar") := NULL]


# extract data from pred model
dfpred <- data.table(Trait = attr(colMeans(pred$Sol), "names"),
                     Value = colMeans(pred$Sol)) 
dfpred$Trait <- gsub(pattern = "s.",replacement = "s-", as.factor(dfpred$Trait)) 
dfpred$Trait <- gsub(pattern = "D.",replacement = "D-", as.factor(dfpred$Trait)) 
dfpred$Trait <- gsub(pattern = ").",replacement = ")-", as.factor(dfpred$Trait)) 
# extract rows of interest - understory
dfpred[1:100]
dfunder <-  dfpred[36:64,]
dfunder[, c('Trait', 'collar' ,'ID') := tstrsplit(Trait, '-', fixed = TRUE)][, c("collar") := NULL]
# extract rows of interest - overstory
dfover <-  dfpred[94:122,]
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

# plot correlations from stoich and predation models 
# choosing to correlate random effects from stoich and predation models 
# instead of global model because global model ranked last in DIC - worst model

# CN understory correlation 
cnunder <- ggplot(ac, aes(VAAN_CN_s, underPCA_s)) +
  geom_point(size = 2, alpha = 0.65) +
  geom_smooth(method = "lm", se = F, color="black") +
  stat_cor(aes(VAAN_CN_s,underPCA_s,label = paste(..r.label.., ..p.label.., sep = "~`, `~")),  label.y = 0.19)+
  stat_regline_equation(aes(VAAN_CN_s,underPCA_s), label.y = 0.22)+
  ylab("Understory Complexity") +
  xlab("Lowbush Blueberry C:N") +
  ggtitle("D")+
  theme(legend.position = 'none',
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.key = element_blank(),
        axis.text=element_text(size=12, color = "black"),
        plot.title = element_text(size = 12),
        axis.title=element_text(size=12),
        axis.title.y = element_blank(),
        strip.text = element_text(size=12,face = "bold"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour="black")
  )

cpunder <- ggplot(ac, aes(VAAN_CP_s, underPCA_s)) +
  geom_point(size = 2, alpha = 0.65) +
  geom_smooth(method = "lm", se = F, color="black") +
  stat_cor(aes(VAAN_CP_s,underPCA_s,label = paste(..r.label.., ..p.label.., sep = "~`, `~")), label.y = 0.19)+
  stat_regline_equation(aes(VAAN_CN_s,underPCA_s), label.y = 0.22)+
  ylab("Understory Complexity") +
  xlab("Lowbush Blueberry C:P") +
  ggtitle("C")+
  theme(legend.position = 'none',
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.key = element_blank(),
        axis.text=element_text(size=12, color = "black"),
        plot.title = element_text(size = 12),
        axis.title=element_text(size=12),
        strip.text = element_text(size=12,face = "bold"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour="black")
  )

cnover <- ggplot(ac, aes(VAAN_CN_s, overPCA_s)) +
  geom_point(size = 2, alpha = 0.65) +
  geom_smooth(method = "lm", se = F, color="black") +
  stat_cor(aes(VAAN_CN_s,overPCA_s,label = paste(..r.label.., ..p.label.., sep = "~`, `~")), label.y = 0.19)+
  stat_regline_equation(aes(VAAN_CN_s,overPCA_s), label.y = 0.22)+
  ylab("Overstory Complexity") +
  xlab("Lowbush Blueberry C:N") +
  ggtitle("B")+
  theme(legend.position = 'none',
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.key = element_blank(),
        axis.text=element_text(size=12, color = "black"),
        plot.title = element_text(size = 12),
        axis.title = element_blank(),
        #axis.title=element_text(size=12),
        strip.text = element_text(size=12,face = "bold"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour="black")
  )

cpover <- ggplot(ac, aes(VAAN_CP_s, overPCA_s)) +
  geom_point(size = 2, alpha = 0.65) +
  geom_smooth(method = "lm", se = F, color="black") +
  stat_cor(aes(VAAN_CP_s,overPCA_s,label = paste(..r.label.., ..p.label.., sep = "~`, `~")), label.y = 0.19)+
  stat_regline_equation(aes(VAAN_CP_s,overPCA_s), label.y = 0.22)+
  ylab("Overstory Complexity") +
  xlab("Lowbush Blueberry C:P") +
  ggtitle("A")+
  theme(legend.position = 'none',
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        legend.key = element_blank(),
        axis.text=element_text(size=12, color = "black"),
        plot.title = element_text(size = 12),
        axis.title=element_text(size=12),
        axis.title.x = element_blank(),
        strip.text = element_text(size=12,face = "bold"),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour="black")
  )

# combine plots 
(cpover|cnover)/(cpunder|cnunder)
ggsave("graphics/slopecorrelationsMCMC_greyscale.pdf", dpi=400)

##################  POPULATION RELATIONSHIPS ##################
# plot the relationship between kUD and explanatory variables with line of best fit 
# assessing the relationship at the population level
cn <- ggplot(full_stack_s2)+
  geom_point(aes(VAAN_CN_s, kUD),colour="grey46")+
  geom_smooth(aes(VAAN_CN_s, kUD), color = "grey3",method = lm)+
  stat_cor(aes(VAAN_CN_s,kUD,label = paste(..r.label..)), label.y = 0.93)+
  stat_regline_equation(aes(VAAN_CN_s,kUD), label.y = 0.99)+
  ggtitle("A")+
  theme(panel.background = element_blank(), axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"), axis.text = element_text(colour="black"))+
  labs(x = "Lowbush Blueberry C:N", y = "Kernel Utilization Distribution")

cp <- ggplot(full_stack_s2)+
  geom_point(aes(VAAN_CP_s, kUD), colour="grey46")+
  geom_smooth(aes(VAAN_CP_s, kUD), color = "grey3",method = lm)+
  stat_cor(aes(VAAN_CP_s,kUD,label = paste(..r.label..)), label.y = 0.93)+
  stat_regline_equation(aes(VAAN_CP_s,kUD), label.y = 0.99)+
  ggtitle("C")+
  theme(panel.background = element_blank(), axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"), axis.text = element_text(colour="black"))+
  labs(x = "Lowbush Blueberry C:P", y = "Kernel Utilization Distribution")

over <- ggplot(full_stack_s2)+
  geom_point(aes(overPCA_s, kUD), colour="grey46")+
  geom_smooth(aes(overPCA_s, kUD), color = "grey3",method = lm)+
  stat_cor(aes(overPCA_s,kUD,label = paste(..r.label..)), label.y = 0.93)+
  stat_regline_equation(aes(overPCA_s, kUD), label.y = 0.99)+
  ggtitle("B")+
  theme(panel.background = element_blank(), axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"), axis.text = element_text(colour="black"))+
  labs(x = "Overstory Habitat Complexity", y = " ")

under <- ggplot(full_stack_s2)+
  geom_point(aes(underPCA_s, kUD), colour="grey46")+
  geom_smooth(aes(underPCA_s, kUD), color = "grey3",method = lm)+
  stat_cor(aes(underPCA_s,kUD,label = paste(..r.label..)), label.y = 0.93)+
  stat_regline_equation(aes(underPCA_s,kUD), label.y = 0.99)+
  ggtitle("D")+
  theme(panel.background = element_blank(), axis.line = element_line(colour="black"),
        axis.ticks = element_line(colour="black"), axis.text = element_text(colour="black"))+
  labs(x = "Understory Habitat Complexity", y = " ")

(cn | over)/(cp | under)
ggsave("graphics/KUD_explanatory_MCMC.pdf", dpi = 400)
