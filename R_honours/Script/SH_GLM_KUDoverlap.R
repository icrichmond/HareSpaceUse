# import dataset

shad_KUD_Overlap <- read_csv("input/PointCollarJoin_shad.csv")
view(shad_KUD_Overlap)

#linear regression 
#scatterplot
scatter.smooth(x=shad_KUD_Overlap$SH_Alpha_D, y=shad_KUD_Overlap$Join_Count_SH) 

#create general linear model, parametric test
# paramteric test assumes normality, heterogeneity and independence of residuals
SH_linearMod <- lm(shad_KUD_Overlap$Join_Count_SH ~ shad_KUD_Overlap$SH_Alpha_D) 
# build linear regression model on full data
print(SH_linearMod)


#test residuals for assumptions
#using super cool package that bella is very excited about :)

install.packages("performance")
library("performance")
# need interaction terms
check_normality(SH_linearMod)
# non - normality detected
check_heteroscedasticity(SH_linearMod)
check_homogeneity(SH_linearMod)
## result data not normal as expected

#checking for normality, and it failed no normal data here
qqnorm(resid(SH_linearMod))
qqline(resid(SH_linearMod))
hist(resid(SH_linearMod))

# we have a Poisson distribution.. maybe lets check

SH_poissonMod <- glm(shad_KUD_Overlap$Join_Count_SH~shad_KUD_Overlap$SH_Alpha_D, family = poisson)
qqnorm(resid(SH_poissonMod))
qqline(resid(SH_poissonMod))

check_normality(SH_poissonMod)

hist(resid(SH_poissonMod))

# binomial?, but this is either 0 or 1 so its not ideal for this data where we are looking at a scale of overlap from 1>N

# try different links, each family has default link fnc,, ie poisson default link is log.... this also didnt work for identity, logit or inverse (family = (link name))


summary(SH_poissonMod)
#no significant relationship between overlap and alpha diversity

anova(SH_poissonMod)

# linear regression scat plot
# GLM, poisson family, (y= mx+b) overlap dependent on alpha diversity plus error
# diagnostic plots for ^ look into dp more, we can have this for glm and show how poisson improved assumptions
#summary table for models
# intercept null model is trend that would occur randomly if no effect/relationship between x and y
# model is compared to null model and it is showing that there is no difference between null and my predicted model, the Pr(>z) value is not less than 0.05 so it is not significant " its kind of like the p model "
# make results and discussion slide from this, why no relationship?
