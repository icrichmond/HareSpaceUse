# import dataset
KUD_Overlap_shad <- read_csv("input/shad.csv")
view(KUD_Overlap_shad)

#linear regression 
#scatterplot
scatter.smooth(x=KUD_Overlap_shad$sh_alpha_d
, y=KUD_Overlap_shad$Join_Count) 

#create general linear model, parametric test
# paramteric test assumes normality, heterogeneity and independence of residuals
linearMod <- lm(KUD_Overlap$Join_Count ~ KUD_Overlap$Alpha_D) 
# build linear regression model on full data
print(linearMod)


#test residuals for assumptions
#using super cool package that bella is very excited about :)

install.packages("performance")
library("performance")
# need interaction terms
check_normality(linearMod)
# non - normality detected
check_heteroscedasticity(linearMod)
check_homogeneity(linearMod)
## result data not normal as expected

#checking for normality, and it failed no normal data here
qqnorm(resid(linearMod))
qqline(resid(linearMod))
hist(resid(linearMod))

# we have a Poisson distribution.. maybe lets check

poissonMod <- glm(KUD_Overlap$Join_Count~KUD_Overlap$Alpha_D, family = poisson)
qqnorm(resid(poissonMod))
qqline(resid(poissonMod))
check_normality(poissonMod)

hist(resid(poissonMod))

# binomial?, but this is either 0 or 1 so its not ideal for this data where we are looking at a scale of overlap from 1>N

# try different links, each family has default link fnc,, ie poisson default link is log.... this also didnt work for identity, logit or inverse (family = (link name))


summary(poissonMod)
#no significant relationship between overlap and alpha diversity

anova(poissonMod)

# linear regression scat plot
# GLM, poisson family, (y= mx+b) overlap dependent on alpha diversity plus error
# diagnostic plots for ^ look into dp more, we can have this for glm and show how poisson improved assumptions
#summary table for models
  # intercept null model is trend that would occur randomly if no effect/relationship between x and y
  # model is compared to null model and it is showing that there is no difference between null and my predicted model, the Pr(>z) value is not less than 0.05 so it is not significant " its kind of like the p model "
  # make results and discussion slide from this, why no relationship?
