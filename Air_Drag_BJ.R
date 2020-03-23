require(rstan)
require(brms)

#####optimize setup for Bayesian Linear Models (rstan/brms)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7')

libraries <- c("cowplot","tidyverse","grid","gridExtra","lme4","lmerTest","ggpubr","apastats","rstatix")
# Where is this script?
invisible(lapply(libraries, function(x) {
  if(!require(x, character.only = T, quietly = T)) {
    install.packages(x)
    require(x, character.only = T)
  }
}
))
rm(libraries)

Where_Am_I <- function(path=T){
  if (path == T){
    dirname(rstudioapi::getSourceEditorContext()$path)
  }
  else {
    rstudioapi::getSourceEditorContext()$path
  }
}

theme_set(theme_cowplot())

setwd(Where_Am_I())

source("Utilities/Funs.R")

collapsed <- read.table(file = "Data/All_Data_Air_Drag.txt", header = T)
Continuous <- read.table(file = "Data/ContinuousData.txt", header = T)

air_drag <-  collapsed %>%
  select(trial,x_max,y_max,t_max,vx,
         x_max_model,t_max_model,y_max_model,
         G,ball,condsize,r,air_drag,label,
         random_x,rtime_timing,rtime_spatial,ball_x_spatial,ball_x_timing,
         id,TTC,condsize,visible) %>%
  # Raw differneces
  mutate(terror = rtime_timing-t_max-0.049,
         xerror = ball_x_spatial - x_max,
         OccludedDuration = t_max-visible,
         OccludedPercentage = visible/t_max,
         OccludedDistance = x_max-(x_max/2+vx*0.8*t_max*(OccludedPercentage-0.5)),
         # Differences with model modelling AD at oclusion
         xerror_model = ball_x_spatial - x_max_model,
         terror_model = rtime_timing-t_max_model-0.049,
         # Ratio Raw
         xerrorratio = (OccludedDistance+xerror) / OccludedDistance,
         terrorratio = (OccludedDuration+terror) / OccludedDuration,
         # Ratio Model -
         xerror_ratio_model = xerror_model / x_max_model,
         terror_ratio_model = terror_model / t_max_model,
         condsize = factor(condsize,levels = c("cong","incongr"),
                            labels = c("Congruent","Incongruent")),
         ball = factor(ball,levels = c("tennis","basket"),
                       labels = c("Tennis","Basket")),
         airdrag = case_when(air_drag == 1 ~ "Airdrag",
                             air_drag == 0 ~ "NoAirdrag"))

nAllTrials = length(air_drag$trial)

air_drag = air_drag %>% 
  group_by(id,label) %>%
  filter(terror > mean(terror)-2.5*sd(terror) & terror < mean(terror)+2.5*sd(terror),
         xerror > mean(xerror)-2.5*sd(xerror) & xerror < mean(xerror)+2.5*sd(xerror)) %>%
  filter(trim(terror, filter = T)) %>%
  filter(trim(xerror, filter = T)) %>%
  mutate(sd_time = sd(terror),
         sd_space = sd(xerror))

#number of excluded trials
nAllTrials - length(air_drag$trial)

#####Look at airdrag and stuff
ggplot(air_drag[air_drag$condsize == "Congruent",], aes(airdrag,terrorratio)) +
  geom_violin()

H1_Temporal_lmer <- lmer(terrorratio ~ airdrag + (1|id), 
                         data = air_drag)
H1_Temporal_lmer_NullModel <- lmer(terrorratio ~ (1|id), 
                                   data = air_drag)
anova(H1_Temporal_lmer,H1_Temporal_lmer_NullModel)
summary(H1_Temporal_lmer)

H1_Temporal_lmer <- lmer(terrorratio-1 ~ 0 + (1|id), 
                         data = air_drag[air_drag$airdrag == "Airdrag",])
H1_Temporal_lmer_NullModel <- lmer(terrorratio-1 ~ (1|id), 
                         data = air_drag[air_drag$airdrag == "Airdrag",])
anova(H1_Temporal_lmer,H1_Temporal_lmer_NullModel)
summary(H1_Temporal_lmer_NullModel)

H1_Temporal_lmer <- lmer(terrorratio-1 ~ 0 + (1|id), 
                         data = air_drag[air_drag$airdrag == "NoAirdrag",])
H1_Temporal_lmer_NullModel <- lmer(terrorratio-1 ~ (1|id), 
                                   data = air_drag[air_drag$airdrag == "NoAirdrag",])
anova(H1_Temporal_lmer,H1_Temporal_lmer_NullModel)
summary(H1_Temporal_lmer_NullModel)



get_prior(terrorratio ~ airdrag + (1|id),
          data = air_drag, family = gaussian())

fit1 <- brm(bf(terrorratio ~ airdrag + (1|id), 
               sigma ~ airdrag + (1|id)),
            data = air_drag, family = gaussian(),prior = set_prior("normal(1,10)", class = "b"))
fit1$prior
fit$prior
plot(Hypotheses)

Hypotheses = hypothesis(fit1,c("airdragNoAirdrag < 0",
                               "sigma_airdragNoAirdrag < 0",
                               "abs(Intercept-1) > abs(Intercept+airdragNoAirdrag-1)"))

fit2 <- brm(bf(xerrorratio ~ airdrag + (1|id), 
               sigma ~ airdrag + (1|id)),
            data = air_drag, family = gaussian())

Hypotheses = hypothesis(fit2,c("airdragNoAirdrag < 0",
                               "sigma_airdragNoAirdrag < 0",
                               "Intercept = 1",
                               "Intercept+airdragNoAirdrag = 1"))

ggplot(air_drag, aes(as.factor(condsize),t_max_model-t_max, color = as.factor(air_drag))) +
  geom_point() +
  facet_grid(.~ball)


H1_Spatial_TestModel <- lmer(xerrorratio ~ airdrag + (1|id), 
                         data = air_drag)
H1_Spatial_NullModel <- lmer(xerrorratio ~ (1|id), 
                                   data = air_drag)
anova(H1_Spatial_TestModel,H1_Spatial_NullModel)
summary(H1_Spatial_TestModel)

H1_Temporal_lmer <- lmer(xerrorratio-1 ~ 0 + (1|id), 
                         data = air_drag[air_drag$airdrag == "Airdrag",])
H1_Temporal_lmer_NullModel <- lmer(xerrorratio-1 ~ (1|id), 
                                   data = air_drag[air_drag$airdrag == "Airdrag",])
anova(H1_Temporal_lmer,H1_Temporal_lmer_NullModel)
summary(H1_Temporal_lmer_NullModel)

H1_Temporal_lmer <- lmer(xerrorratio-1 ~ 0 + (1|id), 
                         data = air_drag[air_drag$airdrag == "NoAirdrag",])
H1_Temporal_lmer_NullModel <- lmer(xerrorratio-1 ~ (1|id), 
                                   data = air_drag[air_drag$airdrag == "NoAirdrag",])
anova(H1_Temporal_lmer,H1_Temporal_lmer_NullModel)
summary(H1_Temporal_lmer_NullModel)


ggplot(air_drag,aes(terrorratio,xerrorratio, color = condsize)) +
  geom_point(alpha=0.2) +
  coord_fixed() +
  geom_smooth() +
  xlim(c(0,2))


######Look at ball size and congruency
##time
ggplot(air_drag, aes(condsize,terror_ratio, color = ball)) +
  geom_violin() +
  binomial_smooth()


SizeCongruency_TestModel_Time <- lmer(terrorratio ~ ball*condsize + (1|id), 
                         data = air_drag[air_drag$airdrag == "Airdrag",])
SizeCongruency_NullModel_Time <- lmer(terrorratio ~ ball + condsize + (1|id), 
                                      data = air_drag[air_drag$airdrag == "Airdrag",])
anova(SizeCongruency_TestModel_Time,SizeCongruency_NullModel_Time)
summary(SizeCongruency_TestModel_Time)
#(lol) (jaja)

##space
ggplot(air_drag, aes(ball,xerrorratio, color = condsize)) +
  geom_violin() +
  geom_boxplot()
SizeCongruency_TestModel_Space <- lmer(xerror_ratio ~ ball*condsize + (1|id), 
                                       data = air_drag[air_drag$airdrag == "Airdrag",])
SizeCongruency_NullModel_Space <- lmer(xerror_ratio ~ ball + condsize + (1|id), 
                                       data = air_drag[air_drag$airdrag == "Airdrag",])
anova(SizeCongruency_TestModel_Space,SizeCongruency_NullModel_Space)
summary(SizeCongruency_TestModel_Space)
#interesting




#####Hypothesis 1
#Hypothesis 1a
fit1 <- brm(bf(terrorratio ~ airdrag + (1|id),
               sigma ~ airdrag + (1|id)),
            data = air_drag, family = gaussian())

Hypotheses = hypothesis(fit1,c("airdragNoAirdrag < 0",
                               "abs(Intercept-1) < abs(Intercept+airdragNoAirdrag-1)"))

#Hypothesis 1a
fit2 <- brm(bf(xerrorratio ~ airdrag + (1|id),
               sigma ~ airdrag + (1|id)),
            data = air_drag, family = gaussian())

Hypotheses = hypothesis(fit2,c("airdragNoAirdrag < 0",
                               "abs(Intercept-1) < abs(Intercept+airdragNoAirdrag-1)"))



#####Hypothesis 2
#Hypothesis 2a
#
fit3 <- brm(bf(terrorratio ~ ball + condsize + (1|id),
               sigma ~ ball + condsize + (1|id)),
            data = air_drag, family = gaussian())

Hypotheses = hypothesis(fit3,c("sigma_ballBasket > 0",
                               "sigma_condsizeIncongruent < 0",
                               "sigma_ballBasket:condsizeIncongruent < 0"))

fit4 <- brm(bf(xerrorratio ~ ball*condsize + (1|id),
               sigma ~ ball*condsize + (1|id)),
            data = air_drag, family = gaussian())
Hypotheses = hypothesis(fit4,c("sigma_ballBasket > 0",
                               "sigma_condsizeIncongruent < 0",
                               "sigma_ballBasket:condsizeIncongruent < 0"))

fit5 <- brm(bf(terrorratio ~ condsize + (1|id),
               sigma ~ condsize + (1|id)),
            data = air_drag[air_drag$airdrag == "Airdrag",], family = gaussian())
Hypotheses = hypothesis(fit5,c("sigma_condsizeIncongruent < 0"))


fit6 <- brm(bf(xerrorratio ~ condsize + (1|id),
               sigma ~ condsize + (1|id)),
            data = air_drag[air_drag$airdrag == "Airdrag",], family = gaussian())
Hypotheses = hypothesis(fit6,c("sigma_condsizeIncongruent < 0"))

fit7 <- brm(bf(terrorratio ~ ball + (1|id),
               sigma ~ ball + (1|id)),
            data = air_drag, family = gaussian())
Hypotheses = hypothesis(fit7,c("sigma_ballTennis < 0"))

fit8 <- brm(bf(xerrorratio ~ ball + (1|id),
               sigma ~ ball + (1|id)),
            data = air_drag, family = gaussian())
Hypotheses = hypothesis(fit8,c("sigma_ballTennis < 0"))


######exploratory impact of air drag on
Hypotheses = hypothesis(fit1,c("sigma_airdragNoAirdrag > 0"))
Hypotheses = hypothesis(fit2,c("sigma_airdragNoAirdrag > 0"))



air_drag %>%
  group_by(ball,condsize) %>%
  mutate(Variability_Time = sd(terror_ratio),
         Variability_Space = sd(xerror_ratio)) %>%
  slice(1) %>%
  select(Variability_Time,Variability_Space)

ggplot(air_drag, aes(terror_ratio,xerror_ratio, color = ball)) +
  geom_point() +
  geom_smooth()




?set_prior


## Not run: 
## define priors
prior <- c(set_prior("normal(0,2)", class = "b"),
           set_prior("student_t(10,0,1)", class = "sigma"),
           set_prior("student_t(10,0,1)", class = "sd"))
?set_prior
## fit a linear mixed effects models
fit <- brm(time ~ age + sex + disease + (1 + age|patient),
           data = kidney, family = lognormal(),
           prior = prior, sample_prior = "yes", 
           control = list(adapt_delta = 0.95))

## perform two-sided hypothesis testing
(hyp1 <- hypothesis(fit, "sexfemale = age + diseasePKD"))
plot(hyp1)
hypothesis(fit, "exp(age) - 3 = 0", alpha = 0.01)
?hypothesis
## perform one-sided hypothesis testing
hypothesis(fit, "diseasePKD + diseaseGN - 3 < 0")

hypothesis(fit, "age < Intercept", 
           class = "sd", group  = "patient")

## test the amount of random intercept variance on all variance
h <- paste("sd_patient__Intercept^2 / (sd_patient__Intercept^2 +",
           "sd_patient__age^2 + sigma^2) = 0")
(hyp2 <- hypothesis(fit, h, class = NULL))
plot(hyp2)

## test more than one hypothesis at once
h <- c("diseaseGN = diseaseAN", "2 * diseaseGN - diseasePKD = 0")
(hyp3 <- hypothesis(fit, h))
plot(hyp3, ignore_prior = TRUE)

## compute hypotheses for all levels of a grouping factor
hypothesis(fit, "age = 0", scope = "coef", group = "patient")

## use the default method
dat <- as.data.frame(fit)
hypothesis(dat, "b_age > 0")

## End(Not run)
roh = 1.225
cd = 0.535
r = 0.12
D = 0.5 * roh * cd * pi * r^2
m = 0.6
vx = 3.5
vy = 6.86
g = 9.81

for (i in seq(0,1.4*0.575,0.01)){
  vt = (vx^2+vy^2)^0.5
  ax = -(D/m)*vt*vx
  ay = g - (D/m)*vt*vy
  
  vx = vx  + ax  * 0.01
  vy = vy + ay  * 0.01
  
  print(i)
}

