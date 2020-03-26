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

require(rstan)
require(brms)
require(Hmisc)

#####optimize setup for Bayesian Linear Models (rstan/brms)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
Sys.setenv(LOCAL_CPPFLAGS = '-march=corei7')

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
         G,ball,cond_size,r,air_drag,label,
         random_x,rtime_timing,rtime_spatial,ball_x_spatial,ball_x_timing,
         id,TTC,visible) %>%
  # Raw differneces
  mutate(terror = rtime_timing-t_max-0.049,
         xerror = ball_x_spatial - x_max,
         OccludedDuration = t_max-visible,
         OccludedPercentage = visible/t_max,
         OccludedDistance = case_when(
           air_drag == 1 ~ x_max-(x_max/2+vx*0.8*t_max*(OccludedPercentage-0.5)), #vx is down to 80% of the original speed in air drag condition
           air_drag == 0 ~ x_max-(x_max/2+vx*t_max*(OccludedPercentage-0.5))),
         # Differences with model modelling AD at oclusion
         xerror_model = ball_x_spatial - x_max_model,
         terror_model = rtime_timing-t_max_model-0.049,
         # Ratio Raw
         xerrorratio = (OccludedDistance+xerror) / OccludedDistance,
         terrorratio = (OccludedDuration+terror) / OccludedDuration,
         # Ratio Model -
         xerror_ratio_model = xerror_model / x_max_model,
         terror_ratio_model = terror_model / t_max_model,
         condsize = factor(cond_size,levels = c("cong","incongr"),
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
         sd_space = sd(xerror),
         Median_terror = median(terror),
         Median_xerror = median(xerror),
         SD_Ratio_t = sd(terror)/(OccludedDuration + Median_terror),
         SD_Ratio_x = sd(terror)/(OccludedDistance + Median_xerror))

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

H1_Spatial_TestModel <- lmer(terrorratio ~ ball + (1|id), 
                             data = air_drag)
H1_Spatial_NullModel <- lmer(terrorratio ~ (1|id), 
                             data = air_drag)
anova(H1_Spatial_TestModel,H1_Spatial_NullModel)
summary(H1_Spatial_TestModel)


H1_Spatial_TestModel <- lmer(xerrorratio ~ ball + (1|id), 
                             data = air_drag)
H1_Spatial_NullModel <- lmer(xerrorratio ~ (1|id), 
                             data = air_drag)
anova(H1_Spatial_TestModel,H1_Spatial_NullModel)
summary(H1_Spatial_TestModel)



air_drag %>%
  group_by(ball,condsize) %>%
  mutate(Variability_Time = sd(terror_ratio),
         Variability_Space = sd(xerror_ratio)) %>%
  slice(1) %>%
  select(Variability_Time,Variability_Space)

ggplot(air_drag, aes(terror_ratio,xerror_ratio, color = ball)) +
  geom_point() +
  geom_smooth()



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







# =============================================================================
# Response variability, errors and conditions
# =============================================================================
air_drag_sum <- air_drag %>%
  group_by(TTC,vx,id,air_drag,ball,cond_size) %>%
  mutate(sd_timing = sd(terror),
         sd_spatial = sd(xerror),
         visible = mean(visible),
         terrorratio = mean(terrorratio),
         xerrorratio = mean(xerrorratio),
         SDratio_t = sd_timing/mean(t_max),
         SDratio_x = sd_spatial/mean(x_max),
         #xerror_t_ratio = mean(ball_x_spatial)/mean(x_max),
         t_max = mean(t_max),
         x_max = mean(x_max),
         terror = mean(terror),
         xerror= mean(xerror)) 

air_drag_participant <- air_drag_sum %>%
  group_by(id) %>%
  summarize_all(.funs = "mean")


air_drag_ratios <- air_drag_sum %>%
  ungroup() %>%
  select(-id,-cond_size,-ball) %>%
  summarize_all(.funs = "mean")


mean_ratio_sd_air_Drag <- air_drag_sum %>%
  group_by(id) %>%
  mutate(max_xerror = mean_cl_normal(SDratio_x)$ymax,
         min_xerror = mean_cl_normal(SDratio_x)$ymin,
         ratio_x = mean_cl_normal(SDratio_x)$y,
         max_terror = mean_cl_normal(SDratio_t)$ymax,
         min_terror = mean_cl_normal(SDratio_t)$ymin,
         ratio_t = mean_cl_normal(SDratio_t)$y)



air_drag_VariabilityvsBias <- air_drag_sum %>%
  group_by(TTC,vx,id,air_drag,ball,cond_size) %>%
  slice(1)

H1_Spatial_TestModel <- lmer(terrorratio ~ SDratio_t + (1|TTC), 
                             data = air_drag_VariabilityvsBias)
H1_Spatial_NullModel <- lmer(terrorratio ~ (1|TTC), 
                             data = air_drag_VariabilityvsBias)
anova(H1_Spatial_TestModel,H1_Spatial_NullModel)
summary(H1_Spatial_TestModel)
plot(H1_Spatial_TestModel)

# =============================================================================
# a) Timing: Variability ratio vs. error ratio 
# =============================================================================
ggplot(air_drag_sum,aes(terror_ratio,ratio_t,fill = id)) + 
  geom_abline(linetype = 2) +
  geom_point(alpha = 0.25, shape = 21) + 
  geom_point(data= air_drag_participant, size = 3, shape = 21) +
  stat_cor(data = air_drag_participant, aes(group = 0)) + 
  labs(x = expression(rt/t[max]),
       y = expression(sigma[t]/t[max]),
       color = NULL) + 
  guides(fill = F) +
  scale_fill_viridis_d()

# =============================================================================
# b) Spatial: Variability ratio vs. error ratio 
# =============================================================================

ggplot(air_drag_sum,aes(xerror_ratio,ratio_x,fill = id)) + 
  geom_abline(linetype = 2) +
  geom_point(alpha = 0.25, shape = 21) + 
  geom_point(data= air_drag_participant, size = 3, shape = 21) +
  stat_cor(data = air_drag_participant, aes(group = 0)) + 
  labs(x = expression(rx/x[max]),
       y = expression(sigma[x]/x[max]),
       color = NULL) + 
  guides(fill = F) +
  scale_fill_viridis_d() 

# =============================================================================
# c) Timing variability ratio vs. Spatial variability ratio
# =============================================================================

ggplot(air_drag_sum,aes(ratio_x,ratio_t,fill = id)) + 
  geom_abline(linetype = 2) + 
  geom_point(alpha = 0.25, shape = 21) + 
  geom_point(data= air_drag_participant, size = 3, shape = 21) +
  ggpubr::stat_cor(aes(group = 0))+ 
  labs(title = paste0("Spatial ratio: ", round( air_drag_ratios$ratio_x[1],3),
                      "\nTemporal ratio: ", round( air_drag_ratios$ratio_t[1],3)),
       x = expression(sigma[x]/x[max]),
       y = expression(sigma[t]/t[max]),
       color = NULL) + 
  scale_fill_viridis_d() +
  guides(fill = FALSE) +
  # geom_errorbarh(data = mean_ratio_sd_air_Drag, 
  #                aes(x = ratio_x,
  #                    xmax = max_xerror,   
  #                    xmin = min_xerror)) +   
  # geom_errorbar(data = mean_ratio_sd_air_Drag,aes(y = ratio_t,
  #                    ymax = max_terror,
  #                    ymin = min_terror)) +
  theme_minimal(12)


# =============================================================================
# c) Terror ratio vs. Xerror ratio
# =============================================================================
ggplot(air_drag_sum  ,aes(xerror_ratio,terror_ratio,fill = id)) + 
  geom_abline(linetype = 2) + 
  geom_point(alpha = 0.25, shape = 21) + 
  geom_point(data= air_drag_participant, size = 3, shape = 21) +
  ggpubr::stat_cor(aes(group = 0))+ 
  labs(x = expression(r[x]/x[max]),
       y = expression(r[t]/t[max]),
       color = NULL) + 
  scale_fill_viridis_d() +
  #guides(fill = FALSE) +
  #facet_wrap(~id) +
  theme_minimal(12)
