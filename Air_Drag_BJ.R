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

air_drag <-  collapsed %>%
  select(trial,x_max,y_max,t_max,vx,
         x_max_model,t_max_model,y_max_model,
         G,ball,cond_size,r,air_drag,label,
         random_x,rtime_timing,rtime_spatial,ball_x_spatial,ball_x_timing,
         id,TTC,cond_size,visible) %>%
  # Raw differneces
  mutate(terror = rtime_timing-t_max-0.049,
         xerror = ball_x_spatial - x_max,
         # Differences with model modelling AD at oclusion
         xerror_model = ball_x_spatial - x_max_model,
         terror_model = rtime_timing-t_max_model-0.049,
         # Ratio Raw
         xerror_ratio = (x_max+xerror) / x_max,
         terror_ratio = (terror+t_max) / t_max,
         # Ratio Model
         xerror_ratio_model = xerror_model / x_max_model,
         terror_ratio_model = terror_model / t_max_model,
         cond_size = factor(cond_size,levels = c("cong","incongr"),
                            labels = c("Congruent","Incongruent")),
         ball = factor(ball,levels = c("tennis","basket"),
                       labels = c("Tennis","Basket")),
         air_drag = ordered(factor(air_drag),levels = c(0,1),
                            labels = c("Gravity","Gravity + Air_Drag"))
  ) %>%
  group_by(id,label) %>%
  filter(abs(terror) < 1,
         abs(xerror) < 2) %>%
  filter(trim(terror, filter = T)) %>%
  filter(trim(xerror, filter = T))


#####Look at airdrag and stuff
ggplot(air_drag, aes(air_drag,terror_ratio)) +
  geom_violin()
H1_Temporal_lmer <- lmer(terror_ratio ~ air_drag + (1|id), 
                         data = air_drag)
H1_Temporal_lmer_NullModel <- lmer(terror_ratio ~ (1|id), 
                         data = air_drag)
anova(H1_Temporal_lmer,H1_Temporal_lmer_NullModel)


ggplot(air_drag, aes(air_drag,xerror_ratio)) +
  geom_violin()
H1_Spatial_TestModel <- lmer(xerror_ratio ~ air_drag + (1|id), 
                         data = air_drag)
H1_Spatial_NullModel <- lmer(xerror_ratio ~ (1|id), 
                                   data = air_drag)
anova(H1_Spatial_TestModel,H1_Spatial_NullModel)



######Look at ball size and congruency
##time
ggplot(air_drag, aes(cond_size,terror_ratio, color = ball)) +
  geom_violin()
SizeCongruency_TestModel_Time <- lmer(terror_ratio ~ ball*cond_size + (1|id), 
                         data = air_drag)
SizeCongruency_NullModel_Time <- lmer(terror_ratio ~ (1|id), 
                                   data = air_drag)
anova(SizeCongruency_TestModel_Time,SizeCongruency_NullModel_Time)
#(lol)

##space
ggplot(air_drag, aes(cond_size,xerror_ratio, color = ball)) +
  geom_violin()
SizeCongruency_TestModel_Space <- lmer(xerror_ratio ~ ball*cond_size + (1|id), 
                                      data = air_drag)
SizeCongruency_NullModel_Space <- lmer(xerror_ratio ~ ball + cond_size + (1|id), 
                                      data = air_drag)
anova(SizeCongruency_TestModel_Space,SizeCongruency_NullModel_Space)
summary(SizeCongruency_TestModel_Space)
#interesting

air_drag %>%
  group_by(ball,cond_size) %>%
  mutate(Variability_Time = sd(terror_ratio),
         Variability_Space = sd(xerror_ratio)) %>%
  slice(1) %>%
  select(Variability_Time,Variability_Space)

ggplot(air_drag, aes(terror_ratio,xerror_ratio, color = ball)) +
  geom_point() +
  geom_smooth()
