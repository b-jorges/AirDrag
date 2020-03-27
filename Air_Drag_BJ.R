libraries <- c("cowplot","tidyverse","grid","gridExtra","lme4","lmerTest",
               "ggpubr","apastats","rstatix", "rstan", "brms", "Hmisc")
# Where is this script?
invisible(lapply(libraries, function(x) {
  if(!require(x, character.only = T, quietly = T)) {
    install.packages(x)
    require(x, character.only = T)
  }
}
))
rm(libraries)

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

#set ggplot theme
theme_set(theme_cowplot())

#set working directory to the folder that contains this file
setwd(Where_Am_I())

#load necessary functions
source("Utilities/Funs.R")

#we pasted all data into one file, which is loaded here
collapsed <- read.table(file = "Data/All_Data_Air_Drag.txt", header = T) 

air_drag <-  collapsed %>% 
  select(trial,x_max,t_max,vx, #take all relevant data from the data file
         x_max_model,t_max_model, ball,cond_size,r,air_drag,label,
         random_x,rtime_timing,rtime_spatial,ball_x_spatial,
         id,TTC,visible) %>%
  #x_max is the point where the ball hit the table
  #t_max is the time at which the ball hit the table
  #x_max_model
  #x_max_model is the x position where participants should believe the ball hit the table, if they had a representation of air drag 
  #t_max_model is the moment participants should believe the ball hit the table, if they had a representation of air drag

  #ball indicates whether the ball had a tennis ball texture or a basketball texture
  #cond_size indicates whether the texture of the ball was congruent with its size and other air drag properties
  #r is the radius of the target
  #air_drag indicates whether airdrag was simulated in the first part of the trajectory or not
  #label is a categorical variable with 48 levels, one for each combination of horizontal velocity, 
       #air drag yes/no, ball type, congruency category and time-to-contact
  #random_x is the initial position of the ball the observers used to give their spatial response
  #rtime_timing is the moment from movement onset that observers pressed the button for the timing task
  #rtime_spatial is the time between appearance of the ball used for the spatial response until they 
      #pressed the button again to indicate they were satisfied with the position of the ball, i. e. the time
      #it took them to give their answer
  #ball_x_spatial is where participants indicated the ball hit the table
  #id are participant ids
  #TTC are the overall flight durations
  #visible denotes the time where the target became invisible


  mutate(
    #timing error with respect to the real time of impact:
    #we subtract 0.049s from the temporal responses because we have seen before that our projectors inrtoduce a delay of 0.049s      
         terror = rtime_timing-t_max-0.049, 
         #spatial error with respect to real point of impact:
         xerror = ball_x_spatial - x_max,
         #how long was the ball occluded:
         OccludedDuration = t_max-visible,
         #what percentage of the trajectory was the ball occluded:
         OccludedPercentage = visible/t_max,
         #for what length was the ball occluded in spatial terms:
         OccludedDistance = case_when(
           air_drag == 1 ~ x_max-(x_max/2+vx*0.8*t_max*(OccludedPercentage-0.5)), #vx is down to 80% of the original speed in air drag condition
           air_drag == 0 ~ x_max-(x_max/2+vx*t_max*(OccludedPercentage-0.5))),
         #the temporal error normalized by the duration of the occlusion:
         terrorratio = (OccludedDuration+terror) / OccludedDuration,
         #the spatial error normalized by the length of the occlusion:
         xerrorratio = (OccludedDistance+xerror) / OccludedDistance,
         
         #neater way of denoting variables:
         condsize = factor(cond_size,levels = c("cong","incongr"),
                            labels = c("Congruent","Incongruent")),
         ball = factor(ball,levels = c("tennis","basket"),
                       labels = c("Tennis","Basket")),
         airdrag = case_when(air_drag == 1 ~ "Airdrag",
                             air_drag == 0 ~ "NoAirdrag"))

#how many data points to we have before getting rid of outliers?
nAllTrials = length(air_drag$trial)


air_drag = air_drag %>% 
  group_by(id,label) %>%
  filter(terrorratio < 4 & terrorratio > 0.25,
         xerrorratio < 4 & xerrorratio > 0.25) %>%
  #exclude outliers that are more than 2.5 standard deviations above or below the mean of the corresponding condition/id
  filter(terrorratio > mean(terrorratio)-2.5*sd(terrorratio) & terrorratio < mean(terrorratio)+2.5*sd(terrorratio),
         xerrorratio > mean(xerrorratio)-2.5*sd(xerrorratio) & xerrorratio < mean(xerrorratio)+2.5*sd(xerrorratio)) %>%
  #trim
  filter(trim(terror, filter = T)) %>%
  filter(trim(xerror, filter = T))

#number of excluded trials
nAllTrials - length(air_drag$trial)


####################################################################
################Confirmatory Analyses###############################
####################################################################

##############Hypothesis 1: Do Humans have a Representation of Air Drag?
###Hypothesis 1: Temporal Error
#quick visualization:
ggplot(air_drag, aes(airdrag,terrorratio,color=airdrag)) +
  geom_violin() +
  geom_hline(yintercept=1) +
  geom_boxplot()

#Are temporal error different between airdrag present/absent?
H1_Temporal <- lmer(terrorratio ~ airdrag + (1|id), 
                         data = air_drag)
H1_Temporal_Null <- lmer(terrorratio ~ (1|id), 
                                   data = air_drag)
anova(H1_Temporal,H1_Temporal_Null)
summary(H1_Temporal)

#are temporal errors in either condition centered around 1?
H1_Temporal_Intercept1 <- lmer(terrorratio-1 ~ 0 + (1|id), 
                         data = air_drag[air_drag$airdrag == "Airdrag",])
H1_Temporal_Intercept1_Null <- lmer(terrorratio-1 ~ (1|id), 
                         data = air_drag[air_drag$airdrag == "Airdrag",])
anova(H1_Temporal_Intercept1,H1_Temporal_Intercept1_Null)
summary(H1_Temporal_Intercept1)
#not for the Airdrag: Present condition

H1_Temporal_Intercept2 <- lmer(terrorratio-1 ~ 0 + (1|id), 
                         data = air_drag[air_drag$airdrag == "NoAirdrag",])
H1_Temporal_Intercept2_Null <- lmer(terrorratio-1 ~ (1|id), 
                                   data = air_drag[air_drag$airdrag == "NoAirdrag",])
anova(H1_Temporal_Intercept2,H1_Temporal_Intercept2_Null)
summary(H1_Temporal_Intercept2)
#not for the Airdrag: Absent condition either



###Hypothesis 1b: Spatial Error
#quick visualization
ggplot(air_drag, aes(airdrag,xerrorratio,color = airdrag)) +
  geom_violin() +
  geom_hline(yintercept=1) +
  geom_boxplot()

#Are temporal error different between airdrag present/absent?
H1_Spatial <- lmer(xerrorratio ~ airdrag + (1|id), 
                         data = air_drag)
H1_Spatial_Null <- lmer(xerrorratio ~ (1|id), 
                                   data = air_drag)
anova(H1_Spatial,H1_Spatial_Null)
summary(H1_Spatial)

#are temporal errors in either condition centered around 1?
H1_Spatial <- lmer(xerrorratio-1 ~ 0 + (1|id), 
                         data = air_drag[air_drag$airdrag == "Airdrag",])
H1_Spatial_Null <- lmer(xerrorratio-1 ~ (1|id), 
                                   data = air_drag[air_drag$airdrag == "Airdrag",])
anova(H1_Spatial,H1_Spatial_Null)
summary(H1_Spatial)

H1_Spatial_Intercept1 <- lmer(xerrorratio-1 ~ 0 + (1|id), 
                         data = air_drag[air_drag$airdrag == "NoAirdrag",])
H1_Spatial_Intercept1_Null <- lmer(xerrorratio-1 ~ (1|id), 
                                   data = air_drag[air_drag$airdrag == "NoAirdrag",])
anova(H1_Spatial_Intercept1,H1_Spatial_Intercept1_Null)
summary(H1_Spatial_Intercept1)


#Hypothesis 1: Bayesian Linear Mixed Models
fit1 <- brm(bf(terrorratio ~ airdrag + (1|id),
               sigma ~ airdrag + (1|id)),
            data = air_drag, family = gaussian())

hypothesis(fit1,c("abs(Intercept-1) < abs(Intercept+airdragNoAirdrag-1)"))


fit2 <- brm(bf(xerrorratio ~ airdrag + (1|id),
               sigma ~ airdrag + (1|id)),
            data = air_drag, family = gaussian())

hypothesis(fit2,c("abs(Intercept-1) < abs(Intercept+airdragNoAirdrag-1)"))



##############Hypothesis 2: Does the texture of the ball have any impact?
###Hypothesis 2: Temporal Error
#quick visualization:
ggplot(air_drag, aes(condsize,terrorratio, color = ball)) +
  geom_violin() +
  geom_boxplot()

H2_Time <- lmer(terrorratio ~ ball*condsize + (1|id), 
                         data = air_drag[air_drag$airdrag == "Airdrag",])
H2_Time_Null <- lmer(terrorratio ~ ball + condsize + (1|id), 
                                      data = air_drag[air_drag$airdrag == "Airdrag",])
anova(H2_Time,H2_Time_Null)
summary(H2_Time)
#no difference

##space
ggplot(air_drag, aes(ball,xerrorratio, color = condsize)) +
  geom_violin() +
  geom_boxplot()

H2_Space <- lmer(xerrorratio ~ ball*condsize + (1|id), 
                                       data = air_drag[air_drag$airdrag == "Airdrag",])
H2_Space_Null <- lmer(xerrorratio ~ ball + condsize + (1|id), 
                                       data = air_drag[air_drag$airdrag == "Airdrag",])
anova(H2_Space,H2_Space_Null)
summary(H2_Space)
#yes difference


####################################################################
################Confirmatory Analyses###############################
####################################################################
###Is precision lower when no air drag is presented in the visible part of the trajectory?

hypothesis(fit1,c("sigma_airdragAirdrag < 0"))

hypothesis(fit2,c("sigma_airdragAirdrag < 0"))


###Does variability in responses differ between congruent and incongruent trials?
fit3 <- brm(bf(terrorratio ~ condsize + (1|id),
               sigma ~ condsize + (1|id)),
            data = air_drag[air_drag$airdrag == "Airdrag",], family = gaussian())
Hypotheses = hypothesis(fit3,c("sigma_condsizeIncongruent < 0"))


fit4 <- brm(bf(xerrorratio ~ condsize + (1|id),
               sigma ~ condsize + (1|id)),
            data = air_drag[air_drag$airdrag == "Airdrag",], family = gaussian())
Hypotheses = hypothesis(fit4,c("sigma_condsizeIncongruent < 0"))


###Do different ballsizes elicit differences in conditions?
#variability
fit5 <- brm(bf(terrorratio ~ ball + (1|id),
               sigma ~ ball + (1|id)),
            data = air_drag, family = gaussian())
Hypotheses = hypothesis(fit5,c("sigma_ballTennis < 0"))

fit6 <- brm(bf(xerrorratio ~ ball + (1|id),
               sigma ~ ball + (1|id)),
            data = air_drag, family = gaussian())
Hypotheses = hypothesis(fit6,c("sigma_ballTennis < 0"))

#biases
Expl_Ballsize_Time_Bias <- lmer(terrorratio ~ ball + (1|id), 
                             data = air_drag)
Expl_Ballsize_Bias_Time_Null <- lmer(terrorratio ~ (1|id), 
                             data = air_drag)
anova(Expl_Ballsize_Time_Bias,Expl_Ballsize_Bias_Time_Null)
summary(Expl_Ballsize_Time_Bias)


Expl_Ballsize_Space_Bias <- lmer(xerrorratio ~ ball + (1|id), 
                             data = air_drag)
Expl_Ballsize_Space_Bias_Null <- lmer(xerrorratio ~ (1|id), 
                             data = air_drag)
anova(Expl_Ballsize_Space_Bias,Expl_Ballsize_Space_Bias_Null)
summary(Expl_Ballsize_Space_Bias)


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


###Relationship between variability and bias
air_drag_sum <- air_drag %>% #make new data frame with means for the following values
  group_by(TTC,vx,id,air_drag,ball,cond_size) %>%
  mutate(sd_timing = sd(terror),
         sd_spatial = sd(xerror),
         visible = mean(visible),
         terrorratio = mean(terrorratio),
         xerrorratio = mean(xerrorratio),
         SDratio_t = sd_timing/mean(OccludedDuration),
         SDratio_x = sd_spatial/mean(OccludedDistance),
         t_max = mean(t_max),
         x_max = mean(x_max),
         terror = mean(terror),
         xerror= mean(xerror)) 

air_drag_VariabilityvsBias <- air_drag_sum %>%
  group_by(TTC,vx,id,air_drag,ball,cond_size) %>%
  slice(1)


#Overall Variability (between and within together)
Expl_BiasVsPrecision_Overall_Time <- lm(terrorratio ~ SDratio_t, 
                             data = air_drag_VariabilityvsBias)
summary(Expl_BiasVsPrecision_Overall_Time)

Expl_BiasVsPrecision_Overall_Space <- lm(xerrorratio ~ SDratio_x, 
                           data = air_drag_VariabilityvsBias)
summary(Expl_BiasVsPrecision_Overall_Space)


#within Variability
Expl_BiasVsPrecision_Within_Time <- lmer(terrorratio ~ SDratio_t + (1|id), 
                             data = air_drag_VariabilityvsBias)
Expl_BiasVsPrecision_Within_Time_Null <- lmer(terrorratio ~  (1|id), 
                             data = air_drag_VariabilityvsBias)
anova(Expl_BiasVsPrecision_Within_Time,Expl_BiasVsPrecision_Within_Time_Null)
summary(Expl_BiasVsPrecision_Within_Time)

Expl_BiasVsPrecision_Within_Space <- lmer(xerrorratio ~ SDratio_x + (1|id), 
                             data = air_drag_VariabilityvsBias)
Expl_BiasVsPrecision_Within_Space_Null <- lmer(terrorratio ~  (1|id), 
                             data = air_drag_VariabilityvsBias)
anova(Expl_BiasVsPrecision_Within_Space,Expl_BiasVsPrecision_Within_Space_Null)
summary(Expl_BiasVsPrecision_Within_Space)


#between Variability
Expl_BiasVsPrecision_Between_Time <- lmer(terrorratio ~ SDratio_t + (1|label), 
                             data = air_drag_VariabilityvsBias)
Expl_BiasVsPrecision_Between_Time_Null <- lmer(terrorratio ~  (1|label), 
                             data = air_drag_VariabilityvsBias) ###singular fit
anova(Expl_BiasVsPrecision_Between_Time,Expl_BiasVsPrecision_Between_Time_Null)
coef(H1_Spatial_NullModel)

Expl_BiasVsPrecision_Between_Space <- lmer(xerrorratio ~ SDratio_x + (1|label), 
                             data = air_drag_VariabilityvsBias) ###singular fit
Expl_BiasVsPrecision_Between_Space_Null <- lmer(xerrorratio ~  (1|label), 
                             data = air_drag_VariabilityvsBias) ###singular fit
anova(Expl_BiasVsPrecision_Between_Space,Expl_BiasVsPrecision_Between_Space_Null)
summary(Expl_BiasVsPrecision_Between_Space)

fit7 <- brm(bf(terrorratio ~ SDratio_t + (1|label),
               sigma ~ SDratio_t + (1|label)),
            data = air_drag_VariabilityvsBias, family = gaussian())
hypothesis(fit7, "SDratio_t > 0")

fit8 <- brm(bf(xerrorratio ~ SDratio_x + (1|label),
               sigma ~ SDratio_x + (1|label)),
            data = air_drag_VariabilityvsBias, family = gaussian())
hypothesis(fit8, "SDratio_t > 0")



####################################################################
############################Plots###################################
####################################################################


Exp_Conds <- expand.grid(Xi = 0,
                         Yi = 0,
                         vx = c(3,3.5),
                         vy = c(1,1.2,1.4)*9.807/2,
                         ball = c("tennis","basket"),
                         cond_size = c("cong","incongr"),
                         air_drag = c(1,0), rho = 1.225,
                         G = 9.807) %>%
  mutate(m = ifelse(ball == "tennis",
                    ifelse(cond_size == "cong",0.06,0.6),
                    ifelse(ball == "basket",
                           ifelse(cond_size == "cong",0.6,0.06),NA)),
         r = ifelse(ball == "tennis",
                    ifelse(cond_size == "cong",0.033,0.12),
                    ifelse(ball == "basket",
                           ifelse(cond_size == "cong",0.12,0.033),NA)),
         TTC = vy * 2 /9.807, 
         v = sqrt(vx^2+vy^2),
         cd = ifelse(air_drag == 0, 0, 0.535),
         id_TTC = paste0(ifelse(air_drag == 1, "Air_Drag_","Gravity_"),TTC),
         id_cond = paste0(id_TTC,"_",ball,"_",cond_size),
         id_exp = paste0(id_cond,"_vx_",vx), # Condition long identifier for experiment
         label = 1:n()) # Condition identifier for experiment

# =============================================================================
# Experimental conditions with max values
# =============================================================================
Exp_Conds <- Exp_Conds %>%
  group_by(air_drag,ball,cond_size,id_cond,m,r,rho,cd,TTC,id_TTC,v,G,vy,vx,Xi,Yi,id_exp,label) %>%
  do(xy_drag(vh = .$vx, vv = .$vy,C = .$cd,
             rho = .$rho,m = .$m,radius = .$r,
             dt = 0.001,g = .$G,vectors = F)) %>%
  mutate(x_max = max(x),
            y_max = max(y),
            t_max = max(t))


# =============================================================================
# Continuous
# =============================================================================
Continuous <- Exp_Conds %>%
  group_by(air_drag,ball,cond_size,id_cond,m,r,rho,cd,TTC,id_TTC,v,G,vy,vx,Xi,Yi,id_exp,label,t_max) %>%
  do(xy_drag_model(vh = vx, vv = vy,C = cd,
                   rho = rho,m = m,radius = r,
                   dt = 0.001,g = G,vectors = T))

# =============================================================================
# See experimental conditions with ggplot
# =============================================================================
# READ ME: if you want to see the proper figure download "setup.png" at root folder
# =============================================================================
try(image_setup <-  png::readPNG("Setup.png"))

Temp_cont <- Continuous %>%
  group_by(id_exp) %>%
  mutate(cond_size = factor(cond_size,levels = c("cong","incongr"),labels = c("Congruent","Incongruent")),
         ball = factor(ball,levels = c("tennis","basket"),labels = c("Tennis","Basket")),
         air_drag = factor(air_drag,levels = c(0,1),labels = c("Gravity","Gravity + Air Drag")),
         id_TTC = factor(id_TTC, levels = c("Air_Drag_1","Air_Drag_1.2", "Air_Drag_1.4",
                                            "Gravity_1", "Gravity_1.2","Gravity_1.4"),
                         labels = c("G + AD @ iTTC:1 (s)", "G + AD @ iTTC:1.2 (s)", "G + AD @ iTTC:1.4 (s)",
                                    "G @ iTTC:1 (s)", "G @ iTTC:1.2 (s)", "G @ iTTC:1.4 (s)")))
Temp_Exp_Conds <- Exp_Conds %>%
  group_by(id_exp) %>%
  mutate(cond_size = factor(cond_size,levels = c("cong","incongr"),labels = c("Congruent","Incongruent")),
         ball = factor(ball,levels = c("tennis","basket"),labels = c("Tennis","Basket")),
         air_drag = factor(air_drag,levels = c(0,1),labels = c("Gravity","Gravity + Air Drag")),
         id_TTC = factor(id_TTC, levels = c("Air_Drag_1","Air_Drag_1.2", "Air_Drag_1.4",
                                            "Gravity_1", "Gravity_1.2","Gravity_1.4"),
                         labels = c("G + AD @ iTTC:1 (s)", "G + AD @ iTTC:1.2 (s)", "G + AD @ iTTC:1.4 (s)",
                                    "G @ iTTC:1 (s)", "G @ iTTC:1.2 (s)", "G @ iTTC:1.4 (s)")))



if (exists("image_setup")){
  Setup_Plot <- plot_grid(
    ggplot(Temp_cont %>% filter(ball == "Basket"),aes(x,y+0.7, group = label, color = v)) +
      annotation_custom(rasterGrob(image_setup,
                                   width = unit(1,"npc"),
                                   height = unit(1,"npc")),
                        -Inf, Inf, -Inf, Inf) +
      geom_line() +
      coord_cartesian(ylim=c(0,max(Temp_cont$y)+0.7),xlim=c(-0.25,5.5)) +
      geom_line(data = Temp_cont %>%
                  filter(round(t,3) > round(t_max*0.3,3), round(t,3) < round(t_max*0.6,3),
                         ball == "Basket"), aes(x,y+0.7), color = "red") + 
      labs(x = "x (m)",
           y = "y (m)",
           color = "Condition") + 
      ggtitle(label =  "Experimental conditions",
              subtitle = "Red section indicates random oclusion range") +
      guides(color=FALSE)+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()),
    plot_shared_legend(
      ggplot(Temp_Exp_Conds,aes(TTC,t_max,color=factor(air_drag),shape=factor(cond_size))) +
        geom_point() +
        facet_wrap(~ball)+
        labs(x = "TTC under no Air Drag",
             y = "Real TTC",
             color = "Acceleration Condition",
             shape = "Ball size") + 
        scale_color_discrete(), 
      
      ggplot(Temp_Exp_Conds,aes(TTC,x_max,color=factor(air_drag),shape=factor(cond_size))) +
        geom_point() +
        facet_wrap(~ball) +
        labs(x = "TTC under no Air Drag",
             y = expression(X[End])) + 
        scale_color_discrete()
    ),
    nrow=2,rel_heights = c(0.6,0.4)
  )
} else{
  Setup_Plot <- plot_grid(
    ggplot(Temp_cont %>% filter(ball == "Basket"),aes(x,y+0.7, group = label, color = v)) +
      geom_line() +
      coord_cartesian(ylim=c(0,max(Temp_cont$y)+0.7),xlim=c(-0.25,5.5)) +
      geom_line(data = Temp_cont %>% 
                  filter(round(t,3) > round(t_max*0.55,3), round(t,3) < round(t_max*0.6,3), ball == "Basket"),
                aes(x,y+0.7), color = "red") + 
      labs(x = "x (m)",
           y = "y (m)",
           color = "Condition") + 
      guides(color=FALSE)+
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank()),
    plot_shared_legend(
      ggplot(Temp_Exp_Conds,aes(TTC,t_max,color=factor(air_drag),shape=factor(cond_size))) +
        geom_point() +
        facet_wrap(~ball)+
        labs(x = "TTC under no Air Drag",
             y = "Real TTC",
             color = "Acceleration Condition",
             shape = "Ball size") + 
        scale_color_discrete(), 
      
      ggplot(Temp_Exp_Conds,aes(TTC,x_max,color=factor(air_drag),shape=factor(cond_size))) +
        geom_point() +
        facet_wrap(~ball) +
        labs(x = "TTC under no Air Drag",
             y = expression(X[End])) + 
        scale_color_discrete()
    ),
    nrow=2,rel_heights = c(0.6,0.4)
  )
}



rm(Temp_cont,Temp_Exp_Conds,image_setup)
Setup_Plot
ggsave("Figure1 Experimental Conditions.jpg", w = 12, h = 8)


######################Air drag###########################
Figure2a = ggplot(air_drag, aes(airdrag,terrorratio,color = factor(airdrag))) +
  geom_hline(linetype = 2, yintercept = 1) + 
  geom_jitter(alpha = 0.025, width = 0.1) +
  geom_flat_violin(size=1) + 
  scale_color_discrete(name = "") + 
  theme(legend.position = "none") + 
#  stat_pvalue_manual(stat_test_t_b, label = "p.adj",color = "color_p") +
  labs(x = NULL,
       y = "Timing Error ratio") 

Figure2b = ggplot(air_drag, aes(airdrag,xerrorratio,color = factor(airdrag))) +
  geom_hline(linetype = 2, yintercept = 1) + 
  geom_jitter(alpha = 0.025, width = 0.1) +
  geom_flat_violin(size=1) + 
  scale_color_discrete(name = "") + 
  theme(legend.position = "none") + 
  #  stat_pvalue_manual(stat_test_t_b, label = "p.adj",color = "color_p") +
  labs(x = NULL,
       y = "Spatial Error ratio") 

plot_grid(Figure2a,Figure2b, labels = "AUTO")
ggsave("Figure2.jpg", w = 12, h = 6)


######################Congruency###########################
air_drag = air_drag %>%
  mutate(BallXCongruency = 
           case_when(
             condsize == "Congruent" & ball == "Tennis" ~ "Tennis, Congruent",
             condsize == "Incongruent" & ball == "Tennis" ~ "Tennis, Incongruent",
             condsize == "Congruent" & ball == "Basket" ~ "Basket, Congruent",
             condsize == "Incongruent" & ball == "Basket" ~ "Basket, Incongruent",
           )
         )

Figure3a = ggplot(air_drag, aes(BallXCongruency,terrorratio,color = as.factor(BallXCongruency))) +
  geom_hline(linetype = 2, yintercept = 1) + 
  geom_jitter(alpha = 0.025, width = 0.1) +
  geom_flat_violin(size=1) + 
  scale_color_manual(name = "", values = c(Red, LightRed, BlauUB, LightBlauUB)) + 
  theme(legend.position = "") + 
  labs(x = NULL,
       y = "Timing Error ratio") 

Figure3b = ggplot(air_drag, aes(BallXCongruency,xerrorratio,color = as.factor(BallXCongruency))) +
  geom_hline(linetype = 2, yintercept = 1) + 
  geom_jitter(alpha = 0.025, width = 0.1) +
  geom_flat_violin(size=1) + 
  scale_color_manual(name = "", values = c(Red, LightRed, BlauUB, LightBlauUB)) + 
  theme(legend.position = "") + 
  labs(x = NULL,
       y = "Spatial Error ratio") 

plot_grid(Figure3a,Figure3b, nrow = 2, labels = "AUTO")
ggsave("Figure3.jpg", w = 12, h = 6)

######################Congruency###########################
Figure4a = ggplot(air_drag, aes(condsize,terrorratio,color = factor(condsize))) +
  geom_hline(linetype = 2, yintercept = 1) + 
  geom_jitter(alpha = 0.025, width = 0.1) +
  geom_flat_violin(size=1) + 
  scale_color_manual(name = "", values = c(Lila, Turquoise)) + 
  theme(legend.position = "none") + 
  #  stat_pvalue_manual(stat_test_t_b, label = "p.adj",color = "color_p") +
  labs(x = NULL,
       y = "Timing Error ratio") 

Figure4b = ggplot(air_drag, aes(condsize,xerrorratio,color = factor(condsize))) +
  geom_hline(linetype = 2, yintercept = 1) + 
  geom_jitter(alpha = 0.025, width = 0.1) +
  geom_flat_violin(size=1) + 
  scale_color_manual(name = "", values = c(Lila, Turquoise)) + 
  theme(legend.position = "none") + 
  #  stat_pvalue_manual(stat_test_t_b, label = "p.adj",color = "color_p") +
  labs(x = NULL,
       y = "Spatial Error ratio") 

plot_grid(Figure4a,Figure4b, labels = "AUTO")
ggsave("Figure4.jpg", w = 12, h = 6)

######################Ball Size###########################
Figure5a = ggplot(air_drag %>% mutate(BallSize = case_when(r == 0.033 ~ "0.033 m", r == 0.12 ~ "0.12 m")), 
                                      aes(BallSize,terrorratio,color = BallSize)) +
  geom_hline(linetype = 2, yintercept = 1) + 
  geom_jitter(alpha = 0.025, width = 0.1) +
  geom_flat_violin(size=1) + 
  scale_color_manual(name = "", values = c("green", "darkgreen")) + 
  theme(legend.position = "none") + 
  #  stat_pvalue_manual(stat_test_t_b, label = "p.adj",color = "color_p") +
  labs(x = NULL,
       y = "Timing Error ratio") 

Figure5b = ggplot(air_drag %>% mutate(BallSize = case_when(r == 0.033 ~ "0.033 m", r == 0.12 ~ "0.12 m")), 
                  aes(BallSize,xerrorratio,color = BallSize)) +
  geom_hline(linetype = 2, yintercept = 1) + 
  geom_jitter(alpha = 0.025, width = 0.1) +
  geom_flat_violin(size=1) + 
  scale_color_manual(name = "", values = c("green", "darkgreen")) + 
  theme(legend.position = "none") + 
  #  stat_pvalue_manual(stat_test_t_b, label = "p.adj",color = "color_p") +
  labs(x = NULL,
       y = "Spatial Error ratio") 

plot_grid(Figure5a,Figure5b, labels = "AUTO")
ggsave("Figure5.jpg", w = 12, h = 6)


# =============================================================================
# a) Timing: Variability ratio vs. error ratio 
# =============================================================================

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

Figure10a = ggplot(air_drag_sum,aes(terrorratio,SDratio_t,fill = id)) + 
  geom_point(alpha = 0.025, shape = 21) + 
  geom_point(data= air_drag_participant, size = 5, shape = 21) +
  labs(x = expression("Error Ratio (s)"),
       y = expression("Normalized SD (s)"),
       color = NULL) +
  guides(fill = F) +
  scale_fill_viridis_d()

# =============================================================================
# b) Spatial: Variability ratio vs. error ratio 
# =============================================================================

Figure10b = ggplot(air_drag_sum,aes(xerrorratio,SDratio_x,fill = id)) + 
  geom_point(alpha = 0.025, shape = 21) + 
  geom_point(data= air_drag_participant, size = 5, shape = 21) +
  labs(x = expression("Error Ratio (m)"),
       y = expression("Normalized SD (m)"),
       color = NULL) +
  guides(fill = F) +
  scale_fill_viridis_d() +
  coord_cartesian(ylim = c(0,0.4), xlim = c(0.4,1.4))

plot_grid(Figure10a,Figure10b, labels = "AUTO")
ggsave("Figure10.jpg", w = 12, h = 6)
# =============================================================================
# c) Timing variability ratio vs. Spatial variability ratio
# =============================================================================

ggplot(air_drag_sum,aes(SDratio_t,SDratio_x,fill = id)) + 
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
