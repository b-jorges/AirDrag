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
set.seed(777)

######################Colorblind safe pallete###########################
# The palette with grey:
cbp1 <- c("#999999", "#E69F00", "#56B4E9", "#009E73",
          "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
########################################################################
#####optimize setup for Bayesian Linear Models (rstan/brms)
########################################################################
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
             dt = 0.001,g = .$G,vectors = T)) %>%
  group_by(label) %>% 
  mutate(x_max = max(x),
         y_max = max(y),
         t_max = max(t),
         occlusion = ifelse(t > t_max * 0.55 & t < t_max * 0.6,"Occlusion period",
                            ifelse(t < t_max * 0.55, "Visible","Occluded"))) 



# =============================================================================
# See experimental conditions with ggplot
# =============================================================================
# READ ME: if you want to see the proper figure download "setup.png" at root folder
# =============================================================================
try(image_setup <-  png::readPNG("Setup.png"))

Temp_cont <- Exp_Conds %>%
  group_by(id_exp) %>%
  mutate(cond_size = factor(cond_size,levels = c("cong","incongr"),labels = c("Congruent","Incongruent")),
         ball = factor(ball,levels = c("tennis","basket"),labels = c("Tennis","Basket")),
         air_drag = factor(air_drag,levels = c(0,1),labels = c("Gravity","Gravity + Air Drag")),
         id_TTC = factor(id_TTC, levels = c("Air_Drag_1","Air_Drag_1.2", "Air_Drag_1.4",
                                            "Gravity_1", "Gravity_1.2","Gravity_1.4"),
                         labels = c("G + AD @ iTTC:1 (s)", "G + AD @ iTTC:1.2 (s)", "G + AD @ iTTC:1.4 (s)",
                                    "G @ iTTC:1 (s)", "G @ iTTC:1.2 (s)", "G @ iTTC:1.4 (s)")))
Temp_Exp_Conds_max <- Exp_Conds %>%
  group_by(id_exp) %>%
  mutate(cond_size = factor(cond_size,levels = c("cong","incongr"),labels = c("Congruent","Incongruent")),
         ball = factor(ball,levels = c("tennis","basket"),labels = c("Tennis","Basket")),
         ball_condition = ifelse(ball == "Tennis" & cond_size == "Congruent", "Tennis Properties",
                                 ifelse(ball == "Tennis" & cond_size == "Incongruent", "Basket Properties",
                                        ifelse(ball == "Basket" & cond_size == "Congruent", "Basket Properties", "Tennis Properties"))),
         air_drag = factor(air_drag,levels = c(0,1),labels = c("Gravity","Gravity + Air Drag")),
         id_TTC = factor(id_TTC, levels = c("Air_Drag_1","Air_Drag_1.2", "Air_Drag_1.4",
                                            "Gravity_1", "Gravity_1.2","Gravity_1.4"),
                         labels = c("G + AD @ iTTC:1 (s)", "G + AD @ iTTC:1.2 (s)", "G + AD @ iTTC:1.4 (s)",
                                    "G @ iTTC:1 (s)", "G @ iTTC:1.2 (s)", "G @ iTTC:1.4 (s)"))) %>%
  filter(round(t,3) == 0.001)

Figure1a <- ggplot(Temp_cont %>% filter(ball == "Basket"),aes(x,y+0.7, group = interaction(id_exp,label),color = occlusion)) +
  annotation_custom(rasterGrob(image_setup,
                               width = unit(1,"npc"),
                               height = unit(1,"npc")),
                    -Inf, Inf, -Inf, Inf) +
  geom_point(size = 1.25) +
  coord_cartesian(ylim=c(0,max(Temp_cont$y)+0.7),xlim=c(-0.25,5.5)) +
  labs(x = "x (m)",
       y = "y (m)") + 
  ggtitle("A") + 
  scale_color_manual(name = NULL, values = cbp1[1:3]) + 
  #guides(color=FALSE) + 
  theme_cowplot(15) + 
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=2)))


Figure1b <- Exp_Conds %>%
  filter(ball == "basket") %>%
  select(-ay) %>%
  #mutate(id = factor(id, unique(id), labels = c("Horizontal acceleration","Vertical acceleration"))) %>% 
  ggplot(.,aes(t,ax, color = occlusion, group = label)) +
  geom_point(size = 0.75) + 
  labs(x = "t (s)",
       y = expression(Horizontal~acceleration~(frac(m,s^2)))) + 
  scale_color_manual(name = NULL, values = cbp1[1:3]) + 
  ggtitle("B") + 
  theme_cowplot(15) + 
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=2)))

Figure1c <- Exp_Conds %>%
  filter(ball == "basket") %>%
  select(-ax) %>%
  #mutate(id = factor(id, unique(id), labels = c("Horizontal acceleration","Vertical acceleration"))) %>% 
  ggplot(.,aes(t,ay, color = occlusion, group = label)) +
  geom_point(size = 0.75) + 
  labs(x = "t (s)",
       y = expression(Vertical~acceleration~(frac(m,s^2)))) + 
  scale_color_manual(name = NULL, values = cbp1[1:3]) + 
  ggtitle("C") + 
  theme_cowplot(15) + 
  theme(legend.position = "none")

Figure1d <- ggplot(Temp_Exp_Conds_max ,aes(TTC,t_max,color=factor(air_drag),shape=factor(ball_condition))) +
  geom_point(size = 3) +
  labs(x = "TTC under no Air Drag (s)",
       y = "Real TTC (s)",
       color = "Air drag conditions:",
       shape = "Ball size:") + 
  scale_color_manual(values = c("red","royalblue1")) + 
  ggtitle("D") +
  theme_cowplot(15)

Figure1e <- ggplot(Temp_Exp_Conds_max,aes(TTC,x_max,color=factor(air_drag),shape=factor(ball_condition))) +
  geom_point(size = 3) +
  #facet_wrap(~ball) +
  labs(x = "TTC under no Air Drag (s)",
       y = "Horizontal point of impact (m)") + 
  scale_color_manual(values = c("red","royalblue1")) + 
  ggtitle("E") + 
  theme_cowplot(15)



Figure1 <- plot_grid(plot_shared_legend(Figure1a),
                     plot_shared_legend(Figure1b,Figure1c,ncol = 2,position = "bottom"),
                     plot_shared_legend(Figure1d,Figure1e,ncol = 2,position = "bottom"),
                     nrow = 3)

ggsave(Figure1, filename = "Figure1 Experimental Conditions.jpg", w = 12, h = 12)



# =============================================================================
# Experimental data
# =============================================================================
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
  #spatial task: exclude trials where targets were placed at one fourth of the occluded distance, or four times the occluded distance
  filter(xerrorratio < 4 & xerrorratio > 0.25) %>% 
  #timing task: exclude trials where subjects pressed the button less than 100ms after disappearance or more than 3s after disappearance 
  #(the minimum occluded duration was 380 and the maximum was 640ms)
  filter(rtime_timing > visible + 0.1 & rtime_timing < 3) %>% 
  group_by(id)

#determine if more than 10% of trials are removed by the previous step for any of the participants
outliers <- air_drag %>% 
  group_by(id) %>% 
  summarise(outliers =n()/960) %>%  
  filter(outliers > 0.9)
1-mean(outliers$outliers)

#remove s02 and s04, who did not meet this criterion
air_drag = air_drag %>% 
  mutate(remove = ifelse(n()/960 < 0.9, 1, 0)) %>% 
  filter(remove == 0)

#outlier removal based on quartiles (removes all trials beyond the "antennae" of a boxplot)
air_drag = air_drag %>% 
  group_by(id,label) %>%
  filter(trim(terror, filter = T)) %>%
  filter(trim(xerror, filter = T))

#how many trials did we lose?
nRemoveN = 20 -  length(as.vector(unique(air_drag$id))) 
nRemove = (20 - nRemoveN) * 960 - nrow(air_drag)
percRemove = nRemove / ((20 - nRemoveN) * 960) *100

####################################################################
################Confirmatory Analyses###############################
####################################################################
#Are temporal error different between airdrag present/absent?
####################################################################
###Hypothesis 1a: Temporal Error
####################################################################

H1_Temporal <- lmer(terrorratio ~ airdrag + (1|id), 
                    data = air_drag)
summary(H1_Temporal)
# confint(H1_Temporal, method = "boot")

H1_Temporal_Null <- lmer(terrorratio ~ (1|id), 
                         data = air_drag)
summary(H1_Temporal_Null)

### Compare models
anova(H1_Temporal,H1_Temporal_Null)
####################################################################
# Are temporal errors in either condition centered around 1?
####################################################################
# For Air Drag: 
H1_Temporal_Intercept1 <- lmer(terrorratio-1 ~ 0 + (1|id), 
                               data = air_drag[air_drag$airdrag == "Airdrag",])
summary(H1_Temporal_Intercept1)

H1_Temporal_Intercept1_Null <- lmer(terrorratio-1 ~ (1|id), 
                                    data = air_drag[air_drag$airdrag == "Airdrag",])
summary(H1_Temporal_Intercept1_Null)
set.seed(777)
confint(H1_Temporal_Intercept1_Null,method = "boot")

### Compare models
anova(H1_Temporal_Intercept1,H1_Temporal_Intercept1_Null)
####################################################################
# For Gravity: 

H1_Temporal_Intercept2 <- lmer(terrorratio-1 ~ 0 + (1|id), 
                               data = air_drag[air_drag$airdrag == "NoAirdrag",])
summary(H1_Temporal_Intercept2)

H1_Temporal_Intercept2_Null <- lmer(terrorratio-1 ~ (1|id), 
                                    data = air_drag[air_drag$airdrag == "NoAirdrag",])
set.seed(777)
confint(H1_Temporal_Intercept2_Null,method = "boot")
summary(H1_Temporal_Intercept2_Null)

### Compare models
anova(H1_Temporal_Intercept2,H1_Temporal_Intercept2_Null)



##################################################################
#Hypothesis 1: Bayesian Linear Mixed Models
####################################################################
#fit1 <- brm(bf(terrorratio ~ airdrag + (1|id), #the fitting of this model can take more than an hour
#                sigma ~ airdrag + (1|id)),
#             data = air_drag, family = gaussian())
#save(fit1,file = paste0(Where_Am_I(),"/fit1.RData")) 

load(paste0(Where_Am_I(),"/fit1.RData")) #you can load the fitted object by using this line, to avoid having to fit the model
fit1
####################################################################
### Translated: Is air_drag closer to zero than non_air_drag?
####################################################################
view(fit1)
h1a <- hypothesis(fit1,c("abs(Intercept-1) < abs(Intercept+airdragNoAirdrag-1)"))
h1a$hypothesis

####################################################################
#Are spatial error different between airdrag present/absent?
####################################################################
###Hypothesis 1b: spatial Error
####################################################################

H1_Spatial <- lmer(xerrorratio ~ airdrag + (1|id), 
                   data = air_drag)
summary(H1_Spatial)

H1_Spatial_Null <- lmer(xerrorratio ~ (1|id), 
                        data = air_drag)
summary(H1_Spatial_Null)

### Compare models
anova(H1_Spatial,H1_Spatial_Null)
####################################################################
# Are spatial errors in either condition centered around 1?
####################################################################
# For Air Drag: 
H1_Spatial_Intercept1 <- lmer(xerrorratio-1 ~ 0 + (1|id), 
                              data = air_drag[air_drag$airdrag == "Airdrag",])
summary(H1_Spatial_Intercept1)


H1_Spatial_Intercept1_Null <- lmer(xerrorratio-1 ~ (1|id), 
                                   data = air_drag[air_drag$airdrag == "Airdrag",])
summary(H1_Spatial_Intercept1_Null)
set.seed(777)
confint(H1_Spatial_Intercept1_Null,method = "boot")

### Compare models
anova(H1_Spatial_Intercept1,H1_Spatial_Intercept1_Null)
####################################################################
# For Gravity: 

H1_Spatial_Intercept2 <- lmer(xerrorratio-1 ~ 0 + (1|id), 
                              data = air_drag[air_drag$airdrag == "NoAirdrag",])
summary(H1_Spatial_Intercept2)

H1_Spatial_Intercept2_Null <- lmer(xerrorratio-1 ~ (1|id), 
                                   data = air_drag[air_drag$airdrag == "NoAirdrag",])
summary(H1_Spatial_Intercept2_Null)
set.seed(777)
confint(H1_Spatial_Intercept2_Null,method = "boot")
?confint
### Compare models
anova(H1_Spatial_Intercept2,H1_Spatial_Intercept2_Null)

####################################################################
#Hypothesis 1: Bayesian Linear Mixed Models
####################################################################
#fit2 <- brm(bf(xerrorratio ~ airdrag + (1|id), #the fitting of this model can take more than an hour
#                sigma ~ airdrag + (1|id)),
#             data = air_drag, family = gaussian())
#save(fit2,file = paste0(Where_Am_I(),"/fit2.RData")) #you can load the fitted object by using this line, to avoid having to fit the model

load(paste0(Where_Am_I(),"/fit2.RData"))
fit2
h1b <- hypothesis(fit2,c("abs(Intercept-1) < abs(Intercept+airdragNoAirdrag-1)"))
h1b$hypothesis

######################Plots Hypothesis 1 ###########################
air_drag = air_drag %>% 
  group_by(airdrag) %>%
  mutate(xerrorratio_mean_by_airdrag = mean(xerrorratio),
         xerrorratio_SE_by_airdrag = sd(xerrorratio)/(length(xerrorratio))^0.5,
         terrorratio_mean_by_airdrag = mean(terrorratio),
         terrorratio_SE_by_airdrag = sd(terrorratio)/(length(terrorratio))^0.5)

(Figure2a = ggplot(air_drag, aes(airdrag,terrorratio,color = factor(airdrag))) +
   geom_hline(linetype = 2, yintercept = 1) + 
   geom_jitter(alpha = 0.025, width = 0.1) +
   geom_flat_violin(size=1) + 
   stat_summary(fun = "median", geom = "point",size = 6, aes(group=id), shape = 95) +
   stat_summary(fun = "median", geom = "point",size = 4, aes(group=0), shape = 16, color = "black") +
   scale_color_manual(name = "", values = c("grey","gold3")) + 
   theme(legend.position = "none") + 
   #facet_wrap(TTC~id) +
   labs(x = NULL,
        y = "Timing Error ratio"))

Figure2a_inset = ggplot(air_drag, aes(x = airdrag,
                                      y = terrorratio_mean_by_airdrag,
                                      ymin = terrorratio_mean_by_airdrag-terrorratio_SE_by_airdrag,
                                      ymax = terrorratio_mean_by_airdrag+terrorratio_SE_by_airdrag,
                                      color = factor(airdrag))) +
  geom_point(size = 3) +
  geom_errorbar(width=0.2, size = 1.5) +
  scale_color_manual(name = "", 
                     values = c("grey","gold3")) +
  scale_x_discrete(name ="", 
                   labels = c("Airdrag" = "AD", "NoAirdrag" = "No AD")) +
  scale_y_continuous(breaks=c(1.01,1.02,1.03)) +
  theme(legend.position = "", 
        axis.text=element_text(size=12),
        axis.title=element_text(size=12)) +
  xlab("") +
  ylab("Error Ratio")

Figure2a_Complete <-
  ggdraw() +
  draw_plot(Figure2a) +
  draw_plot(Figure2a_inset, x = 0.38, y = .6, width = .3, height = .3)


Figure2b = ggplot(air_drag, aes(airdrag,xerrorratio,color = factor(airdrag))) +
  geom_hline(linetype = 2, yintercept = 1) + 
  geom_jitter(alpha = 0.025, width = 0.1) +
  geom_flat_violin(size=1) + 
  stat_summary(fun = "median", geom = "point",size = 6, aes(group=id), shape = 95) +
  stat_summary(fun = "median", geom = "point",size = 4, aes(group=0), shape = 16, color = "black") +
  scale_color_manual(name = "", values = c("grey","gold3")) + 
  theme(legend.position = "none") + 
  #  stat_pvalue_manual(stat_test_t_b, label = "p.adj",color = "color_p") +
  labs(x = NULL,
       y = "Spatial Error ratio") 
  
Figure2b_inset = ggplot(air_drag, aes(x = airdrag,
                     y = xerrorratio_mean_by_airdrag,
                     ymin = xerrorratio_mean_by_airdrag-xerrorratio_SE_by_airdrag,
                     ymax = xerrorratio_mean_by_airdrag+xerrorratio_SE_by_airdrag,
                     color = factor(airdrag))) +
  geom_point(size = 3) +
  geom_errorbar(width=0.2, size = 1.5) +
  scale_color_manual(name = "", 
                     values = c("grey","gold3")) +
  scale_x_discrete(name ="", 
                   labels = c("Airdrag" = "AD", "NoAirdrag" = "No AD")) +
  scale_y_continuous(breaks=c(0.89,0.9,0.91)) +
  theme(legend.position = "", 
        axis.text=element_text(size=12),
        axis.title=element_text(size=12)) +
  xlab("") +
  ylab("Error Ratio")

Figure2b_Complete <-
  ggdraw() +
  draw_plot(Figure2b) +
  draw_plot(Figure2b_inset, x = 0.38, y = .6, width = .3, height = .3)


plot_grid(Figure2a_Complete,Figure2b_Complete, labels = "AUTO")
ggsave("Figure2_with_Insets.jpg", w = 12, h = 6)





####################################################################
#Are errors different for different sizes?
####################################################################
###Hypothesis 2a: Temporal Error
####################################################################

H2_Temporal <- lmer(terrorratio ~ as.factor(r) + (1|id),
                    data = air_drag[air_drag$airdrag == "Airdrag",])
summary(H2_Temporal)

H2_Temporal_Null <- lmer(terrorratio ~ (1|id), 
                         data = air_drag[air_drag$airdrag == "Airdrag",])
summary(H2_Temporal_Null)

### Compare models
anova(H2_Temporal,H2_Temporal_Null)


####################################################################
###Hypothesis 2b: Spatial Error

H2_Spatial <- lmer(xerrorratio ~ as.factor(r) + (1|id), 
                   data = air_drag[air_drag$airdrag == "Airdrag",])
summary(H2_Spatial)

H2_Spatial_Null <- lmer(xerrorratio ~ (1|id), 
                        data = air_drag[air_drag$airdrag == "Airdrag",])
summary(H2_Spatial_Null)

### Compare models
anova(H2_Spatial,H2_Spatial_Null)




######################Plots Hypothesis 2###########################
air_drag = air_drag %>% 
  group_by(r) %>%
  mutate(xerrorratio_mean_by_size = mean(xerrorratio),
         xerrorratio_SE_by_size = sd(xerrorratio)/(length(xerrorratio))^0.5,
         terrorratio_mean_by_size = mean(terrorratio),
         terrorratio_SE_by_size = sd(terrorratio)/(length(terrorratio))^0.5)

Figure5a = ggplot(air_drag %>% mutate(BallSize = case_when(r == 0.033 ~ "0.033 m", r == 0.12 ~ "0.12 m")), 
                  aes(BallSize,terrorratio,color = BallSize)) +
  geom_hline(linetype = 2, yintercept = 1) + 
  geom_jitter(alpha = 0.025, width = 0.1) +
  geom_flat_violin(size=1) + 
  stat_summary(fun = "median", geom = "point",size = 6, aes(group=id), shape = 95) +
  stat_summary(fun = "median", geom = "point",size = 4, aes(group=0), shape = 16, color = "black") +
  scale_color_manual(name = "", values = c("deeppink1","deeppink4")) +
  theme(legend.position = "none") + 
  #  stat_pvalue_manual(stat_test_t_b, label = "p.adj",color = "color_p") +
  labs(x = NULL,
       y = "Timing Error ratio") 

Figure5a_inset = ggplot(air_drag %>% mutate(BallSize = case_when(r == 0.033 ~ "0.033 m", r == 0.12 ~ "0.12 m")), 
                                      aes(x = BallSize,
                                      y = terrorratio_mean_by_size,
                                      ymin = terrorratio_mean_by_size-terrorratio_SE_by_size,
                                      ymax = terrorratio_mean_by_size+terrorratio_SE_by_size,
                                      color = factor(BallSize))) +
  geom_point(size = 3) +
  geom_errorbar(width=0.2, size = 1.5) +
  scale_color_manual(name = "", 
                     values = c("deeppink1","deeppink4")) +
  scale_x_discrete(name ="") + 
  scale_y_continuous(breaks=c(1.015,1.02,1.025)) +
  theme(legend.position = "", 
        axis.text=element_text(size=10),
        axis.title=element_text(size=10)) +
  xlab("") +
  ylab("Error Ratio")

Figure5a_Complete <-
  ggdraw() +
  draw_plot(Figure5a) +
  draw_plot(Figure5a_inset, x = 0.38, y = .6, width = .3, height = .3)


Figure5b = ggplot(air_drag %>% mutate(BallSize = case_when(r == 0.033 ~ "0.033 m", r == 0.12 ~ "0.12 m")), 
                  aes(BallSize,xerrorratio,color = BallSize)) +
  geom_hline(linetype = 2, yintercept = 1) + 
  geom_jitter(alpha = 0.025, width = 0.1) +
  geom_flat_violin(size=1) + 
  stat_summary(fun = "median", geom = "point",size = 6, aes(group=id), shape = 95) +
  stat_summary(fun = "median", geom = "point",size = 4, aes(group=0), shape = 16, color = "black") +
  scale_color_manual(name = "", values = c("deeppink1","deeppink4")) +
  theme(legend.position = "none") + 
  #  stat_pvalue_manual(stat_test_t_b, label = "p.adj",color = "color_p") +
  labs(x = NULL,
       y = "Spatial Error ratio")

  Figure5b_inset = ggplot(air_drag %>% mutate(BallSize = case_when(r == 0.033 ~ "0.033 m", r == 0.12 ~ "0.12 m")), 
                        aes(x = BallSize,
                            y = xerrorratio_mean_by_size,
                            ymin = xerrorratio_mean_by_size-xerrorratio_SE_by_size,
                            ymax = xerrorratio_mean_by_size+xerrorratio_SE_by_size,
                            color = factor(BallSize))) +
  geom_point(size = 3) +
  geom_errorbar(width=0.2, size = 1.5) +
  scale_color_manual(name = "", 
                     values = c("deeppink1","deeppink4")) +
  scale_x_discrete(name ="") + 
  scale_y_continuous(breaks=c(0.89,0.90,0.91)) +
  theme(legend.position = "", 
        axis.text=element_text(size=10),
        axis.title=element_text(size=10)) +
  xlab("") +
  ylab("Error Ratio")

Figure5b_Complete <-
  ggdraw() +
  draw_plot(Figure5b) +
  draw_plot(Figure5b_inset, x = 0.38, y = .6, width = .3, height = .3)


plot_grid(Figure5a_Complete,Figure5b_Complete, labels = "AUTO")
ggsave("Figure3 with Insets.jpg", w = 12, h = 6)

####################################################################
#Are errors dependant on contextual cues?
####################################################################
###Hypothesis 3a: Temporal Error
####################################################################

H3_Temporal_Interaction <- lmer(terrorratio ~ as.factor(r)*ball + (1|id), 
                    data = air_drag)
summary(H3_Temporal_Interaction)

H3_Temporal_Main  <- lmer(terrorratio ~ as.factor(r)+ball + (1|id), 
                         data = air_drag)
summary(H3_Temporal_Main)

### Compare models
anova(H3_Temporal_Interaction,H3_Temporal_Main)


####################################################################
###Hypothesis 3b: Spatial Error
####################################################################

H3_Spatial_Interaction <- lmer(xerrorratio ~ as.factor(r)*ball + (1|id), 
                   data = air_drag)
summary(H3_Spatial_Interaction)

H3_Spatial_Main <- lmer(xerrorratio ~ as.factor(r)+ball +(1|id), 
                        data = air_drag)
summary(H3_Spatial_Main)

### Compare models
anova(H3_Spatial_Interaction,H3_Spatial_Main)



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

air_drag = air_drag %>% 
  group_by(ball,BallXCongruency) %>%
  mutate(xerrorratio_mean_by_texture = mean(xerrorratio),
         xerrorratio_SE_by_texture = sd(xerrorratio)/(length(xerrorratio))^0.5,
         terrorratio_mean_by_texture = mean(terrorratio),
         terrorratio_SE_by_texture = sd(terrorratio)/(length(terrorratio))^0.5)

Figure3a = ggplot(air_drag, aes(BallXCongruency,terrorratio,color = as.factor(BallXCongruency))) +
  geom_hline(linetype = 2, yintercept = 1) + 
  geom_jitter(alpha = 0.025, width = 0.1) +
  geom_flat_violin(size=1) + 
  stat_summary(fun = "median", geom = "point",size = 6, aes(group=id), shape = 95) +
  stat_summary(fun = "median", geom = "point",size = 4, aes(group=0), shape = 16, color = "black") +
  scale_color_manual(name = "", values = c("red","red4","royalblue1","royalblue4")) + 
  theme(legend.position = "") + 
  labs(x = NULL,
       y = "Timing Error ratio") 


Figure3a_inset1 = ggplot(air_drag %>% filter(ball == "Basket"), aes(x = cond_size,
                                      y = terrorratio_mean_by_texture,
                                      ymin = terrorratio_mean_by_texture-terrorratio_SE_by_texture,
                                      ymax = terrorratio_mean_by_texture+terrorratio_SE_by_texture,
                                      color = factor(cond_size))) +
  geom_point(size = 3) +
  geom_errorbar(width=0.2, size = 1.5) +
  scale_color_manual(name = "", 
                     values = c("red","red4")) +
  scale_x_discrete(name ="", 
                   labels = c("cong" = "Cong", "incongr" = "Incong")) +
  scale_y_continuous(breaks=c(1.015,1.025)) +
  theme(legend.position = "", 
        axis.text=element_text(size=10),
        axis.title=element_text(size=10)) +
  xlab("") +
  ylab("Error Ratio")

Figure3a_inset2 = ggplot(air_drag %>% filter(ball == "Tennis"), aes(x = cond_size,
                                                                    y = terrorratio_mean_by_texture,
                                                                    ymin = terrorratio_mean_by_texture-terrorratio_SE_by_texture,
                                                                    ymax = terrorratio_mean_by_texture+terrorratio_SE_by_texture,
                                                                    color = factor(cond_size))) +
  geom_point(size = 3) +
  geom_errorbar(width=0.2, size = 1.5) +
  scale_color_manual(name = "", 
                     values = c("royalblue1","royalblue4")) +
  scale_x_discrete(name ="", 
                   labels = c("cong" = "Cong", "incongr" = "Incong")) +
  scale_y_continuous(breaks=c(1.015,1.025)) +
  theme(legend.position = "", 
        axis.text=element_text(size=10),
        axis.title=element_text(size=10)) +
  xlab("") +
  ylab("Error Ratio")

Figure3a_Complete <-
  ggdraw() +
  draw_plot(Figure3a) +
  draw_plot(Figure3a_inset1, x = 0.2, y = .5, width = .15, height = .4) +
  draw_plot(Figure3a_inset2, x = 0.66, y = .5, width = .15, height = .4)


Figure3b = ggplot(air_drag, aes(BallXCongruency,xerrorratio,color = as.factor(BallXCongruency))) +
  geom_hline(linetype = 2, yintercept = 1) + 
  geom_jitter(alpha = 0.025, width = 0.1) +
  geom_flat_violin(size=1) + 
  stat_summary(fun = "median", geom = "point",size = 6, aes(group=id), shape = 95) +
  stat_summary(fun = "median", geom = "point",size = 4, aes(group=0), shape = 16, color = "black") +
  scale_color_manual(name = "", values = c("red","red4","royalblue1","royalblue4")) + 
  theme(legend.position = "") + 
  labs(x = NULL,
       y = "Spatial Error ratio") 

Figure3b_inset1 = ggplot(air_drag %>% filter(ball == "Basket"), aes(x = cond_size,
                                                                    y = xerrorratio_mean_by_texture,
                                                                    ymin = xerrorratio_mean_by_texture-xerrorratio_SE_by_texture,
                                                                    ymax = xerrorratio_mean_by_texture+xerrorratio_SE_by_texture,
                                                                    color = factor(cond_size))) +
  geom_point(size = 3) +
  geom_errorbar(width=0.2, size = 1.5) +
  scale_color_manual(name = "", 
                     values = c("red","red4")) +
  scale_x_discrete(name ="", 
                   labels = c("cong" = "Cong", "incongr" = "Incong")) +
  scale_y_continuous(breaks=c(0.89,0.9,0.91)) +
  theme(legend.position = "", 
        axis.text=element_text(size=10),
        axis.title=element_text(size=10)) +
  xlab("") +
  ylab("Error Ratio")

Figure3b_inset2 = ggplot(air_drag %>% filter(ball == "Tennis"), aes(x = cond_size,
                                                                    y = xerrorratio_mean_by_texture,
                                                                    ymin = xerrorratio_mean_by_texture-xerrorratio_SE_by_texture,
                                                                    ymax = xerrorratio_mean_by_texture+xerrorratio_SE_by_texture,
                                                                    color = factor(cond_size))) +
  geom_point(size = 3) +
  geom_errorbar(width=0.2, size = 1.5) +
  scale_color_manual(name = "", 
                     values = c("royalblue1","royalblue4")) +
  scale_x_discrete(name ="", 
                   labels = c("cong" = "Cong", "incongr" = "Incong")) +
  scale_y_continuous(breaks=c(0.89,0.9,0.91)) +
  theme(legend.position = "", 
        axis.text=element_text(size=10),
        axis.title=element_text(size=10)) +
  xlab("") +
  ylab("Error Ratio")

Figure3b_Complete <-
  ggdraw() +
  draw_plot(Figure3b) +
  draw_plot(Figure3b_inset1, x = 0.2, y = .5, width = .15, height = .4) +
  draw_plot(Figure3b_inset2, x = 0.66, y = .5, width = .15, height = .4)



plot_grid(Figure3a_Complete,Figure3b_Complete, nrow = 2, labels = "AUTO")
ggsave("Figure4 with insets.jpg", w = 12, h = 8)




#######################################per subject distributions
for (i in 1:length(unique(air_drag$id))){
  if(i < 10){
    air_drag$id2[air_drag$id == unique(air_drag$id)[i]] = paste0("s0",i)    
  }
  
  else{
    air_drag$id2[air_drag$id == unique(air_drag$id)[i]] = paste0("s",i)    
  }

}
unique(air_drag$id)
unique(air_drag$id2)


Figure7a = ggplot(air_drag, aes(airdrag,terrorratio,color = factor(airdrag))) +
  geom_hline(linetype = 2, yintercept = 1) + 
  geom_jitter(alpha = 0.025, width = 0.1) +
  geom_flat_violin(size=1) + 
  scale_color_manual(name = "", values = c("grey","gold3")) + 
  scale_x_discrete(labels=c("AD","No AD")) +
  theme(legend.position = "none") + 
  #  stat_pvalue_manual(stat_test_t_b, label = "p.adj",color = "color_p") +
  labs(x = NULL,
       y = "Timing Error ratio") +
  facet_wrap(.~id2)
Figure7b = ggplot(air_drag, aes(airdrag,xerrorratio,color = factor(airdrag))) +
  geom_hline(linetype = 2, yintercept = 1) + 
  geom_jitter(alpha = 0.025, width = 0.1) +
  geom_flat_violin(size=1) + 
  scale_color_manual(name = "", values = c("grey","gold3")) + 
  scale_x_discrete(labels=c("AD","No AD")) +
  theme(legend.position = "none") + 
  #  stat_pvalue_manual(stat_test_t_b, label = "p.adj",color = "color_p") +
  labs(x = NULL,
       y = "Spatial Error ratio") +
  facet_wrap(.~id2)
plot_grid(Figure7a,Figure7b, labels = "AUTO")
ggsave("Complementary Figure 1.jpg", w = 12, h = 9)
