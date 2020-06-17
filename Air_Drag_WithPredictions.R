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

Temp_cont = Temp_cont %>%
  mutate(OcclusionX = case_when(
    occlusion == "Occluded" ~ "Occluded",
    occlusion == "Visible" ~ "Before occlusion",
    occlusion == "Occlusion period" ~ "Occlusion window")
  )

Figure1a <- ggplot(Temp_cont %>% filter(ball == "Basket"),aes(x,y+0.7, group = interaction(id_exp,label),color = OcclusionX)) +
  annotation_custom(rasterGrob(image_setup,
                               width = unit(1,"npc"),
                               height = unit(1,"npc")),
                    -Inf, Inf, -Inf, Inf) +
  geom_point(size = 1.25) +
  coord_cartesian(ylim=c(0,max(Temp_cont$y)+0.7),xlim=c(-0.25,5.5)) +
  labs(x = "x (m)",
       y = "y (m)") + 
  ggtitle("A") + 
  scale_color_manual(name = NULL, 
                     values = c(cbp1[3],alpha(cbp1[3],0.01),"red"),
                     ) + 
  #guides(color=FALSE) + 
  theme_cowplot(15) + 
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=2)))

Exp_Conds = Exp_Conds %>%
  mutate(OcclusionXsize = case_when(
    occlusion == "Occluded" & r == 0.12 ~ "Basketball, before occlusion",
    occlusion == "Visible" & r == 0.12 ~ "Basketball, visible",
    occlusion == "Occluded" & r == 0.033 ~ "Tennis ball, occluded",
    occlusion == "Visible" & r == 0.033 ~ "Tennis ball, before occlusion",
    occlusion == "Occlusion period" ~ "Occlusion window")
  )


Figure1b <- ggplot(Exp_Conds %>% filter(air_drag == 1),aes(t,vx, color = OcclusionXsize)) +
  geom_point(size = 0.75) + 
  labs(x = "t (s)",
       y = expression(V[x*","*AD]~(m/s))) +
  scale_color_manual(name = NULL, 
                     values = c(alpha(cbp1[1],0.01),
                                cbp1[1], 
                                "red",
                                cbp1[2],
                                alpha(cbp1[2],0.01)) 
  ) + 
  ggtitle("B") + 
  theme_cowplot(15) + 
  geom_hline(yintercept = c(3.0,3.5), linetype = 17) +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=2)))

Figure1c <- ggplot(Exp_Conds %>% filter(air_drag == 1),aes(t-0.001,vy-(TTC*G/2-G*(t-0.001)), color = OcclusionXsize)) +
  geom_point(size = 0.75) + 
  labs(x = "t (s)",
       y = expression(V[y*","*AD]-V[y*","*No*" "*AD]~(m/s))) +
  scale_color_manual(name = NULL, 
                     values = c(alpha(cbp1[1],0.01),
                                cbp1[1], 
                                "red",
                                cbp1[2],
                                alpha(cbp1[2],0.01)) 
                                ) + 
  ggtitle("C") + 
  geom_hline(yintercept = 0, linetype = 17) +
  theme_cowplot(15) + 
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(size=2)))

((c(0,Exp_Conds$vy[Exp_Conds$air_drag == 1]) - c(Exp_Conds$vy[Exp_Conds$air_drag == 1],0))/0.001)[length(Exp_Conds$vy[Exp_Conds$air_drag == 1])]


Temp_Exp_Conds_max = Temp_Exp_Conds_max %>% 
  mutate(airdrag = case_when(
    air_drag == "Gravity" ~ "Absent",
    air_drag == "Gravity + Air Drag" ~ "Present"),
    ballcondition = case_when(
      ball_condition == "Basket Properties" ~ "Basketball size",
      ball_condition == "Tennis Properties" ~ "Tennis ball size")
  )

Figure1d <- ggplot(Temp_Exp_Conds_max ,aes(TTC,t_max,color=factor(airdrag),shape=factor(ballcondition))) +
  geom_point(size = 3) +
  labs(x = "TTC under no Air Drag (s)",
       y = "Real TTC (s)",
       color = "Air drag:",
       shape = "Ball size:") + 
  scale_color_manual(values = c(cbp1[4],cbp1[5])) + 
  ggtitle("D") +
  theme_cowplot(15)

Figure1e <- ggplot(Temp_Exp_Conds_max,aes(TTC,x_max,color=factor(airdrag),shape=factor(ballcondition))) +
  geom_point(size = 3) +
  #facet_wrap(~ball) +
  labs(x = "TTC under no Air Drag (s)",
       y = "Horizontal point of impact (m)") + 
  scale_color_manual(values = c(cbp1[4],cbp1[5])) + 
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
confint(H1_Temporal, method = "boot")

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

(Figure2a = ggplot(air_drag,
                   aes(airdrag,terrorratio,color = factor(airdrag))) +
   geom_hline(linetype = 2, yintercept = 1) + 
   geom_jitter(alpha = 0.025, width = 0.1) +
   geom_flat_violin(size=1) + 
   stat_summary(fun = "mean", geom = "point",size = 6, aes(group=id), shape = 95) +
   stat_summary(fun = "mean", geom = "point",size = 4, aes(group=0), shape = 16, color = "black") +
   scale_color_manual(name = "", values = c("grey","gold3")) + 
   theme(legend.position = "none") + 
   #facet_wrap(TTC~id) +
   labs(x = NULL,
        y = "Timing Error ratio"))

(Figure2a_inset = ggplot(air_drag, aes(airdrag,
                                      terrorratio,
                                      color = factor(airdrag))) +
    geom_hline(yintercept = 1, linetype = 2, alpha = 0.5) +
    stat_summary(geom="point") +
    stat_summary(geom="errorbar", width=0.2) +
    #geom_errorbar(width=0.2, size = 1.5) +
    scale_color_manual(name = "", 
                     values = c("grey","gold3")) +
  scale_x_discrete(name ="", 
                   labels = c("Airdrag" = "AD", "NoAirdrag" = "No AD")) +
  scale_y_continuous(breaks=c(1.01,1.02,1.03)) +
  theme(legend.position = "", 
        axis.text=element_text(size=12),
        axis.title=element_text(size=12)) +
  xlab("") +
  ylab("Error Ratio"))

Figure2a_Complete <-  ggdraw() +
  draw_plot(Figure2a) +
  draw_plot(Figure2a_inset, x = 0.38, y = .6, width = .3, height = .3)


Figure2b = ggplot(air_drag, aes(airdrag,xerrorratio,color = factor(airdrag))) +
  geom_hline(linetype = 2, yintercept = 1) + 
  geom_jitter(alpha = 0.025, width = 0.1) +
  geom_flat_violin(size=1) + 
  stat_summary(fun = "mean", geom = "point",size = 6, aes(group=id), shape = 95) +
  stat_summary(fun = "mean", geom = "point",size = 4, aes(group=0), shape = 16, color = "black") +
  scale_color_manual(name = "", values = c("grey","gold3")) + 
  theme(legend.position = "none") + 
  #  stat_pvalue_manual(stat_test_t_b, label = "p.adj",color = "color_p") +
  labs(x = NULL,
       y = "Spatial Error ratio") 
  
Figure2b_inset = ggplot(air_drag, aes(x = airdrag,
                     y = xerrorratio,
                     color = factor(airdrag))) +
  #geom_hline(yintercept = 1, linetype = 2, alpha = 0.5) +
  stat_summary(geom="point") +
  stat_summary(geom="errorbar", width=0.2) +
  #geom_errorbar(width=0.2, size = 1.5) +
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

Figure2b_Complete <-   ggdraw() +
  draw_plot(Figure2b) +
  draw_plot(Figure2b_inset, x = 0.38, y = .6, width = .3, height = .3)


plot_grid(Figure2a_Complete,Figure2b_Complete, labels = "AUTO")
ggsave("Figure2_with_Insets.jpg", w = 12, h = 6)





####################################################################
#Are errors different for different sizes?
####################################################################
###Hypothesis 2a: Temporal Error
####################################################################

H2_Temporal <- lmer(terrorratio ~ as.factor(r)*airdrag + (1|id),
                    data = air_drag)
summary(H2_Temporal)

H2_Temporal_Null <- lmer(terrorratio ~ as.factor(r) + airdrag + (1|id), 
                         data = air_drag)
summary(H2_Temporal_Null)

### Compare models
anova(H2_Temporal,H2_Temporal_Null)


####################################################################
###Hypothesis 2b: Spatial Error

H2_Spatial <- lmer(xerrorratio ~ as.factor(r)*airdrag + (1|id), 
                   data = air_drag)
summary(H2_Spatial)

H2_Spatial_Null <- lmer(xerrorratio ~ as.factor(r) + airdrag + (1|id), 
                        data = air_drag)
summary(H2_Spatial_Null)

### Compare models
anova(H2_Spatial,H2_Spatial_Null)




######################Plots Hypothesis 2###########################
air_drag = air_drag %>%
  mutate(AirdragXSize = 
           case_when(
             airdrag == "Airdrag" & r == 0.12 ~ "Airdrag, r = 0.12m",
             airdrag == "NoAirdrag" & r == 0.033 ~ "No Airdrag, r = 0.033m",
             airdrag == "Airdrag" & r == 0.033 ~ "Airdrag, r = 0.033m",
             airdrag == "NoAirdrag" & r == 0.12 ~ "No Airdrag, r = 0.12m")
  )

air_drag = air_drag %>% 
  group_by(r,airdrag) %>%
  mutate(xerrorratio_mean_by_size = mean(xerrorratio),
         xerrorratio_SE_by_size = sd(xerrorratio)/(length(xerrorratio))^0.5,
         terrorratio_mean_by_size = mean(terrorratio),
         terrorratio_SE_by_size = sd(terrorratio)/(length(terrorratio))^0.5)

Figure3a = ggplot(air_drag, 
                  aes(AirdragXSize,terrorratio,color = AirdragXSize)) +
  geom_hline(linetype = 2, yintercept = 1) + 
  geom_jitter(alpha = 0.025, width = 0.1) +
  geom_flat_violin(size=1) + 
  stat_summary(fun = "mean", geom = "point",size = 6, aes(group=id), shape = 95) +
  stat_summary(fun = "mean", geom = "point",size = 4, aes(group=0), shape = 16, color = "black") +
  scale_color_manual(name = "", values = c("deeppink1","deeppink4","green","darkgreen")) +
  theme(legend.position = "none") + 
  #  stat_pvalue_manual(stat_test_t_b, label = "p.adj",color = "color_p") +
  labs(x = NULL,
       y = "Timing Error ratio") +
  ggtitle("A")

Figure3a_inset1 = ggplot(air_drag %>% 
                          mutate(BallSize = case_when(r == 0.033 ~ "0.033 m", r == 0.12 ~ "0.12 m")), 
                                      aes(AirdragXSize,terrorratio,color = AirdragXSize)) +
    stat_summary(geom = "point") +
    stat_summary(width=0.2, geom = "errorbar") +
    # geom_errorbar(width=0.2, size = 1.5) +
  scale_color_manual(name = "", 
                     values = c("deeppink1","deeppink4","green","darkgreen")) +
    scale_x_discrete(name ="",
                     labels = c("AD \n0.033 m","AD \n0.12 m", "No AD \n0.033 m","No AD \n0.12 m")) + 
    #scale_y_continuous(breaks=c(0.9,0.915,0.93)) +
    theme(legend.position = "none", 
          axis.text= element_text(size = rel(0.7)),
          axis.title = element_text(size = rel(0.7))) +
  xlab("") +
  ylab("Error Ratio") +
  ggtitle("B")

#Figure5a_Complete <-  ggdraw() +
#  draw_plot(Figure5a)  +
#  draw_plot(Figure5a_inset1, x = 0.6, y = .5, width = .25, height = .5) 
Figure3a_Complete = plot_grid(Figure3a,Figure3a_inset1, rel_widths = c(0.75,0.25))


Figure3b = ggplot(air_drag, 
                  aes(AirdragXSize,xerrorratio,color = AirdragXSize)) +
  geom_hline(linetype = 2, yintercept = 1) + 
  geom_jitter(alpha = 0.025, width = 0.1) +
  geom_flat_violin(size=1) + 
  stat_summary(fun = "median", geom = "point",size = 6, aes(group=id), shape = 95) +
  stat_summary(fun = "median", geom = "point",size = 4, aes(group=0), shape = 16, color = "black") +
  scale_color_manual(name = "", values = c("deeppink1","deeppink4","green","darkgreen")) +
  theme(legend.position = "none") + 
  labs(x = NULL,
       y = "Spatial Error ratio") +
  coord_cartesian(ylim=c(0.3,2.2)) +
  ggtitle("C")

Figure3b_inset1 = ggplot(air_drag %>% 
                           mutate(BallSize = case_when(r == 0.033 ~ "0.033 m", r == 0.12 ~ "0.12 m")), 
                         aes(AirdragXSize,
                             y = xerrorratio,
                             color = factor(AirdragXSize))) +
  stat_summary(geom = "point") +
  stat_summary(width=0.2, geom = "errorbar") +
  scale_color_manual(name = "", 
                     values = c("deeppink1","deeppink4","green","darkgreen")) +
  scale_x_discrete(name ="",
                   labels = c("AD \n0.033 m","AD \n0.12 m", "No AD \n0.033 m","No AD \n0.12 m")) + 
    theme(legend.position = "none", 
          axis.text= element_text(size = rel(0.7)),
          axis.title = element_text(size = rel(0.7))) +
  xlab("") +
  ylab("Error Ratio") +
  ggtitle("D")
Figure3b_Complete = plot_grid(Figure3b,Figure3b_inset1, rel_widths = c(0.75,0.25))

#Figure5b_Complete <-  ggdraw() +
#  draw_plot(Figure5b) +
#  draw_plot(Figure5b_inset1, x = 0.6, y = .5, width = .25, height = .5) 

plot_grid(Figure3a_Complete,Figure3b_Complete, nrow = 2, labels = "AUTO")
ggsave("Figure3.jpg", w = 12, h = 8)




####################################################################
#Are errors dependant on contextual cues?
####################################################################
###Hypothesis 3a: Temporal Error
####################################################################
H3_Temporal_Interaction <- lmer(terrorratio ~ as.factor(r)*ball*airdrag + (1|id), 
                    data = air_drag)
summary(H3_Temporal_Interaction)

H3_Temporal_Main  <- lmer(terrorratio ~ as.factor(r)*airdrag + ball*airdrag + (1|id), 
                         data = air_drag)
summary(H3_Temporal_Main)

### Compare models
anova(H3_Temporal_Interaction,H3_Temporal_Main)


####################################################################
###Hypothesis 3b: Spatial Error
####################################################################

H3_Spatial_Interaction <- lmer(xerrorratio ~ as.factor(r)*ball*airdrag + (1|id), 
                   data = air_drag)
summary(H3_Spatial_Interaction)

H3_Spatial_Main <- lmer(xerrorratio ~ as.factor(r)*airdrag + ball*airdrag +(1|id), 
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

Figure4a = ggplot(air_drag, aes(BallXCongruency,terrorratio,color = as.factor(BallXCongruency))) +
  geom_hline(linetype = 2, yintercept = 1) + 
  geom_jitter(alpha = 0.025, width = 0.1) +
  geom_flat_violin(size=1) + 
  stat_summary(fun = "median", geom = "point",size = 6, aes(group=id), shape = 95) +
  stat_summary(fun = "median", geom = "point",size = 4, aes(group=0), shape = 16, color = "black") +
  scale_color_manual(name = "", values = c("red","red4","royalblue1","royalblue4")) + 
  theme(legend.position = "") + 
  labs(x = NULL,
       y = "Timing Error Ratio")  +
  ggtitle("A") +
  facet_wrap(air_drag~.)


Figure4a_inset1 = ggplot(air_drag , aes(x = BallXCongruency,
                                      y = terrorratio,
                                      color = factor(BallXCongruency))) +
  stat_summary(geom="point") +
  stat_summary(width = 0.2, geom="errorbar") +
  scale_color_manual(name = "", 
                     values = c("red","red4","royalblue1","royalblue4")) +
  scale_x_discrete(name ="", 
                   labels = paste0(rep(c("Basket","Tennis"),each=2), c("\nCong.","\nIncong."))) +
  #scale_y_continuous(breaks=c(0.89,0.9,0.91)) +
  theme(legend.position = "none",
        axis.text= element_text(size = rel(0.7)),
        axis.title = element_text(size = rel(0.7))) +
  xlab("") +
  ylab("Timing Error Ratio") +
  ggtitle("B") +
  facet_wrap(air_drag~.)

Figure4aa = plot_grid(Figure4a,Figure4a_inset1,rel_widths = c(0.75,0.25))

#Figure4a_Complete <- ggdraw() +
#  draw_plot(Figure3a) +
#  draw_plot(Figure3a_inset1, x = 0.6, y = .5, width = .25, height = .5) 

Figure4b = ggplot(air_drag, aes(BallXCongruency,xerrorratio,color = as.factor(BallXCongruency))) +
  geom_hline(linetype = 2, yintercept = 1) + 
  geom_jitter(alpha = 0.025, width = 0.1) +
  geom_flat_violin(size=1) + 
  stat_summary(fun = "mean", geom = "point",size = 6, aes(group=id), shape = 95) +
  stat_summary(fun = "mean", geom = "point",size = 4, aes(group=0), shape = 16, color = "black") +
  scale_color_manual(name = "", values = c("red","red4","royalblue1","royalblue4")) + 
  theme(legend.position = "") + 
  labs(x = NULL,
       y = "Spatial Error Ratio")  +
  ggtitle("C") +
  facet_wrap(air_drag~.)

Figure4b_inset1 = ggplot(air_drag,
                         aes(x = BallXCongruency,
                             y = xerrorratio,
                             color = factor(BallXCongruency))) +
  #geom_hline(yintercept = 1, linetype = 2, alpha = 0.5) +
  stat_summary(geom="point") +
  stat_summary(width = 0.2, geom="errorbar") +
  #geom_errorbar(width=0.2, size = 1.5) +
  scale_color_manual(name = "", 
                     values = c("red","red4","royalblue1","royalblue4")) +
  scale_x_discrete(name ="", 
                   labels = paste0(rep(c("Basket","Tennis"),each=2), c("\nCong.","\nIncong."))) +
  #scale_y_continuous(breaks=c(0.89,0.9,0.91)) +
  theme(legend.position = "none", 
        axis.text= element_text(size = rel(0.7)),
        axis.title = element_text(size = rel(0.7))) +
  xlab("") +
  ylab("Spatial Error Ratio") +
  ggtitle("D") +
  facet_wrap(air_drag~.)

#Figure3b_Complete <- ggdraw() +
#  draw_plot(Figure3b) +
#  draw_plot(Figure3b_inset1, x = 0.6, y = .5, width = .25, height = .5) 

Figure4bb = plot_grid(Figure4b,Figure4b_inset1,rel_widths = c(0.75,0.25))

plot_grid(Figure4aa,Figure4bb, nrow = 2)
ggsave("Figure4.jpg", w = 12, h = 8)


####################################################################
#############   Mean effect of AD         ##########################
####################################################################

# Exp_Conds %>%
#   filter(air_drag == 1) %>%
#   #select(-vx) %>% 
#   ungroup() %>% 
#   summarise(mean_ay = mean(ay,na.rm=T)+9.807,
#             mean_ax = mean(ax,na.rm=T))
# 
# `mean(ay, na.rm = T) + 9.807` `mean(ax, na.rm = T)`
# <dbl>                 <dbl>
#   1                      -0.00989                -0.296
####################################################################




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





###################################
###################################
######Predictions for Hypotheses###
#(Seguro que tu lo harías todo en dos líneas de código)
#(en principio solo lo hice para confirmar una mini-cosa y no quería pasar mucho tiempo haciendolo, 
#pero luego ... bueno ... es que con unas predicciones tan complicadas pues mejor hacerlo con modelos)
#Y mira, aqui está el código feito



###expected differences Hypothesis 1
#get predictions for extrapolated distance:
###first I add a couple of values that arent in air_drag, then I take one line for each combination of relevance
GetModelResponses_AirdragModel = air_drag %>%
  mutate(cd = case_when(
    airdrag == "Airdrag" ~ 0.535,
    airdrag == "NoAirdrag" ~ 0),
    vy = case_when( #initial vertical velocity
      TTC == 1.0 ~ 4.9035,
      TTC == 1.2 ~ 5.8842,
      TTC == 1.4 ~ 6.8649),
    m = case_when( #mass of object
      r == 0.033 ~ 0.06,
      r == 0.12 ~ 0.6),
    A = pi * r^2, #area
    rho = 1.225, #rho, whatever that is
    D = rho*cd*A/2, #resistance factor, is zero when there is no AD simulated, and some positive value when AD is simulated
    D_WithAD = rho*0.535*A/2) %>% #resistance factor for second part of trajectory, air drag assumption
  group_by(airdrag,TTC,r,vx,cond_size) %>%
  slice(1)

i = 1
dt = 0.001
GetModelResponses_AirdragModel$Last_x = 0
GetModelResponses_AirdragModel$Last_y = 0
GetModelResponses_AirdragModel$Last_t = 0
GetModelResponses_AirdragModel$xAtDisappearance = 0
GetModelResponses_AirdragModel$Last_vx = 0
GetModelResponses_AirdragModel$Last_vy = 0
GetModelResponses_AirdragModel$vxAtDisappearance = 0
GetModelResponses_AirdragModel$vyAtDisappearance = 0

####the loop simulates the trajectories for each of the 48 stimulus combinations
#and saves the several values in the same dataframe
#this version assumes an AD model, that is, after disappearance (after 0.575% of the trajectory)
#the rest of the trajectory is computed according to gravity 
for (i in (1:length(GetModelResponses_AirdragModel$trial))){

  D = GetModelResponses_AirdragModel[i,]$D
  D2 = GetModelResponses_AirdragModel[i,]$D_WithAD
  vx = GetModelResponses_AirdragModel[i,]$vx
  vy = GetModelResponses_AirdragModel[i,]$vy
  m = GetModelResponses_AirdragModel[i,]$m
  TTC = GetModelResponses_AirdragModel[i,]$TTC
  OccludedPercentage = 0.575
  
  g = 9.807
  t = 0
  y = 0
  x = 0
  Got_X = 0
  
  for (j in seq(0,1.5,dt)){
    
    if (y >= 0){
      
      if (t > OccludedPercentage*TTC){
        
        D = D2
        
        if (Got_X == 0){
          xAtDisappearance = x #get lastx before disappearance
          vxAtDisappearance = vx #get last horizontal velocity before disappearance (only for checking)
          vyAtDisappearance = vy #get last vertical velocity before disappearance (only for checking)
          Got_X = 1 #
        }
      }
      
      v = sqrt(vx^2+vy^2)
      
      ax <- -(D/m)*v*vx
      ay <- -g - (D/m)*v*vy
      
      vx = vx + ax * dt 
      vy = vy + ay * dt 
      
      x = x + vx * dt + 0.5*ax*dt^2
      y = y + vy * dt + 0.5*ay*dt^2

      t = t + dt
    }
  }
  GetModelResponses_AirdragModel$Last_x[i] = x #x at response given
  GetModelResponses_AirdragModel$Last_y[i] = y #y at response given
  GetModelResponses_AirdragModel$Last_t[i] = t #t at response given
  GetModelResponses_AirdragModel$Last_vx[i] = vx
  GetModelResponses_AirdragModel$Last_vy[i] = vy
  GetModelResponses_AirdragModel$xAtDisappearance[i] = xAtDisappearance
  GetModelResponses_AirdragModel$vxAtDisappearance[i] = vxAtDisappearance
  GetModelResponses_AirdragModel$vyAtDisappearance[i] = vyAtDisappearance
}

#here we get another couple of values and save everything into a new dataframe with the results under the assumption of a (correct) AD model
GetModelResponses_AirdragModel = GetModelResponses_AirdragModel %>%
  select(airdrag,t_max,vy,r,vx,Last_x,Last_t,Last_vx,Last_vy,xAtDisappearance,vxAtDisappearance,vyAtDisappearance,condsize,x_max) %>%
  mutate(t_max = t_max-0.001, #t max is taken from original dataframe and should be fine. 
         #Subtracted 0.001 because the for loops used to simulated the trajectories above loop once more than hthey should (I think)
         tAtDisappearance = TTC*0.575,
         ErrorSpace = Last_x-x_max,
         ErrorTime = Last_t-t_max,
         ExtrapolatedSpace_Model = Last_x - xAtDisappearance,
         ExtrapolatedTime_Model = Last_t - tAtDisappearance,
         ErrorRatioSpace = (ErrorSpace+ExtrapolatedSpace_Model)/ExtrapolatedSpace_Model,
         ErrorRatioTime = (ErrorTime+ExtrapolatedTime_Model)/ExtrapolatedTime_Model)


###heres the really ugly part: I make new dataframes for each of the response conditions
##get spatial error for "Air Drag: Present" and "Air Drag: Absent"
ADModel_Space_Ballsize = data.frame(label = c("(2) No AD",
                                             "(1) AD"),
                                   value = c(mean((GetModelResponses_AirdragModel %>% 
                                                     filter(airdrag == "NoAirdrag"))$ErrorRatioSpace),
                                             mean((GetModelResponses_AirdragModel %>% 
                                                     filter(airdrag == "Airdrag"))$ErrorRatioSpace)))

ADModel_Space_Ballsize$Model = "Air Drag Model"
ADModel_Space_Ballsize$Modality = "Space"

##get timing error for "Air Drag: Present" and "Air Drag: Absent"
ADModel_Time_Ballsize = data.frame(label = c("(2) No AD",
                                              "(1) AD"),
                                    value = c(mean((GetModelResponses_AirdragModel %>% 
                                                      filter(airdrag == "NoAirdrag"))$ErrorRatioTime),
                                              mean((GetModelResponses_AirdragModel %>% 
                                                      filter(airdrag == "Airdrag"))$ErrorRatioSpace)))

ADModel_Time_Ballsize$Model = "Air Drag Model"
ADModel_Time_Ballsize$Modality = "Time"


###expected differences Hypothesis 1
#get predictions for extrapolated distance with gravity only (o sea, inverse of what Borja did in the original code)
#as above, but drag coefficient changes to 0 at disappearance for all trajectories, for a no air drag model
GetModelResponses_NoAirdragModel = air_drag %>%
  mutate(cd = case_when(
    airdrag == "Airdrag" ~ 0.535,
    airdrag == "NoAirdrag" ~ 0),
    vy = case_when(
      TTC == 1.0 ~ 4.9035,
      TTC == 1.2 ~ 5.8842,
      TTC == 1.4 ~ 6.8649),
    m = case_when(
      r == 0.033 ~ 0.06,
      r == 0.12 ~ 0.6),
    A = pi * r^2,
    rho = 1.225,
    D = rho*cd*A/2,
    D_WithAD = rho*0.535*A/2) %>% 
  group_by(airdrag,TTC,r,vx,condsize) %>% 
  slice(1)

i = 1
dt = 0.001
GetModelResponses_NoAirdragModel$Last_x = 0
GetModelResponses_NoAirdragModel$Last_y = 0
GetModelResponses_NoAirdragModel$Last_t = 0
GetModelResponses_NoAirdragModel$xAtDisappearance = 0
GetModelResponses_NoAirdragModel$Last_vx = 0
GetModelResponses_NoAirdragModel$Last_vy = 0
GetModelResponses_NoAirdragModel$vxAtDisappearance = 0
GetModelResponses_NoAirdragModel$vyAtDisappearance = 0

for (i in (1:length(GetModelResponses_NoAirdragModel$trial))){
  
  D = GetModelResponses_NoAirdragModel[i,]$D
  D2 = GetModelResponses_NoAirdragModel[i,]$D_WithAD
  vx = GetModelResponses_NoAirdragModel[i,]$vx
  vy = GetModelResponses_NoAirdragModel[i,]$vy
  m = GetModelResponses_NoAirdragModel[i,]$m
  TTC = GetModelResponses_NoAirdragModel[i,]$TTC
  OccludedPercentage = 0.575
  
  g = 9.807
  t = 0
  y = 0
  x = 0
  Got_X = 0
  
  for (j in seq(0,1.5,dt)){
    
    if (y >= 0){
      
      if (t > OccludedPercentage*TTC){
        
        D = 0
        
        if (Got_X == 0){
          xAtDisappearance = x
          vxAtDisappearance = vx
          vyAtDisappearance = vy
          Got_X = 1
        }
      }
      
      
      v = sqrt(vx^2+vy^2)
      
      ax <- -(D/m)*v*vx
      ay <- -g - (D/m)*v*vy
      
      vx = vx + ax * dt 
      vy = vy + ay * dt 
      
      x = x + vx * dt + 0.5*ax*dt^2
      y = y + vy * dt + 0.5*ay*dt^2
      
      t = t + dt
    }
  }
  GetModelResponses_NoAirdragModel$Last_x[i] = x
  GetModelResponses_NoAirdragModel$Last_y[i] = y
  GetModelResponses_NoAirdragModel$Last_t[i] = t
  GetModelResponses_NoAirdragModel$Last_vx[i] = vx
  GetModelResponses_NoAirdragModel$Last_vy[i] = vy
  GetModelResponses_NoAirdragModel$xAtDisappearance[i] = xAtDisappearance
  GetModelResponses_NoAirdragModel$vxAtDisappearance[i] = vxAtDisappearance
  GetModelResponses_NoAirdragModel$vyAtDisappearance[i] = vyAtDisappearance
}

GetModelResponses_NoAirdragModel = GetModelResponses_NoAirdragModel %>%
  select(airdrag,TTC,vy,r,vx,Last_x,Last_t,Last_vx,Last_vy,xAtDisappearance,vxAtDisappearance,vyAtDisappearance,x_max,t_max,condsize,t_max,x_max) %>%
  mutate(t_max = t_max-0.001,
         tAtDisappearance = TTC*0.575,
         ErrorSpace = Last_x-x_max,
         ErrorTime = Last_t-t_max,
         ExtrapolatedSpace_Model = Last_x - xAtDisappearance,
         ExtrapolatedTime_Model = Last_t - tAtDisappearance,
         ErrorRatioSpace = (ErrorSpace+ExtrapolatedSpace_Model)/ExtrapolatedSpace_Model,
         ErrorRatioTime = (ErrorTime+ExtrapolatedTime_Model/ExtrapolatedTime_Model)
  )


GravityModel_Space_Ballsize = data.frame(label = c("(2) No AD",
                                                    "(1) AD"),
                                          value = c(mean((GetModelResponses_NoAirdragModel %>% 
                                                            filter(airdrag == "NoAirdrag"))$ErrorRatioSpace),
                                                    mean((GetModelResponses_NoAirdragModel %>% 
                                                            filter(airdrag == "Airdrag"))$ErrorRatioSpace)))

GravityModel_Space_Ballsize$Model = "No Air Drag Model"
GravityModel_Space_Ballsize$Modality = "Space"

GravityModel_Timing_Ballsize = data.frame(label = c("(2) No AD",
                                                    "(1) AD"),
                                          value = c(mean((GetModelResponses_NoAirdragModel %>% 
                                                            filter(airdrag == "NoAirdrag"))$ErrorRatioTime),
                                                    mean((GetModelResponses_NoAirdragModel %>% 
                                                           filter(airdrag == "Airdrag"))$ErrorRatioTime)))
GravityModel_Timing_Ballsize$Model = "No Air Drag Model"
GravityModel_Timing_Ballsize$Modality = "Time"

###put all of the above dataframes together for a full dataframe with all relevant (?) values for Hypothesis 1
DF_Predictions_H1 = rbind(ADModel_Space_Ballsize,
                          ADModel_Time_Ballsize,
                          GravityModel_Space_Ballsize,
                          GravityModel_Timing_Ballsize)




###expected differences Hypothesis 2
#get predictions for extrapolated distance:
#as above, but drag coefficient changes to a mean value for both objects at disappearance
#so extrapolation would happen the same for both big and small targets, with the mean size and mean mass between both;
#if only size changes, we get much bigger differences,
#and it might actually be a better alternative hypothesis than the other one
#but Im not sure

GetModelResponses_MeanR = air_drag %>%
  mutate(cd = case_when(
    airdrag == "Airdrag" ~ 0.535,
    airdrag == "NoAirdrag" ~ 0),
    vy = case_when(
      TTC == 1.0 ~ 4.9035,
      TTC == 1.2 ~ 5.8842,
      TTC == 1.4 ~ 6.8649),
    m = case_when(
      r == 0.033 ~ 0.06,
      r == 0.12 ~ 0.6),
    m2 = (0.06+0.6)/2, #mean mass
    A = pi * r^2,
    A2 = pi * ((0.033+0.12)/2)^2, #mean area
    rho = 1.225,
    D = rho*cd*A/2,
    D_WithAD_MeanR = rho*0.535*A2/2) %>% 
  group_by(airdrag,TTC,r,vx) %>% 
  slice(1)

i = 1
dt = 0.001
GetModelResponses_MeanR$Last_x = 0
GetModelResponses_MeanR$Last_y = 0
GetModelResponses_MeanR$Last_t = 0
GetModelResponses_MeanR$xAtDisappearance = 0
GetModelResponses_MeanR$Last_vx = 0
GetModelResponses_MeanR$Last_vy = 0
GetModelResponses_MeanR$vxAtDisappearance = 0
GetModelResponses_MeanR$vyAtDisappearance = 0

for (i in (1:length(GetModelResponses_MeanR$trial))){
  
  D = GetModelResponses_MeanR[i,]$D
  D2 = GetModelResponses_MeanR[i,]$D_WithAD_MeanR
  vx = GetModelResponses_MeanR[i,]$vx
  vy = GetModelResponses_MeanR[i,]$vy
  m = GetModelResponses_MeanR[i,]$m
  m2 = GetModelResponses_MeanR[i,]$m2
  TTC = GetModelResponses_MeanR[i,]$TTC
  OccludedPercentage = 0.575
  
  g = 9.807
  t = 0
  y = 0
  x = 0
  Got_X = 0
  
  for (j in seq(0,1.5,dt)){
    
    if (y >= 0){
      
      if (t > OccludedPercentage*TTC){
        
        D = D2
        m = m2
        
        if (Got_X == 0){
          xAtDisappearance = x
          vxAtDisappearance = vx
          vyAtDisappearance = vy
          Got_X = 1
        }
      }
      
      v = sqrt(vx^2+vy^2)
      
      ax <- -(D/m)*v*vx
      ay <- -g - (D/m)*v*vy
      
      vx = vx + ax * dt 
      vy = vy + ay * dt 
      
      x = x + vx * dt + 0.5*ax*dt^2
      y = y + vy * dt + 0.5*ay*dt^2
      
      t = t + dt
    }
  }
  GetModelResponses_MeanR$Last_x[i] = x
  GetModelResponses_MeanR$Last_y[i] = y
  GetModelResponses_MeanR$Last_t[i] = t
  GetModelResponses_MeanR$Last_vx[i] = vx
  GetModelResponses_MeanR$Last_vy[i] = vy
  GetModelResponses_MeanR$xAtDisappearance[i] = xAtDisappearance
  GetModelResponses_MeanR$vxAtDisappearance[i] = vxAtDisappearance
  GetModelResponses_MeanR$vyAtDisappearance[i] = vyAtDisappearance
}

GetModelResponses_MeanR = GetModelResponses_MeanR %>%
  select(airdrag,TTC,vy,r,vx,Last_x,Last_t,Last_vx,Last_vy,xAtDisappearance,vxAtDisappearance,vyAtDisappearance,condsize,t_max,x_max) %>%
  mutate(t_max = t_max-0.001,
         tAtDisappearance = TTC*0.575,
         ErrorSpace = Last_x-x_max,
         ErrorTime = Last_t-t_max,
         ExtrapolatedSpace_Model = Last_x - xAtDisappearance,
         ExtrapolatedTime_Model = Last_t - tAtDisappearance,
         ErrorRatioSpace = (ErrorSpace+ExtrapolatedSpace_Model)/ExtrapolatedSpace_Model,
         ErrorRatioTime = (ErrorTime+ExtrapolatedTime_Model)/ExtrapolatedTime_Model)

#spatial errors if we have extrapolate motion with these mean values
ADModel_Space_Ballsize_MeanR = data.frame(label = c("(4) 0.12 No AD",
                                              "(3) 0.033 No AD",
                                              "(2) 0.12 AD",
                                              "(1) 0.033 AD"),
                                    value = c(mean((GetModelResponses_MeanR %>% 
                                                      filter(airdrag == "NoAirdrag" & r == 0.12))$ErrorRatioSpace),
                                              mean((GetModelResponses_MeanR %>% 
                                                      filter(airdrag == "NoAirdrag" & r == 0.033))$ErrorRatioSpace),
                                              mean((GetModelResponses_MeanR %>% 
                                                      filter(airdrag == "Airdrag" & r == 0.12))$ErrorRatioSpace),
                                              mean((GetModelResponses_MeanR %>% 
                                                      filter(airdrag == "Airdrag" & r == 0.033))$ErrorRatioSpace))
)

ADModel_Space_Ballsize_MeanR$Model = "Air Drag Model"
ADModel_Space_Ballsize_MeanR$Modality = "Space"
ADModel_Space_Ballsize_MeanR$Condition = "Mean Size"

#time errors if we have extrapolate motion with these mean values AND we use an air drag model
ADModel_Time_Ballsize_MeanR = data.frame(label = c("(4) 0.12 No AD",
                                             "(3) 0.033 No AD",
                                             "(2) 0.12 AD",
                                             "(1) 0.033 AD"),
                                   value = c(mean((GetModelResponses_MeanR %>% 
                                                     filter(airdrag == "NoAirdrag" & r == 0.12))$ErrorRatioTime),
                                             mean((GetModelResponses_MeanR %>% 
                                                     filter(airdrag == "NoAirdrag" & r == 0.033))$ErrorRatioTime),
                                             mean((GetModelResponses_MeanR %>% 
                                                     filter(airdrag == "Airdrag" & r == 0.12))$ErrorRatioTime),
                                             mean((GetModelResponses_MeanR %>% 
                                                     filter(airdrag == "Airdrag" & r == 0.033))$ErrorRatioTime))
                                   )
ADModel_Time_Ballsize_MeanR$Model = "Air Drag Model"
ADModel_Time_Ballsize_MeanR$Modality = "Time"
ADModel_Time_Ballsize_MeanR$Condition = "Mean Size"

#spatial errors if we use a no air drag model and use the actual values
GravityModel_Space_Ballsize = data.frame(label = c("(4) 0.12 No AD",
                                                        "(3) 0.033 No AD",
                                                        "(2) 0.12 AD",
                                                        "(1) 0.033 AD"),
                                              value = c(mean((GetModelResponses_NoAirdragModel %>% 
                                                                filter(airdrag == "NoAirdrag" & r == 0.12))$ErrorRatioSpace),
                                                        mean((GetModelResponses_NoAirdragModel %>% 
                                                                filter(airdrag == "NoAirdrag" & r == 0.033))$ErrorRatioSpace),
                                                        mean((GetModelResponses_NoAirdragModel %>% 
                                                                filter(airdrag == "Airdrag" & r == 0.12))$ErrorRatioSpace),
                                                        mean((GetModelResponses_NoAirdragModel %>% 
                                                                filter(airdrag == "Airdrag" & r == 0.033))$ErrorRatioSpace))
)
GravityModel_Space_Ballsize$Model = "No Air Drag Model"
GravityModel_Space_Ballsize$Modality = "Space"
GravityModel_Space_Ballsize$Condition = "Correct Size"

#time errors if we use a no air drag model and use the actual values
GravityModel_Timing_Ballsize = data.frame(label = c("(4) 0.12 No AD",
                                                         "(3) 0.033 No AD",
                                                         "(2) 0.12 AD",
                                                         "(1) 0.033 AD"),
                                               value = c(mean((GetModelResponses_NoAirdragModel %>% 
                                                                 filter(airdrag == "NoAirdrag" & r == 0.12))$ErrorRatioTime),
                                                         mean((GetModelResponses_NoAirdragModel %>% 
                                                                 filter(airdrag == "NoAirdrag" & r == 0.033))$ErrorRatioTime),
                                                         mean((GetModelResponses_NoAirdragModel %>% 
                                                                 filter(airdrag == "Airdrag" & r == 0.12))$ErrorRatioTime),
                                                         mean((GetModelResponses_NoAirdragModel %>% 
                                                                 filter(airdrag == "Airdrag" & r == 0.033))$ErrorRatioTime))
                                               )
GravityModel_Timing_Ballsize$Model = "No Air Drag Model"
GravityModel_Timing_Ballsize$Modality = "Time"
GravityModel_Timing_Ballsize$Condition = "Correct Size"

#spatial errors if we use an air drag model and use the actual values
ADModel_Space_Ballsize_CorrectR = data.frame(label = c("(4) 0.12 No AD",
                                                             "(3) 0.033 No AD",
                                                             "(2) 0.12 AD",
                                                             "(1) 0.033 AD"),
                                                   value = c(mean((GetModelResponses_AirdragModel %>% 
                                                                     filter(airdrag == "NoAirdrag" & r == 0.12))$ErrorRatioSpace),
                                                             mean((GetModelResponses_AirdragModel %>% 
                                                                     filter(airdrag == "NoAirdrag" & r == 0.033))$ErrorRatioSpace),
                                                             mean((GetModelResponses_AirdragModel %>% 
                                                                     filter(airdrag == "Airdrag" & r == 0.12))$ErrorRatioSpace),
                                                             mean((GetModelResponses_AirdragModel %>% 
                                                                     filter(airdrag == "Airdrag" & r == 0.033))$ErrorRatioSpace))
)
ADModel_Space_Ballsize_CorrectR$Model = "Air Drag Model"
ADModel_Space_Ballsize_CorrectR$Modality = "Space"
ADModel_Space_Ballsize_CorrectR$Condition = "Correct Size"

#time errors if we use an air drag model and use the actual values
ADModel_Timing_Ballsize_CorrectR = data.frame(label = c("(4) 0.12 No AD",
                                                              "(3) 0.033 No AD",
                                                              "(2) 0.12 AD",
                                                              "(1) 0.033 AD"),
                                                    value = c(mean((GetModelResponses_AirdragModel %>% 
                                                                      filter(airdrag == "NoAirdrag" & r == 0.12))$ErrorRatioTime),
                                                              mean((GetModelResponses_AirdragModel %>% 
                                                                      filter(airdrag == "NoAirdrag" & r == 0.033))$ErrorRatioTime),
                                                              mean((GetModelResponses_AirdragModel %>% 
                                                                      filter(airdrag == "Airdrag" & r == 0.12))$ErrorRatioTime),
                                                              mean((GetModelResponses_AirdragModel %>% 
                                                                      filter(airdrag == "Airdrag" & r == 0.033))$ErrorRatioTime))
                                                  )
ADModel_Timing_Ballsize_CorrectR$Model = "Air Drag Model"
ADModel_Timing_Ballsize_CorrectR$Modality = "Time"
ADModel_Timing_Ballsize_CorrectR$Condition = "Correct Size"

###put all of the above dataframes together for a full dataframe with all relevant (?) values for Hypothesis 2
DF_Predictions_H2 = rbind(ADModel_Space_Ballsize_MeanR,
                          ADModel_Time_Ballsize_MeanR,
                          GravityModel_Space_Ballsize,
                          GravityModel_Timing_Ballsize,
                          ADModel_Space_Ballsize_CorrectR,
                          ADModel_Timing_Ballsize_CorrectR)

#DF_Predictions_H2_MeanR = DF_Predictions_H2_MeanR %>% 
#  mutate(R = case_when(label == "(4) 0.12 No AD" ~ 0.12,
#                       label == "(3) 0.033 No AD" ~ 0.033,
#                       label == "(2) 0.12 AD" ~ 0.12,
#                       label == "(1) 0.033 AD" ~ 0.033),
#         Airdrag = case_when(label == "(4) 0.12 No AD" ~ "Air Drag: Absent",
#                             label == "(3) 0.033 No AD" ~ "Air Drag: Absent",
#                             label == "(2) 0.12 AD" ~ "Air Drag: Present",
#                             label == "(1) 0.033 AD" ~ "Air Drag: Present"))



###expected differences Hypothesis 3
#get predictions for extrapolated distance:
#here, we assume that an incongruent texture biases extrapolation towards the texture 
#(so basketball texture biases towards basketball properties by a bit)
GetModelResponses_Texture = air_drag %>%
  mutate(cd = case_when(
    airdrag == "Airdrag" ~ 0.535,
    airdrag == "NoAirdrag" ~ 0),
    vy = case_when(
      TTC == 1.0 ~ 4.9035,
      TTC == 1.2 ~ 5.8842,
      TTC == 1.4 ~ 6.8649),
    m = case_when(
      r == 0.033 ~ 0.06,
      r == 0.12 ~ 0.6),
    r2 = case_when( #biased size percept
      r == 0.033 & condsize == "Incongruent" ~ 0.05,
      r == 0.12 & condsize == "Incongruent" ~ 0.09,
      r == 0.033 & condsize == "Congruent"   ~ 0.033,
      r == 0.12 & condsize == "Congruent" ~ 0.12),
    m2 = case_when( #biased size percept
      r == 0.033 & condsize == "Incongruent" ~ 0.2,
      r == 0.12 & condsize == "Incongruent" ~ 0.4,
      r == 0.033 & condsize == "Congruent"   ~ 0.06,
      r == 0.12 & condsize == "Congruent" ~ 0.6),
    A = pi * r^2,
    A2 = pi * r2^2,
    rho = 1.225,
    D = rho*cd*A/2, #correct coefficient
    D2 = rho*0.535*A2/2) %>% #texture-biased coefficient 
  group_by(airdrag,TTC,r,vx,condsize) %>% 
  slice(1)

i = 1
dt = 0.001
GetModelResponses_Texture$Last_x = 0
GetModelResponses_Texture$Last_y = 0
GetModelResponses_Texture$Last_t = 0
GetModelResponses_Texture$xAtDisappearance = 0
GetModelResponses_Texture$Last_vx = 0
GetModelResponses_Texture$Last_vy = 0
GetModelResponses_Texture$vxAtDisappearance = 0
GetModelResponses_Texture$vyAtDisappearance = 0

for (i in (1:length(GetModelResponses_Texture$trial))){
  
  D = GetModelResponses_Texture[i,]$D
  D2 = GetModelResponses_Texture[i,]$D2
  vx = GetModelResponses_Texture[i,]$vx
  vy = GetModelResponses_Texture[i,]$vy
  m = GetModelResponses_Texture[i,]$m
  m2 = GetModelResponses_Texture[i,]$m2
  TTC = GetModelResponses_Texture[i,]$TTC
  OccludedPercentage = 0.575
  
  g = 9.807
  t = 0
  y = 0
  x = 0
  Got_X = 0
  
  for (j in seq(0,1.5,dt)){
    
    if (y >= 0){
      
      if (t > OccludedPercentage*TTC){
        
        D = D2
        m = m2
        
        if (Got_X == 0){
          xAtDisappearance = x
          vxAtDisappearance = vx
          vyAtDisappearance = vy
          Got_X = 1
        }
      }
      
      v = sqrt(vx^2+vy^2)
      
      ax <- -(D/m)*v*vx
      ay <- -g - (D/m)*v*vy
      
      vx = vx + ax * dt 
      vy = vy + ay * dt 
      
      x = x + vx * dt + 0.5*ax*dt^2
      y = y + vy * dt + 0.5*ay*dt^2
      
      t = t + dt
    }
  }
  GetModelResponses_Texture$Last_x[i] = x
  GetModelResponses_Texture$Last_y[i] = y
  GetModelResponses_Texture$Last_t[i] = t
  GetModelResponses_Texture$Last_vx[i] = vx
  GetModelResponses_Texture$Last_vy[i] = vy
  GetModelResponses_Texture$xAtDisappearance[i] = xAtDisappearance
  GetModelResponses_Texture$vxAtDisappearance[i] = vxAtDisappearance
  GetModelResponses_Texture$vyAtDisappearance[i] = vyAtDisappearance
}

GetModelResponses_Texture = GetModelResponses_Texture %>%
  select(airdrag,TTC,vy,r,vx,Last_x,Last_t,Last_vx,Last_vy,xAtDisappearance,vxAtDisappearance,vyAtDisappearance,condsize,t_max,x_max) %>%
  mutate(t_max = t_max-0.001,
         tAtDisappearance = TTC*0.575,
         ErrorSpace = Last_x-x_max,
         ErrorTime = Last_t-t_max,
         ExtrapolatedSpace_Model = Last_x - xAtDisappearance,
         ExtrapolatedTime_Model = Last_t - tAtDisappearance,
         ErrorRatioSpace = (ErrorSpace+ExtrapolatedSpace_Model)/ExtrapolatedSpace_Model,
         ErrorRatioTime = (ErrorTime+ExtrapolatedTime_Model)/ExtrapolatedTime_Model)

#air drag model AND effect of texture
ADModel_Space_Texture_Effect = data.frame(label = c("(4) 0.12 No AD",
                                              "(3) 0.033 No AD",
                                              "(4) 0.12 No AD",
                                              "(3) 0.033 No AD",
                                              "(2) 0.12 AD",
                                              "(1) 0.033 AD",
                                              "(2) 0.12 AD",
                                              "(1) 0.033 AD"),
                                    value = c(mean((GetModelResponses_Texture %>% 
                                                      filter(airdrag == "NoAirdrag" & r == 0.12 & condsize == "Congruent"))$ErrorRatioSpace),
                                              mean((GetModelResponses_Texture %>% 
                                                      filter(airdrag == "NoAirdrag" & r == 0.033 & condsize == "Congruent"))$ErrorRatioSpace),
                                              mean((GetModelResponses_Texture %>% 
                                                      filter(airdrag == "NoAirdrag" & r == 0.12 & condsize == "Incongruent"))$ErrorRatioSpace),
                                              mean((GetModelResponses_Texture %>% 
                                                      filter(airdrag == "NoAirdrag" & r == 0.033 & condsize == "Incongruent"))$ErrorRatioSpace),
                                              mean((GetModelResponses_Texture %>% 
                                                      filter(airdrag == "Airdrag" & r == 0.12 & condsize == "Congruent"))$ErrorRatioSpace),
                                              mean((GetModelResponses_Texture %>% 
                                                      filter(airdrag == "Airdrag" & r == 0.033 & condsize == "Congruent"))$ErrorRatioSpace),
                                              mean((GetModelResponses_Texture %>% 
                                                      filter(airdrag == "Airdrag" & r == 0.12 & condsize == "Incongruent"))$ErrorRatioSpace),
                                              mean((GetModelResponses_Texture %>% 
                                                      filter(airdrag == "Airdrag" & r == 0.033 & condsize == "Incongruent"))$ErrorRatioSpace)),
                                    Congruent = c("Congruent","Congruent","Incongruent","Incongruent",
                                                  "Congruent","Congruent","Incongruent","Incongruent")
                                    )
ADModel_Space_Texture_Effect$Model = "Air Drag Model"
ADModel_Space_Texture_Effect$Modality = "Space"
ADModel_Space_Texture_Effect$Effect = "Effect"

#air drag model AND effect of texture
ADModel_Time_Texture_Effect = data.frame(label = c("(4) 0.12 No AD",
                                             "(3) 0.033 No AD",
                                             "(4) 0.12 No AD",
                                             "(3) 0.033 No AD",
                                             "(2) 0.12 AD",
                                             "(1) 0.033 AD",
                                             "(2) 0.12 AD",
                                             "(1) 0.033 AD"),
                                    value = c(mean((GetModelResponses_Texture %>% 
                                                      filter(airdrag == "NoAirdrag" & r == 0.12 & condsize == "Congruent"))$ErrorRatioTime),
                                              mean((GetModelResponses_Texture %>% 
                                                      filter(airdrag == "NoAirdrag" & r == 0.033 & condsize == "Congruent"))$ErrorRatioTime),
                                              mean((GetModelResponses_Texture %>% 
                                                      filter(airdrag == "NoAirdrag" & r == 0.12 & condsize == "Incongruent"))$ErrorRatioTime),
                                              mean((GetModelResponses_Texture %>% 
                                                      filter(airdrag == "NoAirdrag" & r == 0.033 & condsize == "Incongruent"))$ErrorRatioTime),
                                              mean((GetModelResponses_Texture %>% 
                                                      filter(airdrag == "Airdrag" & r == 0.12 & condsize == "Congruent"))$ErrorRatioTime),
                                              mean((GetModelResponses_Texture %>% 
                                                      filter(airdrag == "Airdrag" & r == 0.033 & condsize == "Congruent"))$ErrorRatioTime),
                                              mean((GetModelResponses_Texture %>% 
                                                      filter(airdrag == "Airdrag" & r == 0.12 & condsize == "Incongruent"))$ErrorRatioTime),
                                              mean((GetModelResponses_Texture %>% 
                                                      filter(airdrag == "Airdrag" & r == 0.033 & condsize == "Incongruent"))$ErrorRatioTime)),
                                   Congruent = c("Congruent","Congruent","Incongruent","Incongruent",
                                                 "Congruent","Congruent","Incongruent","Incongruent")
)

ADModel_Time_Texture_Effect$Model = "Air Drag Model"
ADModel_Time_Texture_Effect$Modality = "Time"
ADModel_Time_Texture_Effect$Effect = "Effect"

#no air drag model
GravityModel_Space_Texture = data.frame(label = c("(4) 0.12 No AD",
                                                   "(3) 0.033 No AD",
                                                   "(4) 0.12 No AD",
                                                   "(3) 0.033 No AD",
                                                   "(2) 0.12 AD",
                                                   "(1) 0.033 AD",
                                                   "(2) 0.12 AD",
                                                   "(1) 0.033 AD"),
                                          value = c(mean((GetModelResponses_NoAirdragModel %>% 
                                                            filter(airdrag == "NoAirdrag" & r == 0.12 & condsize == "Congruent"))$ErrorRatioSpace),
                                                    mean((GetModelResponses_NoAirdragModel %>% 
                                                            filter(airdrag == "NoAirdrag" & r == 0.033 & condsize == "Congruent"))$ErrorRatioSpace),
                                                    mean((GetModelResponses_NoAirdragModel %>% 
                                                            filter(airdrag == "NoAirdrag" & r == 0.12 & condsize == "Incongruent"))$ErrorRatioSpace),
                                                    mean((GetModelResponses_NoAirdragModel %>% 
                                                            filter(airdrag == "NoAirdrag" & r == 0.033 & condsize == "Incongruent"))$ErrorRatioSpace),
                                                    mean((GetModelResponses_NoAirdragModel %>% 
                                                            filter(airdrag == "Airdrag" & r == 0.12 & condsize == "Congruent"))$ErrorRatioSpace),
                                                    mean((GetModelResponses_NoAirdragModel %>% 
                                                            filter(airdrag == "Airdrag" & r == 0.033 & condsize == "Congruent"))$ErrorRatioSpace),
                                                    mean((GetModelResponses_NoAirdragModel %>% 
                                                            filter(airdrag == "Airdrag" & r == 0.12 & condsize == "Incongruent"))$ErrorRatioSpace),
                                                    mean((GetModelResponses_NoAirdragModel %>% 
                                                            filter(airdrag == "Airdrag" & r == 0.033 & condsize == "Incongruent"))$ErrorRatioSpace)),
                                         Congruent = c("Congruent","Congruent","Incongruent","Incongruent",
                                                       "Congruent","Congruent","Incongruent","Incongruent")
                                         )
GravityModel_Space_Texture$Model = "No Air Drag Model"
GravityModel_Space_Texture$Modality = "Space"
GravityModel_Space_Texture$Effect = "No Effect"

#no air drag model
GravityModel_Timing_Texture = data.frame(label = c("(4) 0.12 No AD",
                                                   "(3) 0.033 No AD",
                                                   "(4) 0.12 No AD",
                                                   "(3) 0.033 No AD",
                                                   "(2) 0.12 AD",
                                                   "(1) 0.033 AD",
                                                   "(2) 0.12 AD",
                                                   "(1) 0.033 AD"),
                                         value = c(mean((GetModelResponses_NoAirdragModel %>% 
                                                           filter(airdrag == "NoAirdrag" & r == 0.12 & condsize == "Congruent"))$ErrorRatioTime),
                                                   mean((GetModelResponses_NoAirdragModel %>% 
                                                           filter(airdrag == "NoAirdrag" & r == 0.033 & condsize == "Congruent"))$ErrorRatioTime),
                                                   mean((GetModelResponses_NoAirdragModel %>% 
                                                           filter(airdrag == "NoAirdrag" & r == 0.12 & condsize == "Incongruent"))$ErrorRatioTime),
                                                   mean((GetModelResponses_NoAirdragModel %>% 
                                                           filter(airdrag == "NoAirdrag" & r == 0.033 & condsize == "Incongruent"))$ErrorRatioTime),
                                                   mean((GetModelResponses_NoAirdragModel %>% 
                                                           filter(airdrag == "Airdrag" & r == 0.12 & condsize == "Congruent"))$ErrorRatioTime),
                                                   mean((GetModelResponses_NoAirdragModel %>% 
                                                           filter(airdrag == "Airdrag" & r == 0.033 & condsize == "Congruent"))$ErrorRatioTime),
                                                   mean((GetModelResponses_NoAirdragModel %>% 
                                                           filter(airdrag == "Airdrag" & r == 0.12 & condsize == "Incongruent"))$ErrorRatioTime),
                                                   mean((GetModelResponses_NoAirdragModel %>% 
                                                           filter(airdrag == "Airdrag" & r == 0.033 & condsize == "Incongruent"))$ErrorRatioTime)),
                                         Congruent = c("Congruent","Congruent","Incongruent","Incongruent",
                                                       "Congruent","Congruent","Incongruent","Incongruent")
                                        )


GravityModel_Timing_Texture$Model = "No Air Drag Model"
GravityModel_Timing_Texture$Modality = "Time"
GravityModel_Timing_Texture$Effect = "No Effect"

#air drag model, but extrapolation is correct, according to the actual size and mass
ADModel_Space_TextureNoEffect = data.frame(label = c("(4) 0.12 No AD",
                                                  "(3) 0.033 No AD",
                                                  "(4) 0.12 No AD",
                                                  "(3) 0.033 No AD",
                                                  "(2) 0.12 AD",
                                                  "(1) 0.033 AD",
                                                  "(2) 0.12 AD",
                                                  "(1) 0.033 AD"),
                                        value = c(mean((GetModelResponses_AirdragModel %>% 
                                                          filter(airdrag == "NoAirdrag" & r == 0.12 & condsize == "Congruent"))$ErrorRatioSpace),
                                                  mean((GetModelResponses_AirdragModel %>% 
                                                          filter(airdrag == "NoAirdrag" & r == 0.033 & condsize == "Congruent"))$ErrorRatioSpace),
                                                  mean((GetModelResponses_AirdragModel %>% 
                                                          filter(airdrag == "NoAirdrag" & r == 0.12 & condsize == "Incongruent"))$ErrorRatioSpace),
                                                  mean((GetModelResponses_AirdragModel %>% 
                                                          filter(airdrag == "NoAirdrag" & r == 0.033 & condsize == "Incongruent"))$ErrorRatioSpace),
                                                  mean((GetModelResponses_AirdragModel %>% 
                                                          filter(airdrag == "Airdrag" & r == 0.12 & condsize == "Congruent"))$ErrorRatioSpace),
                                                  mean((GetModelResponses_AirdragModel %>% 
                                                          filter(airdrag == "Airdrag" & r == 0.033 & condsize == "Congruent"))$ErrorRatioSpace),
                                                  mean((GetModelResponses_AirdragModel %>% 
                                                          filter(airdrag == "Airdrag" & r == 0.12 & condsize == "Incongruent"))$ErrorRatioSpace),
                                                  mean((GetModelResponses_AirdragModel %>% 
                                                          filter(airdrag == "Airdrag" & r == 0.033 & condsize == "Incongruent"))$ErrorRatioSpace)),
                                        Congruent = c("Congruent","Congruent","Incongruent","Incongruent",
                                                      "Congruent","Congruent","Incongruent","Incongruent")
                                        )
ADModel_Space_TextureNoEffect$Model = "Air Drag Model"
ADModel_Space_TextureNoEffect$Modality = "Space"
ADModel_Space_TextureNoEffect$Effect = "No Effect"

#air drag model, but extrapolation is correct, according to the actual size and mass
ADModel_Timing_TextureNoEffect = data.frame(label = c("(4) 0.12 No AD",
                                                      "(3) 0.033 No AD",
                                                      "(4) 0.12 No AD",
                                                      "(3) 0.033 No AD",
                                                      "(2) 0.12 AD",
                                                      "(1) 0.033 AD",
                                                      "(2) 0.12 AD",
                                                      "(1) 0.033 AD"),
                                            value = c(mean((GetModelResponses_AirdragModel %>% 
                                                              filter(airdrag == "NoAirdrag" & r == 0.12 & condsize == "Congruent"))$ErrorRatioTime),
                                                      mean((GetModelResponses_AirdragModel %>% 
                                                              filter(airdrag == "NoAirdrag" & r == 0.033 & condsize == "Congruent"))$ErrorRatioTime),
                                                      mean((GetModelResponses_AirdragModel %>% 
                                                              filter(airdrag == "NoAirdrag" & r == 0.12 & condsize == "Incongruent"))$ErrorRatioTime),
                                                      mean((GetModelResponses_AirdragModel %>% 
                                                              filter(airdrag == "NoAirdrag" & r == 0.033 & condsize == "Incongruent"))$ErrorRatioTime),
                                                      mean((GetModelResponses_AirdragModel %>% 
                                                              filter(airdrag == "Airdrag" & r == 0.12 & condsize == "Congruent"))$ErrorRatioTime),
                                                      mean((GetModelResponses_AirdragModel %>% 
                                                              filter(airdrag == "Airdrag" & r == 0.033 & condsize == "Congruent"))$ErrorRatioTime),
                                                      mean((GetModelResponses_AirdragModel %>% 
                                                              filter(airdrag == "Airdrag" & r == 0.12 & condsize == "Incongruent"))$ErrorRatioTime),
                                                      mean((GetModelResponses_AirdragModel %>% 
                                                              filter(airdrag == "Airdrag" & r == 0.033 & condsize == "Incongruent"))$ErrorRatioTime)),
                                            Congruent = c("Congruent","Congruent","Incongruent","Incongruent",
                                                          "Congruent","Congruent","Incongruent","Incongruent")
                                            )
ADModel_Timing_TextureNoEffect$Model = "Air Drag Model"
ADModel_Timing_TextureNoEffect$Modality = "Time"
ADModel_Timing_TextureNoEffect$Effect = "No Effect"

###put all of the above dataframes together for a full dataframe with all relevant (?) values for Hypothesis 3
DF_Predictions_H3 = rbind(ADModel_Space_Texture_Effect,
                          ADModel_Time_Texture_Effect,
                          GravityModel_Space_Texture,
                          GravityModel_Timing_Texture,
                          ADModel_Space_TextureNoEffect,
                          ADModel_Timing_TextureNoEffect)

#Plot Hypothesis 1
Hypothesis1 = ggplot(DF_Predictions_H1,aes(label,value)) +
  geom_point(size = 5) +
  facet_grid(Modality~Model) +
  ylab("Error Ratio") +
  scale_x_discrete(labels = c("AD:Pres.",
                              "AD:Abs."),
                   name = "") +
  ggtitle("A. Hypothesis 1") +
  theme(axis.text.x= element_text(size = rel(0.75)),
        axis.title.x = element_text(size = rel(0.75))) +
  geom_hline(yintercept = 1,linetype = 2)

#Plot Hypothesis 2
Hypothesis2 = ggplot(DF_Predictions_H2, aes(label,value,color = Condition)) +
  geom_point(size = 5) +
  facet_grid(Modality~Model) +
  scale_x_discrete(labels = c("0.033 m\nAD:Pres.",
                              "0.12 m\nAD:Pres.",
                              "0.033 m\nAD:Abs.",
                              "0.12 m\nAD:Abs."),
                   name = "") +
  scale_color_manual(values =  c(BlauUB,LightBlauUB),
                     name = "") +
  ylab("Error Ratio") +
  theme(axis.text.x= element_text(size = rel(0.75)),
        axis.title.x = element_text(size = rel(0.75)),
        legend.position = "bottom") +
  ggtitle("B. Hypothesis 2") +
  geom_hline(yintercept = 1,linetype = 2)

#Plot Hypothesis 3
Hypothesis3 = ggplot(DF_Predictions_H3,
       aes(label,value,color = Effect,shape=Congruent)) +
  geom_point(size = 5) +
  facet_grid(Modality~Model) +
    scale_x_discrete(labels = c("0.033 m\nAD:Pres.",
                                "0.12 m\nAD:Pres.",
                                "0.033 m\nAD:Abs.",
                                "0.12 m\nAD:Abs."),
                     name = "") +
  scale_color_manual(values =  c(Red,LightRed),
                     name = "") +
  scale_shape_discrete(name = "") +
  ylab("Error Ratio") +
  xlab("") +
  theme(legend.position = "bottom",
        axis.text.x= element_text(size = rel(0.75)),
        axis.title.x = element_text(size = rel(0.75))) +
  geom_hline(yintercept = 1,linetype = 2) +
  ggtitle("C. Hypothesis 3")

#everything together
H123 = plot_grid(Hypothesis1, Hypothesis2, Hypothesis3, rel_widths = c(0.7,1,1), nrow=1)
ggsave("Figures/AllHypotheses.jpg",w = 16, h = 6)

