# analysis and visualisation of whole-brain average values

########### PREP ####

# load packages
library(dplyr)
library(ggpubr)
library(lme4)
library(ggplot2)
#install.packages("merTools")
library(merTools)
library(gridExtra)
#install.packages("ggprism")
library(ggprism)
library(patchwork)
library(magrittr)
#install.packages("olsrr")
library(olsrr)
library(ggeffects)
library(nlme)
#install.packages("ggpmisc")
library(ggpmisc)
# set working directory where the data files are situated
setwd("/home/max/Documents/Projects/Diffusion/UKBlong/data_export/")
# load data (without DRAD extra outliers removed as there were many outliers)
T1 = read.csv("T1_noDRAD.csv")
T2 = read.csv("T2_noDRAD.csv")
# hence, we have to exclude DRAD extra from our analysis, and focus on global average scores here
T1 = T1 %>% dplyr::select(contains("mean"), age, sex, site_t3, eid) %>% dplyr::select(!Drad_extra_Mean)
T2 = T2 %>% dplyr::select(contains("mean"), age, sex, site_t3, eid) %>% dplyr::select(!Drad_extra_Mean)
# 
# standardise each diffusion metrics (Z-scoring) for later comparison
for (i in 1:27){
  T1[,i] = (T1[,i]-(mean(T1[,i])))/(sd(T1[,i]))
  T2[,i] = (T2[,i]-(mean(T2[,i])))/(sd(T2[,i]))
}

# create time point dummy TP
T1$TP = replicate(nrow(T1),0)
T2$TP = replicate(nrow(T1),1)
T2$time = T2$age - T1$age
T1$time = 0
# make single df
data = rbind(T1, T2)
data$TP = factor(data$TP)
data$sex = factor(data$sex)
data$eid = factor(data$eid)
data$site_t3 = factor(data$site_t3)
data$age = as.numeric(data$age)
# columns which are supposed to be numeric
i = c(1:27)
# then make them numeric
data[ , i] <- apply(data[ , i], 2,          
                    function(x) as.numeric(as.character(x)))
# build a mixed linear effects model & loop over all the dependent vars being the diffusion metrics
## (we have 27 diffusion whole-brain average metrics, which sit in the first 1-27 columns)
outcome_vars = names(data)[1:27]

# we make a copy of the original data for later use
df_orig = data

# test for heteroscedasticity
# names1 = names(data)
# het_list = list()
# for (i in 1:length(names1)){
#   f = formula(paste(names1[i],"~age*sex + age*AgeDiff + site_t3"))
#   model = lm(f,data=data)
#   het_list[[i]] = ols_test_breusch_pagan(model)
#   print(het_list[[i]]$p)
# }
# >> most metrics are homeoscedastic in the model that we plan to use. 
# make list for all the models (linear and non-linear)
# linear_models = list()
# #non_linear_models = list()
# for (i in outcome_vars){
#   f = formula(paste(i,"~ age*sex + TP + (1|site_t3)+(1|eid)"))
#   linear_models[[i]] = lmer(f, data=data)
#   #non_linear_models[[i]] = nlmer(f, data=data)
# }
# 
# # now, we predict from those models and make a list (predictions)containing the predictions for each metric 
# linear_predictions = list()
# non_linear_predictions = list()
# age = data$age
# ID = data$eid
# TP = data$TP
# data1 = data.frame(ID, age, TP)
# data2 = data1
# for (i in 1:length(outcome_vars)){
#   data1$fit = predict(linear_models[[i]])
#   #data2$fit = predict(non_linear_models[[i]])
#   linear_predictions[[i]] = data1
#   #non_linear_predictions[[i]] = data2
# }

# give the prediction list the right names
names(linear_predictions) = outcome_vars
#names(non_linear_predictions) = outcome_vars

# make a list of names for the y-Axis
y_lab = c("BRIA - V intra", "BRIA - V extra", "BRIA - V csf", "BRIA - micro RD","BRIA - micro FA",
          "BRIA - micro AX", "BRIA - micro ADC", "BRIA - DAX intra", "BRIA - DAX extra",
          "DKI - RK", "DKI - AK", "DKI -MK",
          "DTI - FA", "DTI - MD", "DTI - RD", "DTI - AD",
          "SMT - FA", "SMT - long", "SMT - MD", "SMT - trans",
          "SMTmc - intra","SMTmc - extra MD", "SMTmc - extra trans",  "SMTmc - diff",
          "WMTI - axEAD", "WMTI - AWF", "WMTI - radEAD")

######## PLOTTING ####
# first run gams and lms for each time point to get R squared values
lm_r.sq1 = c()
lm_r.sq2 = c()
gam_r.sq1 = c()
gam_r.sq2 = c()

for (i in 1:27){
  T1_lm_model = lm(T1[[i]]~(age)+age:sex+sex, data = T1)
  lm_r.sq1[i] = summary(T1_lm_model)$r.sq
  T2_lm_model = lm(T2[[i]]~(age)+age:sex+sex, data = T2)
  lm_r.sq2[i] = summary(T2_lm_model)$r.sq
  T1_gam_model = gam(T1[[i]]~s(age)+age:sex+sex, data = T1)
  gam_r.sq1[i] = summary(T1_gam_model)$r.sq
  T2_gam_model = gam(T2[[i]]~s(age)+age:sex+sex, data = T2)
  gam_r.sq2[i] = summary(T2_gam_model)$r.sq
}
# #for an overview of the values:
# data.frame(lm_r.sq1,lm_r.sq2,gam_r.sq1,gam_r.sq2)

# SCATTER PLOT OF PREDICTED VALUES >> DIFFERENTIATING BETWEEN TIME POINTS ####
# now, plot each of the predictions in a scatter plot with two time points and a fitted line
lm_plot = list()
for (i in 1:length(y_lab)){
  lm_plot[[i]] = ggplot(data=linear_predictions[[i]],aes(x=age, y=fit,color = TP))+ 
  #geom_line(aes(y=fit,group=ID),alpha=.4,lwd=.5) + #,col=TP
  #geom_point(aes(y=fit,shape=factor(TP), col=TP),alpha=.6) +#
  #geom_ribbon(aes(group=TP,ymin=lwr,ymax=upr),alpha=0.3,fill="grey") + 
  geom_smooth(method = "lm") +
    # stat_regline_equation(
    #   aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
    #   formula = formula(paste("data"[i]"~age")), data =data) +
    #geom_smooth(method = "lm") +
    #stat_regline_equation() +
    #stat_cor()+
  xlab("Age")+
  ylab (paste(y_lab[i]))+ 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  theme(legend.position = "none")+
    labs(title = paste("TP1 R2 = ",round(lm_r.sq1[i],digits = 4), "      ","TP2 R2 = ", round(lm_r.sq2[i],digits = 4)))
}

# we can do the same for non-linear models
smooth_plot = list()
for (i in 1:length(y_lab)){
  smooth_plot[[i]] = ggplot(data=linear_predictions[[i]],aes(x=age, y=fit,color = TP))+ 
    geom_smooth(method = "gam") +
    xlab("Age")+
    ylab (paste(y_lab[i]))+ 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    theme(legend.position = "none") +
    labs(title = paste("TP1 R2 = ",round(gam_r.sq1[i],digits = 4), "      ","TP2 R2 = ", round(gam_r.sq2[i],digits = 4)))
}

lm_plots = do.call("grid.arrange", c(lm_plot,ncol = 3))
smooth_plots = do.call("grid.arrange", c(smooth_plot,ncol = 3))
ggsave("/home/max/Documents/Projects/Diffusion/UKBlong/Results/ADJ_smooth_ageing_plots.pdf",smooth_plots, width = 15, height = 20)
ggsave("/home/max/Documents/Projects/Diffusion/UKBlong/Results/ADJ_lm_ageing_plots.pdf",lm_plots, width = 15, height = 20)

# SCATTER PLOT OF PREDICTED VALUES >> NOT DIFFERENTIATING BETWEEN TIME POINTS ####
# now, plot each of the predictions in a scatter plot with two time points and a fitted line
lm_plot = list()
for (i in 1:length(y_lab)){
  lm_plot[[i]] = ggplot(data=linear_predictions[[i]],aes(x=age, y=fit,color = TP))+ 
    #geom_line(aes(y=fit,group=ID),alpha=.4,lwd=.5) + #,col=TP
    #geom_point(aes(y=fit,shape=factor(TP), col=TP),alpha=.6) +#
    #geom_ribbon(aes(group=TP,ymin=lwr,ymax=upr),alpha=0.3,fill="grey") + 
    geom_smooth(method = "lm") +
    # stat_regline_equation(
    #   aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
    #   formula = formula(paste("data"[i]"~age")), data =data) +
    #geom_smooth(method = "lm") +
    #stat_regline_equation() +
    #stat_cor()+
    xlab("Age")+
    ylab (paste(y_lab[i]))+ 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    theme(legend.position = "none")+
    labs(title = paste("TP1 R2 = ",round(lm_r.sq1[i],digits = 4), "      ","TP2 R2 = ", round(lm_r.sq2[i],digits = 4)))
}

# we can do the same for non-linear models
smooth_plot = list()
for (i in 1:length(y_lab)){
  smooth_plot[[i]] = ggplot(data=linear_predictions[[i]],aes(x=age, y=fit,color = TP))+ 
    geom_smooth(method = "gam") +
    xlab("Age")+
    ylab (paste(y_lab[i]))+ 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    theme(legend.position = "none") +
    labs(title = paste("TP1 R2 = ",round(gam_r.sq1[i],digits = 4), "      ","TP2 R2 = ", round(gam_r.sq2[i],digits = 4)))
}

lm_plots = do.call("grid.arrange", c(lm_plot,ncol = 3))
smooth_plots = do.call("grid.arrange", c(smooth_plot,ncol = 3))
ggsave("/home/max/Documents/Projects/Diffusion/UKBlong/Results/ONE_LINE_ADJ_smooth_ageing_plots.pdf",smooth_plots, width = 15, height = 20)
ggsave("/home/max/Documents/Projects/Diffusion/UKBlong/Results/ONE_LINE_ADJ_lm_ageing_plots.pdf",lm_plots, width = 15, height = 20)
