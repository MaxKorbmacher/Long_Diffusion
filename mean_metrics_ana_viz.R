# analysis and visualisation of whole-brain average values

########### PREP ####

# load packages
library(dplyr)
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

# load data (without DRAD extra outliers removed as there were many outliers)
T1 = read.csv("/cluster/projects/p33/users/maxk/UKB/longitudinal/data/T1_noDRAD.csv")
T2 = read.csv("/cluster/projects/p33/users/maxk/UKB/longitudinal/data/T2_noDRAD.csv")
# hence, we have to exclude DRAD extra from our analysis, and focus on global average scores here
T1 = T1 %>% dplyr::select(contains("mean"), age, sex, site_t3, eid) %>% dplyr::select(!Drad_extra_Mean)
T2 = T2 %>% dplyr::select(contains("mean"), age, sex, site_t3, eid) %>% dplyr::select(!Drad_extra_Mean)
# 
# # standardise each diffusion metrics (Z-scoring) for later comparison
# for (i in 1:27){
#   T1[,i] = (T1[,i]-(mean(T1[,i])))/(sd(T1[,i]))
#   T2[,i] = (T2[,i]-(mean(T2[,i])))/(sd(T2[,i]))
# }
#
# create time point dummy TP
T1$TP = replicate(nrow(T1),0)
T2$TP = replicate(nrow(T1),1)
T2$time = T2$age - T1$age
T1$time = 0
# make single df
data = rbind(T1, T2)

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

# make list for all the models
models = list()
for (i in outcome_vars){
  f = formula(paste(i,"~ age*sex + time + (1|site_t3) + (1|eid)"))
  models[[i]] = lmer(f, data=data)
}


### PREDICT ####

# now, we predict from those models and make a list (predictions)containing the predictions for each metric 
predictions = list()
age = data$age
ID = data$eid
TP = data$TP
data1 = data.frame(ID, age, TP)
for (i in 1:length(outcome_vars)){
  data1$fit = predict(models[[i]])
  predictions[[i]] = data1
}

# give the prediction list the right names
names(predictions) = outcome_vars

# make a list of names for the y-Axis
y_lab = c("BRIA - V intra", "BRIA - V extra", "BRIA - V csf", "BRIA - micro RD","BRIA - micro FA",
          "BRIA - micro AX", "BRIA - micro ADC", "BRIA - DAX intra", "BRIA - DAX extra",
          
          "DKI - RK", "DKI - AK", "DKI -MK",
          
          "DTI - FA", "DTI - MD", "DTI - RD", "DTI - AD",
          
          "SMT - FA", "SMT - long", "SMT - MD", "SMT - trans",
          
          "SMTmc - intra","SMTmc - extra MD", "SMTmc - extra trans",  "SMTmc - diff",
          
          "WMTI - axEAD", "WMTI - AWF", "WMTI - radEAD")


######## PLOTTING ####

# SCATTER ####
# now, plot each of the predictions in a scatter plot with two time points and a fitted line
lm_plot = list()
for (i in 1:length(y_lab)){
  lm_plot[[i]] = ggplot(data=predictions[[i]],aes(x=Age, y=fit))+ 
  geom_line(aes(y=fit,group=ID),alpha=.4,lwd=.5) + #,col=TP
  geom_point(aes(y=fit,shape=factor(TP), col=TP),alpha=.6) +#
  #geom_ribbon(aes(group=TP,ymin=lwr,ymax=upr),alpha=0.3,fill="grey") + 
  geom_smooth(data=predictions[[i]], colour="red", aes(y=fit),lwd=1.1, method = "lm") +
  xlab("Age")+
  ylab (paste(y_lab[i]))+ 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  theme(legend.position = "none")
}
# we can do the same for a smooth
smooth_plot = list()
for (i in 1:length(y_lab)){
  smooth_plot[[i]] = ggplot(data=predictions[[i]],aes(x=Age, y=fit))+ 
    geom_line(aes(y=fit,group=ID),alpha=.4,lwd=.5) + #,col=TP
    geom_point(aes(y=fit,shape=factor(TP), col=TP),alpha=.6) +#
    geom_smooth(data=predictions[[i]], colour="red", aes(y=fit),lwd=1.1) +
    xlab("Age")+
    ylab (paste(y_lab[i]))+ 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    theme(legend.position = "none")
}
lm_plots = do.call("grid.arrange", c(lm_plot,ncol = 3))
smooth_plots = do.call("grid.arrange", c(smooth_plot,ncol = 3))
ggsave("/cluster/projects/p33/users/maxk/UKB/export/smooth_ageing_plots.pdf",smooth_plots, width = 15, height = 20)
ggsave("/cluster/projects/p33/users/maxk/UKB/export/lm_ageing_plots.pdf",lm_plots, width = 15, height = 20)

# For visual inspection on the cluster:
#ggsave("/tsd/p33/home/p33-maxk/test/test2.pdf",test2, width = 15, height = 20)

# inspect the effect of time point in the model
# model inspections: print t-values for the effect of time point in the model
tval = c()
for (i in 1:length(models)){
  tval[i] = summary(models[[i]])$coefficients[4,3]
}
sig = c()
for (i in 1:length(models)){
  sig[i] = ifelse(abs(tval[i]) > 1.96, "*", "non-sig")
}
vals = data.frame(y_lab, tval,sig)
vals

# VIOLINS ####
# we re-run the mixed models described above without controlling for time point, since we are interested in the difference

# make list for all the models, but this time WITHOUT including time point & we use the original data
models = list()
for (i in outcome_vars){
  f = formula(paste(i,"~ sex + (1|site_t3) + (1|eid)"))
  models[[i]] = lmer(f, data=df_orig)
}

# now, we predict from those models and make a list (predictions)containing the predictions for each metric 
predictions = list()
age = df_orig$age
ID = df_orig$eid
TP = df_orig$TP
data1 = data.frame(ID, age, TP)
for (i in 1:length(outcome_vars)){
  data1$fit = predict(models[[i]])
  predictions[[i]] = data1
}

# give the prediction list the right names
names(predictions) = outcome_vars

# calculate p-values for paired-samples t-tests (seemingly not sensitive enough: no TP difference detected)
result = c()
for (i in 1:length(y_lab)){
  result = t.test(fit ~ TP, data = predictions[[i]])$p.value
  result = signif(result, digits = 3)
  print(result)
}

# plot time point differences
plot2 = list()
for (i in 1:length(y_lab)){
  plot2[[i]] = ggplot(predictions[[i]], aes(x=factor(TP), y=fit, fill = factor(TP))) + 
    geom_violin() + xlab("Time Point") + ylab(paste(y_lab[i])) + 
    theme_classic()+scale_fill_manual(values=c("#FFB531","#BC211A"))+ 
    stat_summary(fun.y=mean,  geom="point", color="black")+ 
    theme(legend.position="none")+ 
    theme(aspect.ratio=1) +
    add_pvalue()
}
test3 =  do.call("grid.arrange", c(plot2,ncol = 3))
# print or safe
#test3
#ggsave("/tsd/p33/home/p33-maxk/test/test3.pdf",test3, width = 15, height = 20)
