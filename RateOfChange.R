# rate of change in white matter metrics
# Max Korbmacher (max.korbmacher@gmail.com), 19.07.2023

## 1) PREPARATIONS ####
#### 1.1) PREP ENV ####

# load packages and install if not already installed
if (!require("pacman")) install.packages("pacman")
pacman::p_load(lme4, nlme, ggplot2, tidyverse, lm.beta, remotes, ggpubr, grid, lmtest, car, lmtest,lmeInfo,lmerTest,sjstats,effsize, sjPlot, gridExtra)

#### 1.2) PREP DATA ####
# set working directory where the data files are situated
setwd("/home/max/Documents/Projects/Diffusion/UKBlong/data_export/")
# load data (without DRAD extra outliers removed as there were many outliers)
T1 = read.csv("T1_noDRAD.csv")
T2 = read.csv("T2_noDRAD.csv")
# we have to exclude DRAD extra from our analysis, as this metric did poor in the QC
T1 = T1 %>% dplyr::select(contains("mean"), age, sex, site_t3, eid) %>% dplyr::select(!Drad_extra_Mean)
T2 = T2 %>% dplyr::select(contains("mean"), age, sex, site_t3, eid) %>% dplyr::select(!Drad_extra_Mean)

# now, we can substract one time point from the other
ROC = T2[1:28]-T1[1:28]
lm_list = list()
age = T1$age
sex = T1$sex
site = T1$site_t3
for (i in 1:27){
  tmpROC = scale(ROC[i]/abs(ROC[28]), center = FALSE)
  lm_list[[i]] = data.frame(age, sex, site, tmpROC)
  colnames(lm_list[[i]]) = c("age","sex","site","ROC")
}
# make name list of outcome variables in as in data frame
outcome_vars = colnames(ROC[1:27])

# make name labels for the outcome_vars (more presentable)
var_labels = c("BRIA - V intra", "BRIA - V extra", "BRIA - V csf", "BRIA - micro RD","BRIA - micro FA",
               "BRIA - micro AX", "BRIA - micro ADC", "BRIA - DAX intra", "BRIA - DAX extra",
               
               "DKI - RK", "DKI - AK", "DKI -MK",
               
               "DTI - FA", "DTI - MD", "DTI - RD", "DTI - AD",
               
               "SMT - FA", "SMT - long", "SMT - MD", "SMT - trans",
               
               "SMTmc - intra","SMTmc - extra MD", "SMTmc - extra trans",  "SMTmc - diff",
               
               "WMTI - axEAD", "WMTI - AWF", "WMTI - radEAD")

## 2) Global Annual Change Plots ####
lm_plot = list()
for (i in 1:length(var_labels)){
  lm_plot[[i]] = ggplot(data=lm_list[[i]],aes(x=age, y=ROC))+ 
    #geom_line(aes(y=fit,group=ID),alpha=.4,lwd=.5) + #,col=TP
    #geom_point(aes(y=fit,shape=factor(TP), col=TP),alpha=.6) +#
    #geom_ribbon(aes(group=TP,ymin=lwr,ymax=upr),alpha=0.3,fill="grey") + 
    #geom_smooth(method = "lm", color= "#E69F00") +
    geom_smooth(method = "gam", color = "#56B4E9") +
    geom_line(aes(y=0),linetype="dotted")+
    # stat_regline_equation(
    #   aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
    #   formula = formula(paste("data"[i]"~age")), data =data) +
    #geom_smooth(method = "lm") +
    #stat_regline_equation() +
    #stat_cor()+
    xlab("Age")+
    ylab (paste(var_labels[i]))+ 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    theme(legend.position = "none")+
    stat_cor(label.y = 1)
    #labs(title = (paste("linear R2 = ",round(lm_r.sq1[i],digits = 4), "      ","non-linear R2 = ", round(gam_r.sq1[i],digits = 4))))
}
plot= do.call("grid.arrange", c(lm_plot,ncol = 3))
ggsave("/home/max/Documents/Projects/Diffusion/UKBlong/Results/annual_change_plots.pdf",plot, width = 15, height = 20)

## 3) Regional Annual Changes ####
# load data (without DRAD extra outliers removed as there were many outliers)
T1 = read.csv("T1_noDRAD.csv")
T2 = read.csv("T2_noDRAD.csv")

# now, we can substract one time point from the other
ROC = T2[2:1934]-T1[4:1936]
ROC = ROC %>% select(-eid)
lm_list = list()
age = T1$age
sex = T1$sex
site = T1$site_t3
age_diff = T2$age - T1$age
for (i in 1:1932){
  tmpROC = scale(ROC[[i]]/age_diff, center = FALSE)
  lm_list[[i]] = data.frame(age, sex, site, tmpROC)
  colnames(lm_list[[i]]) = c("age","sex","site","ROC")
}

# the many scatter plots visualisation is messy, we summarise the findings using correlations between age and absolute annual change:
cors = c()
for (i in 1:1932){
  cors[i]=cor(abs(lm_list[[i]]$ROC), lm_list[[i]]$age)
}

cor_df = data.frame(col_names,cors)

# the relationship between age and absolute annual change can be summarized in a density plot
plot(density(cors))

# across all metrics the absolute annual rate of change - age relationship is > 0.
t.test(cor_df$cors,mu=0)
# however, this does not tell us whether the rate of change trends towards 0, or away from it.
mean(cor_df$cors)
sd(cor_df$cors)


## 4) Fornix Annual Changes ####
# load data (without DRAD extra outliers removed as there were many outliers)
# T1 = read.csv("T1_noDRAD.csv")
# T2 = read.csv("T2_noDRAD.csv")
# we have to exclude DRAD extra from our analysis, as this metric did poor in the QC
T1 = T1 %>% dplyr::select(contains("fornix"), age, sex, site_t3, eid)# %>% dplyr::select(!Drad_extra_Mean)
T2 = T2 %>% dplyr::select(contains("fornix"), age, sex, site_t3, eid)# %>% dplyr::select(!Drad_extra_Mean)

# now, we can substract one time point from the other
ROC = T2[1:85]-T1[1:85]
lm_list = list()
age = T1$age
sex = T1$sex
site = T1$site_t3
for (i in 1:84){
  tmpROC = scale(ROC[i]/abs(ROC[85]), center = FALSE)
  lm_list[[i]] = data.frame(age, sex, site, tmpROC)
  colnames(lm_list[[i]]) = c("age","sex","site","ROC")
}
lm_plot = list()
col_names = str_replace(col_names, "Striaterminalis", "St")
for (i in 1:84){
  lm_plot[[i]] = ggplot(data=lm_list[[i]],aes(x=age, y=ROC))+ 
    #geom_line(aes(y=fit,group=ID),alpha=.4,lwd=.5) + #,col=TP
    #geom_point(aes(y=fit,shape=factor(TP), col=TP),alpha=.6) +#
    #geom_ribbon(aes(group=TP,ymin=lwr,ymax=upr),alpha=0.3,fill="grey") + 
    #geom_smooth(method = "lm", color= "#E69F00") +
    geom_smooth(method = "gam", color = "#56B4E9") +
    geom_line(aes(y=0),linetype="dotted")+
    # stat_regline_equation(
    #   aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
    #   formula = formula(paste("data"[i]"~age")), data =data) +
    #geom_smooth(method = "lm") +
    #stat_regline_equation() +
    #stat_cor()+
    xlab("Age")+
    ylab (paste(col_names[i]))+ 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    theme(legend.position = "none")
    #stat_cor(label.y = 1)
  #labs(title = (paste("linear R2 = ",round(lm_r.sq1[i],digits = 4), "      ","non-linear R2 = ", round(gam_r.sq1[i],digits = 4))))
}
plot= do.call("grid.arrange", c(lm_plot))
ggsave("/home/max/Documents/Projects/Diffusion/UKBlong/Results/FORNIX_annual_change_plots.pdf",plot, width = 20, height = 20)

# the many scatter plots visualisation is messy, we summarise the findings using correlations:
cors = c()
for (i in 1:84){
  cors[i]=cor(abs(lm_list[[i]]$ROC), lm_list[[i]]$age)
}

cor_df = data.frame(col_names,cors)
cor_df %>% filter(abs(cors)>0.1)
plot(density(cors))

# across all metrics in the fornix, age associations of the annual rate of change > 0
t.test(cor_df$cors,mu=0)

# we test this for all data as well.
