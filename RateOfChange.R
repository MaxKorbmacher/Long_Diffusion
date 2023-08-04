# rate of change in white matter metrics
# Max Korbmacher (max.korbmacher@gmail.com), 19.07.2023

## 1) PREPARATIONS ####
#### 1.1) PREP ENV ####

# load packages and install if not already installed
if (!require("pacman")) install.packages("pacman")
pacman::p_load(lme4, nlme, ggplot2, tidyverse, lm.beta, remotes, ggpubr, grid, lmtest, car, lmtest,lmeInfo,lmerTest,sjstats,effsize, sjPlot, gridExtra,ggsignif,vioplot)

#### 1.2) PREP DATA ####
# set working directory where the data files are situated
setwd("/home/max/Documents/Projects/Diffusion/UKBlong/data_export/")
# load data (without DRAD extra outliers removed as there were many outliers)
T1 = read.csv("T1_noDRAD.csv")
T2 = read.csv("T2_noDRAD.csv")
outliers = read.csv("/home/max/Documents/Projects/Diffusion/UKBlong/data_export/outliers.csv")

# exclude outliers / impossible values
T1 = T1[!T1$eid %in% outliers,]
T2 = T2[!T2$eid %in% outliers,]

# we have to exclude DRAD extra from our analysis, as this metric did poor in the QC on whole skeleton average values
# the same QC concern is true for axEAD estimated in the Fornix. We hence exclude the metric.
# we also extract only mean / whole skeleton average values in this step of the analysis
T1 = T1 %>% dplyr::select(contains("mean"), age, sex, site_t3, eid) %>% dplyr::select(!(contains("Drad_extra")))%>% dplyr::select(!(contains("axEAD")))
T2 = T2 %>% dplyr::select(contains("mean"), age, sex, site_t3, eid) %>% dplyr::select(!(contains("Drad_extra")))%>% dplyr::select(!(contains("axEAD")))

## 2) Age-stratified Global Annual Change ####

# now, we can subtract one time point from the other, the last column is age, for which we estimate the delta
ROC = T2[1:27]-T1[1:27]
lm_list = list()
age = T1$age
sex = T1$sex
site = T1$site_t3
for (i in 1:26){
  tmpROC = scale(ROC[i]/abs(ROC[27]), center = FALSE)
  lm_list[[i]] = data.frame(age, sex, site, tmpROC)
  lm_list[[i]]$age_group = findInterval(age, c(0, 55, 65, 75, 200))
  colnames(lm_list[[i]]) = c("age","sex","site","ROC", "group")
  lm_list[[i]]$group = factor(lm_list[[i]]$group)
  levels(lm_list[[i]]$group) = c(">55","55-65","65-75","75+")
  lm_list[[i]]$sex = factor(lm_list[[i]]$sex)
  levels(lm_list[[i]]$sex) = c("Female","Male")
}
# we make a list of all the metric names
# make name labels for the outcome_vars (more presentable)
var_labels = c("BRIA - V intra", "BRIA - V extra", "BRIA - V csf", "BRIA - micro RD","BRIA - micro FA",
               "BRIA - micro AX", "BRIA - micro ADC", "BRIA - DAX intra", "BRIA - DAX extra",
               "DKI - RK", "DKI - AK", "DKI -MK",
               "DTI - FA", "DTI - MD", "DTI - RD", "DTI - AD",
               "SMT - FA", "SMT - long", "SMT - MD", "SMT - trans",
               "SMTmc - intra","SMTmc - extra MD", "SMTmc - extra trans",  "SMTmc - diff",
               "WMTI - AWF", "WMTI - radEAD")
# then create an empty list to be filled with plots
plotlist = list()
# and loop over all the metrics
for (i in 1:length(lm_list)){
  # estimate mean by group
  mean <- lm_list[[i]] %>% 
    group_by(group) %>% 
    summarize(average = mean(ROC)) %>%
    ungroup()
  # plot age stratified rate of change
  plotlist[[i]] = ggplot2::ggplot(lm_list[[i]], aes(x = group, y=ROC, fill = group)) +
    ggplot2::geom_violin(alpha = 4/10) +   geom_boxplot(width=0.1,   outlier.shape = 0,outlier.size = 0,outlier.stroke = 0)+
    geom_signif(comparisons = list(c(">55", "55-65")), map_signif_level=c((0.001/(6*length(lm_list))),0.01/(6*length(lm_list)),0.05/(6*length(lm_list))), y_position = 5, tip_length = .01, vjust = 0.2, color="black") +
    geom_signif(comparisons = list(c(">55", "65-75")), map_signif_level=c((0.001/(6*length(lm_list))),0.01/(6*length(lm_list)),0.05/(6*length(lm_list))), y_position = 6, tip_length = .01, vjust = 0.2, color="black") +
    geom_signif(comparisons = list(c(">55", "75+")), map_signif_level=c((0.001/(6*length(lm_list))),0.01/(6*length(lm_list)),0.05/(6*length(lm_list))), y_position = 8, tip_length = .01, vjust = 0.2, color="black") +
    geom_signif(comparisons = list(c("55-65", "65-75")), map_signif_level=c((0.001/(6*length(lm_list))),0.01/(6*length(lm_list)),0.05/(6*length(lm_list))), y_position = 4.8, tip_length = .01, vjust = 0.2, color="black") +
    geom_signif(comparisons = list(c("65-75", "75+")), map_signif_level=c((0.001/(6*length(lm_list))),0.01/(6*length(lm_list)),0.05/(6*length(lm_list))), y_position = 4.6, tip_length = .01, vjust = 0.2, color="black") +
    geom_signif(comparisons = list(c("55-65", "75+")), map_signif_level=c((0.001/(6*length(lm_list))),0.01/(6*length(lm_list)),0.05/(6*length(lm_list))), y_position = 7, tip_length = .01, vjust = 0.2, color="black") +
    geom_line(data = mean,mapping = aes(x = group, y = average, group=1),color="red", linewidth=1.4)+
    labs(x = "Age Groups", y = paste("Annual change in",var_labels[i],sep=" "))+ theme_bw() + theme(legend.position="none")
}
violins = ggarrange(plotlist = plotlist)
ggsave("/home/max/Documents/Projects/Diffusion/UKBlong/Results/violins.pdf",violins, width = 15, height = 20)

## CREATE THE SAME PLOT STRATIFYING FOR MALES AND FEMALES
# and loop over all the metrics
for (i in 1:length(lm_list)){
  # estimate mean by group
  mean <- lm_list[[i]] %>% 
    group_by(group) %>% 
    summarize(average = mean(ROC)) %>%
    ungroup()
  # plot age stratified rate of change
  plotlist[[i]] = ggplot2::ggplot(lm_list[[i]], aes(x = group, y=ROC, fill = sex)) +
    ggplot2::geom_violin(position=position_dodge(1),alpha = 4/10) +   geom_boxplot(position=position_dodge(1),width=0.1,   outlier.shape = 0,outlier.size = 0,outlier.stroke = 0)+
    # geom_signif(comparisons = list(c("Male", "Female")), map_signif_level=c((0.001/(6*length(lm_list))),0.01/(6*length(lm_list)),0.05/(6*length(lm_list))), y_position = 5, tip_length = .01, vjust = 0.2, color="black") +
    # geom_signif(comparisons = list(c(">55", "65-75")), map_signif_level=c((0.001/(6*length(lm_list))),0.01/(6*length(lm_list)),0.05/(6*length(lm_list))), y_position = 6, tip_length = .01, vjust = 0.2, color="black") +
    # geom_signif(comparisons = list(c(">55", "75+")), map_signif_level=c((0.001/(6*length(lm_list))),0.01/(6*length(lm_list)),0.05/(6*length(lm_list))), y_position = 8, tip_length = .01, vjust = 0.2, color="black") +
    # geom_signif(comparisons = list(c("55-65", "65-75")), map_signif_level=c((0.001/(6*length(lm_list))),0.01/(6*length(lm_list)),0.05/(6*length(lm_list))), y_position = 4.8, tip_length = .01, vjust = 0.2, color="black") +
    # geom_signif(comparisons = list(c("65-75", "75+")), map_signif_level=c((0.001/(6*length(lm_list))),0.01/(6*length(lm_list)),0.05/(6*length(lm_list))), y_position = 4.6, tip_length = .01, vjust = 0.2, color="black") +
    # geom_signif(comparisons = list(c("55-65", "75+")), map_signif_level=c((0.001/(6*length(lm_list))),0.01/(6*length(lm_list)),0.05/(6*length(lm_list))), y_position = 7, tip_length = .01, vjust = 0.2, color="black") +
    # geom_line(data = mean,mapping = aes(x = group, y = average, group=1),color="red", linewidth=1.4)+
    labs(x = "Age Groups", y = paste("Annual change in",var_labels[i],sep=" "))+ theme_bw() + theme(legend.position="none")
}
violins1 = ggarrange(plotlist = plotlist)
# females red, males blue
ggsave("/home/max/Documents/Projects/Diffusion/UKBlong/Results/SEX_violins.pdf",violins1, width = 15, height = 20)

# make name list of outcome variables in as in data frame
outcome_vars = colnames(ROC[1:26])

## 3) Age-stratified Regional Annual Change ####
# this hasn't been done, as the trends can only be visually inspected, and inspecting and evaluating more than 1k trends into a simple answer was deemed iunfruitful. 

## 4) Global Annual Change AGE ASSOCIATIONS ####
# now, we can subtract one time point from the other, the last column is age, for which we estimate the delta
ROC = T2[1:27]-T1[1:27]
lm_list = list()
age = T1$age
sex = T1$sex
site = T1$site_t3
for (i in 1:26){
  tmpROC = scale(ROC[i]/abs(ROC[27]), center = FALSE)
  lm_list[[i]] = data.frame(age, sex, site, tmpROC)
  tmp_model = lmer(tmpROC ~ age:sex + sex + (1|site), data = lm_list[[i]])
  lm_list[[i]]$fit = predict(tmp_model)
  colnames(lm_list[[i]]) = c("age","sex","site","ROC", "fit")
}
# make name list of outcome variables in as in data frame
outcome_vars = colnames(ROC[1:26])
# loop for plots
lm_plot = list()
for (i in 1:length(var_labels)){
  lm_plot[[i]] = ggplot(data=lm_list[[i]],aes(x=age, y=fit))+ 
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
    ylab(paste(var_labels[i]))+ 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())+
    theme(legend.position = "none")+
    stat_cor(label.y = 1)
    #labs(title = (paste("linear R2 = ",round(lm_r.sq1[i],digits = 4), "      ","non-linear R2 = ", round(gam_r.sq1[i],digits = 4))))
}
plot= do.call("grid.arrange", c(lm_plot,ncol = 3))
ggsave("/home/max/Documents/Projects/Diffusion/UKBlong/Results/annual_change_plots.pdf",plot, width = 15, height = 20)

## 5) Regional Annual Changes AGE ASSOCIATIONS ####
# load data (without DRAD extra outliers removed as there were many outliers)
T1 = read.csv("T1_noDRAD.csv")
T2 = read.csv("T2_noDRAD.csv")
# exclude outliers / impossible values
T1 = T1[!T1$eid %in% outliers,]
T2 = T2[!T2$eid %in% outliers,]
T1 = T1 %>% dplyr::select(!(contains("Drad_extra"))) %>% dplyr::select(everything(), eid, age, sex, site_t3, site_t4)
T2 = T2 %>% dplyr::select(!(contains("Drad_extra"))) %>% dplyr::select(everything(), eid, age, sex, site_t3, site_t4)
# now, we can subtract one time point from the other
ROC = T2[2:1865]-T1[4:1867]
ROC = ROC %>% select(-eid)
lm_list = list()
age = T1$age
sex = T1$sex
site = T1$site_t3
age_diff = T2$age - T1$age
for (i in 1:1863){
  tmpROC = scale(ROC[i]/abs(ROC[28]), center = FALSE)
  lm_list[[i]] = data.frame(age, sex, site, tmpROC)
  tmp_model = lmer(tmpROC ~ age:sex + sex + (1|site), data = lm_list[[i]])
  lm_list[[i]]$fit = predict(tmp_model)
  colnames(lm_list[[i]]) = c("age","sex","site","ROC", "fit")
}
# plotting nearly 2000 scatter plots is messy. Hence, we summarize the findings using correlations between age and annual change:
cors = c()
ps = c()
col_names = colnames(T1[4:1863])
for (i in 1:length(col_names)){
  cors[i]=cor((lm_list[[i]]$fit), lm_list[[i]]$age)
  ps[i] = cor.test((lm_list[[i]]$fit), lm_list[[i]]$age)$p.value
}
cor_df = data.frame(col_names,cors, ps)
# we keep only the significant associations after Bonferroni correction
cor_df = cor_df %>% filter(ps < .05/length(col_names))
pos_cors = cor_df %>% filter(cors > 0)
neg_cors = cor_df %>% filter(cors < 0)

# the relationship between age and absolute annual change can be summarized in density plots
abs_cors1 = abs(cor_df$cors)
abs_cors = data.frame(abs_cors1,cor_df$col_names)
colnames(abs_cors) = c("cors", "col_names")

# Basic density
p1 = ggplot(abs_cors, aes(x=cors)) +
  geom_density(fill="gray")+
  geom_vline(aes(xintercept=0), color="blue",linetype="dashed")+
  labs(title="Absolute RoC-age relationship",x="RoC-age relationship", y = "Density")+
  theme_classic()
p1.1 = ggplot(cor_df, aes(x=cors)) +
  geom_density(fill="gray")+
  geom_vline(aes(xintercept=0), color="blue",linetype="dashed")+
  labs(title="All RoC-age relationships",x="RoC-age relationship", y = "Density")+
  theme_classic()
p2 = ggplot(neg_cors, aes(x=cors)) +
  geom_density(fill="gray")+
  geom_vline(aes(xintercept=0), color="blue",linetype="dashed")+
  labs(title="Negative RoC-age relationship",x="RoC-age relationship", y = "Density")+
  theme_classic()
p3 = ggplot(pos_cors, aes(x=cors)) +
  geom_density(fill="gray")+
  geom_vline(aes(xintercept=0), color="blue",linetype="dashed")+
  labs(title="Positive RoC-age relationship",x="RoC-age relationship", y = "Density")+
  theme_classic()
p4 = ggarrange(p1.1,p1,p2,p3,ncol = 2, nrow = 2)
ggsave("/home/max/Documents/Projects/Diffusion/UKBlong/Results/density_RoC_age_correlations.pdf",p4,height = 7, width = 7)
# across all metrics the absolute annual rate of change - age relationship is > 0.
# we run t-tests and estimate Cohen's d
t.test(abs(cor_df$cors),mu=0)
mean(abs(cor_df$cors))/sd(abs(cor_df$cors))
t.test(pos_cors$cors,mu=0)
mean(abs(pos_cors$cors))/sd(abs(neg_cors$cors))

t.test(neg_cors$cors,mu=0)
mean(abs(neg_cors$cors))/sd(abs(pos_cors$cors))

# (however, this does not tell us whether the rate of change trends towards 0, or away from it!)

# # if interested in the fornix, specifically, uncomment the following
# 
# 
## 6) Fornix Annual Changes AGE ASSOCIATION ####
# # load data (without DRAD extra outliers removed as there were many outliers)
# # T1 = read.csv("T1_noDRAD.csv")
# # T2 = read.csv("T2_noDRAD.csv")
# # we have to exclude DRAD extra from our analysis, as this metric did poor in the QC
# T1 = T1 %>% dplyr::select(!(contains("Drad_extra")))%>% dplyr::select(contains("fornix"), age, sex, site_t3, eid)# %>% dplyr::select(!Drad_extra_Mean)
# T2 = T2 %>% dplyr::select(!(contains("Drad_extra")))%>% dplyr::select(contains("fornix"), age, sex, site_t3, eid)# %>% dplyr::select(!Drad_extra_Mean)
# 
# # now, we can substract one time point from the other
# ROC = T2[1:81]-T1[1:81]
# lm_list = list()
# age = T1$age
# sex = T1$sex
# site = T1$site_t3
# for (i in 1:ncol(ROC)){
#   tmpROC = scale(ROC[i]/abs(ROC[ncol(ROC)]), center = FALSE)
#   lm_list[[i]] = data.frame(age, sex, site, tmpROC)
#   colnames(lm_list[[i]]) = c("age","sex","site","ROC")
# }
# lm_plot = list()
# col_names = colnames(ROC)
# col_names = str_replace(col_names, "Striaterminalis", "St")
# for (i in 1:ncol(ROC)){
#   lm_plot[[i]] = ggplot(data=lm_list[[i]],aes(x=age, y=ROC))+ 
#     #geom_line(aes(y=fit,group=ID),alpha=.4,lwd=.5) + #,col=TP
#     #geom_point(aes(y=fit,shape=factor(TP), col=TP),alpha=.6) +#
#     #geom_ribbon(aes(group=TP,ymin=lwr,ymax=upr),alpha=0.3,fill="grey") + 
#     #geom_smooth(method = "lm", color= "#E69F00") +
#     geom_smooth(method = "gam", color = "#56B4E9") +
#     geom_line(aes(y=0),linetype="dotted")+
#     # stat_regline_equation(
#     #   aes(label =  paste(..eq.label.., ..adj.rr.label.., sep = "~~~~")),
#     #   formula = formula(paste("data"[i]"~age")), data =data) +
#     #geom_smooth(method = "lm") +
#     #stat_regline_equation() +
#     #stat_cor()+
#     xlab("Age")+
#     ylab (paste(col_names[i]))+ 
#     theme_bw() + 
#     theme(panel.grid.major = element_blank(), 
#           panel.grid.minor = element_blank())+
#     theme(legend.position = "none")
#     #stat_cor(label.y = 1)
#   #labs(title = (paste("linear R2 = ",round(lm_r.sq1[i],digits = 4), "      ","non-linear R2 = ", round(gam_r.sq1[i],digits = 4))))
# }
# plot= do.call("grid.arrange", c(lm_plot))
# ggsave("/home/max/Documents/Projects/Diffusion/UKBlong/Results/FORNIX_annual_change_plots.pdf",plot, width = 20, height = 20)
# 
# # the many scatter plots visualisation is messy, we summarise the findings using correlations:
# cors = c()
# for (i in 1:ncol(ROC)){
#   cors[i]=cor(abs(lm_list[[i]]$ROC), lm_list[[i]]$age)
# }
# length(col_names)
# ncol(ROC)
# cor_df = data.frame(col_names,cors)
# cor_df
# cor_df %>% filter(abs(cors)>0.1)
# plot(density(cors))
# 
# # across all metrics in the fornix, age associations of the annual rate of change > 0
# t.test(cor_df$cors,mu=0)