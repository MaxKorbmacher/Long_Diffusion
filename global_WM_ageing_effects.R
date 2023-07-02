# Primer for longitudinal data analysis with 2 time points nested in different clusters (e.g., scanner sites) + covariates
# Max Korbmacher, 13.02.2022
#
################################################# #
# Purpose:                                        #
## analyse longitudinal data using three models:  #
### a) vanilla linear regression                  #
### b) random intercept model                     #
### c) marginal models                            #
################################################# #
# Contents / Parts:                               #
## 1) PREP ENV                                    #
## 2) LOAD DATA                                   #
## 3) ANALYSE                                     #
#### 3.1) run models                              #
#### 3.2) get betas                               #
#### 3.3) diagnostics                             #
###### 3.3.1 HETEROSCEDASTICITY                   #
###### 3.3.2 MODEL FIT                            #
## 4) PLOT                                        #
################################################# #
## 1) PREPARATIONS ####
#### 1.1) PREP ENV ####

# load packages and install if not already installed
remove.packages("sjPlot")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(lme4, nlme, ggplot2, tidyverse, lm.beta, remotes, ggpubr, grid, lmtest, car, lmtest,lmeInfo,lmerTest,sjstats,effsize)
#install_version("sjstats","0.17.7")
library(sjstats)

#### 1.2) PREP DATA ####
# set working directory where the data files are situated
setwd("/home/max/Documents/Projects/Diffusion/UKBlong/data_export/")
# load data (without DRAD extra outliers removed as there were many outliers)
T1 = read.csv("T1_noDRAD.csv")
T2 = read.csv("T2_noDRAD.csv")
# we have to exclude DRAD extra from our analysis, as this metric did poor in the QC
T1 = T1 %>% dplyr::select(contains("mean"), age, sex, site_t3, eid) %>% dplyr::select(!Drad_extra_Mean)
T2 = T2 %>% dplyr::select(contains("mean"), age, sex, site_t3, eid) %>% dplyr::select(!Drad_extra_Mean)

# make unstandardized data frame prior to standardization
T1_unstd = T1
T2_unstd = T2
data_unstd = rbind(T1, T2)

# standardize each diffusion metrics (Z-scoring) for later comparison
for (i in 1:length(outcome_vars)){
  T1[,i] = (T1[,i]-(mean(T1[,i])))/(sd(T1[,i]))
  T2[,i] = (T2[,i]-(mean(T2[,i])))/(sd(T2[,i]))
}

# 
# create time point dummy TP
T1$TP = replicate(nrow(T1),0)
T2$TP = replicate(nrow(T1),1)
T2$time = T2$age - T1$age
T1$time = 0

# make single df
data = rbind(T1, T2)
data$TP = factor(data$TP)
levels(data$TP) = c("baseline","repeat")
data_unstd$TP = factor(data$TP)
levels(data_unstd$TP) = c("baseline","repeat")

# make name list of outcome variables in as in data frame
outcome_vars = colnames(data[1:27])

# make name labels for the outcome_vars (more presentable)
var_labels = c("BRIA - V intra", "BRIA - V extra", "BRIA - V csf", "BRIA - micro RD","BRIA - micro FA",
               "BRIA - micro AX", "BRIA - micro ADC", "BRIA - DAX intra", "BRIA - DAX extra",
               
               "DKI - RK", "DKI - AK", "DKI -MK",
               
               "DTI - FA", "DTI - MD", "DTI - RD", "DTI - AD",
               
               "SMT - FA", "SMT - long", "SMT - MD", "SMT - trans",
               
               "SMTmc - intra","SMTmc - extra MD", "SMTmc - extra trans",  "SMTmc - diff",
               
               "WMTI - axEAD", "WMTI - AWF", "WMTI - radEAD")

## 2) DEMOGRAPHICS ####
# N
nrow(T1)
# age
min(T1$age)
max(T1$age)
mean(T1$age)
sd(T1$age)
min(T2$age)
max(T2$age)
mean(T2$age)
sd(T2$age)
min(T2$age-T1$age)
max(T1$age)
mean(T2$age-T1$age)
sd(T2$age-T1$age)

# sex
table(T1$sex)
prop.table(table(data$sex))

# site
## check whether site variable is the same for both time points (yes, they are)
setequal(T2$site_t4,T1$site_t3)
## check the proportions
prop.table(table(T1$site_t3))
## and raw numbers, if wanted
#table(T1$site_t3)

## 3) ANALYSE ####

#################### #
#### 3.0) t.tests for time point differences #### 
#################### #
p.val = c()
t.val = c()
CohensD = c()
mean_TP1 = c()
sd_TP1 = c()
mean_TP2 = c()
sd_TP2 = c()
d = c()
d_upper = c()
d_lower = c()
for (i in outcome_vars){
  f = formula(paste("data_unstd$",i,"~data_unstd$TP", sep = ""))
  p.val[i] = t.test(f, paired = TRUE)$p.value
  t.val[i] = t.test(f, paired = TRUE)$statistic
}
for (i in 1:length(outcome_vars)){
  mean_TP1[i] = mean(T1_unstd[,i])
  sd_TP1[i] = sd(T1_unstd[,i])
  mean_TP2[i] = mean(T2_unstd[,i])
  sd_TP2[i] = sd(T2_unstd[,i])
  d[i] = cohen.d(T2_unstd[,i], T1_unstd[,i], paired=TRUE)$estimate
  d_lower[i] = as.numeric(cohen.d(T2_unstd[,i], T1_unstd[,i], paired=TRUE)$conf.int[1])
  d_upper[i] = as.numeric(cohen.d(T2_unstd[,i], T1_unstd[,i], paired=TRUE)$conf.int[2])
}
p.adj = p.adjust(p.val,method = "bonferroni")
paired_t_out = data.frame(var_labels, mean_TP1, sd_TP1, mean_TP2, sd_TP2, t.val, p.val, p.adj, d,d_lower, d_upper)
paired_t_out = data.frame(lapply(paired_t_out,function(x) if(is.numeric(x)) round(x, 3) else x))
colnames(paired_t_out) = c("Metric", "Mean TP1","SD TP1",  "Mean TP2","SD TP1", "T", "p", "adjusted p", "Cohens's d", "Lower 95% CI", "Upper 95% CI")
write.csv(paired_t_out, "/home/max/Documents/Projects/Diffusion/UKBlong/Results/paired_t_results.csv")

# make a simple figure summarizing the effect sizes
colnames(paired_t_out) = c("Metric", "Mean_TP1","SD_TP1", "Mean_TP2","SD_TP1", "T", "p", "adjusted_p", "Cohens_d", "LL", "UL")
colnames(paired_t_out) = make.unique(names(paired_t_out))
paired_t_out$Index = seq(1:27)
plot1 <- ggplot(paired_t_out, aes(y = Index, x = Cohens_d)) +
  geom_point(shape = 18, size = 5) +  
  geom_errorbarh(aes(xmin = LL, xmax = UL), height = 0.25) +
  geom_vline(xintercept = 0, color = "red", linetype = "dashed", cex = 1, alpha = 0.5) +
  scale_y_continuous(name = "", breaks=1:27, labels = paired_t_out$Metric, trans = "reverse") +
  xlab("Cohen's d (95% CI)") + 
  ylab(" ") + 
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.text.x.bottom = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"))
ggsave(plot1, filename = "/home/max/Documents/Projects/Diffusion/UKBlong/Results/TP_differences_CohensD_plot.pdf", width = 6, height = 9)

#################### #
#### 3.1) run models #### 
#################### #

# make sure that TP2 -TP1 to see the right direction of developments 
data$TP = ordered(data$TP, levels =  c("repeat","baseline"))

# make empty list to be filled
RI = list()
# loop over outcome_vars and run the same LM, RI and MM models for each outcome var.
for (i in outcome_vars){
  f = formula(paste(i,"~age*sex + TP+(1|site_t3/eid)"))
  RI[[i]] = lmer(f, data = data)
}

#################### ###
#### 3.2) get betas  #### 
#################### ###

### STANDARDIZED BETAS

# compile standardized beta values for a variable specified in ROW_NUMBER, (age ROW_NUMBER = 1, sex = 2, time point = 3, age-sex interaction = 4)
# change the ROW_NUMBER value accordingly!

# empty vectors to be filled
RIp = c()
RIstde = c()
RIb = c()
betas = list()
diffusion_models = c(replicate(9,"BRIA"), replicate(3, "DKI"), replicate(4, "DTI"), replicate(4, "SMT"), replicate(4, "SMTmc"), replicate(3, "WMTI"))
# estimate betas, standard errors, t-vals, p-vals
for (ROW_NUMBER in 1:4){
  for (i in 1:length(outcome_vars)){
  RIb[i] = (summary(RI[[i]])$coefficients[ROW_NUMBER+1,1])
  RIstde[i] = summary(RI[[i]])$coefficients[ROW_NUMBER+1,2]
  RIp[i] = summary(RI[[i]])$coefficients[ROW_NUMBER+1,5]
}
  RI_betas = data.frame(var_labels, RIb, RIstde, RIp)
  # correct labelling of column names
  colnames(RI_betas) =  c("Outcome","Std.Beta", "Std.Err", "p")
  # Bonferroni correction
  RI_betas$p.adj = p.adjust(RI_betas$p,method = "bonferroni")
  RI_betas$Model = diffusion_models
  betas[[ROW_NUMBER]] = data.frame(lapply(RI_betas,function(x) if(is.numeric(x)) round(x, 3) else x))
}

# save  
write.csv(betas[[1]], "/home/max/Documents/Projects/Diffusion/UKBlong/Results/effect_of_age.csv")
write.csv(betas[[2]], "/home/max/Documents/Projects/Diffusion/UKBlong/Results/effect_of_sex.csv")
write.csv(betas[[3]], "/home/max/Documents/Projects/Diffusion/UKBlong/Results/effect_of_timepoint.csv")
write.csv(betas[[4]], "/home/max/Documents/Projects/Diffusion/UKBlong/Results/effect_of_age_sex_interaction.csv")

# 
# ##################### #
# #### 3.3) diagnostics #### 
# ##################### #
# # this section can be un-commented, if interested in the results (the short story is: there is heteroscedasticity across time, site and sex)
# ###### 3.3.1 HETEROSCEDASTICITY ####
# 
# # 1a) inspect variability across metrics for the two measures [TIME]
# df01 = filter(data, TP==0)
# df02 = filter(data, TP==1)
# df01 = data.frame(apply(df01,2,sd))
# colnames(df01) = c("SD")
# df02 = data.frame(apply(df02,2,sd))
# colnames(df02) = c("SD")
# df03 = df01-df02
# rownames(df03) = rownames(df02)
# # display difference: variability increases
# df03
# print("Ergo, the variability in diffusion scores increases over time")
# 
# # using an F-test to compare variability across time points (not actual time)
# fres = c()
# for (i in outcome_vars){
#   f = formula(paste(i,"~TP"))
#   fres[i] = var.test(f,data=data)$p.val
# }
# SDdiff = cbind(df01, df02)
# colnames(SDdiff) = c("TP1", "TP2")
# SDdiff = SDdiff[1:length(outcome_vars),]
# SDdiff$SameVar = fres
# # F-test says that the following values have statistically the same variance
# SDdiff %>% filter(SameVar > 0.05)
# # and those variances differ:
# SDdiff %>% filter(SameVar < 0.05)
# 
# # 1b) repeat for [SEX & TIME]
# df01.0 = filter(data, sex==0 & TP == 0)
# df01.1 = filter(data, sex==0 & TP == 1)
# df02.0 = filter(data, sex==1 & TP == 0)
# df02.1 = filter(data, sex==1 & TP == 1)
# varlist = list(df01.0, df01.1, df02.0, df02.1)
# varres = list()
# for (i in 1:length(varlist)){
#   varres[[i]] = data.frame(apply(varlist[[i]],2,sd))
#   colnames(varres[[i]]) = ("SD")
# }
# SDdiff = data.frame(varres[[1]]-varres[[2]], varres[[3]]-varres[[4]])
# colnames(SDdiff) = c("female", "male")
# a = rownames(SDdiff %>% filter(female > 0 & male > 0))
# paste("The variability in",a, "score DEcreases over time and across sex")
# a = rownames(SDdiff %>% filter(female < 0 & male < 0))
# paste("The variability in",a, "score INcreases over time and across sex")
# print("Ergo, the variability in diffusion scores increases over time and across sexes")
# 
# # looking only at sex
# df01 = filter(data, sex==0)
# df02 = filter(data, sex==1)
# varlist = list(df01, df02)
# varres = list()
# for (i in 1:length(varlist)){
#   varres[[i]] = data.frame(apply(varlist[[i]],2,sd))
#   colnames(varres[[i]]) = ("SD")
# }
# SDdiff = data.frame(varres[[1]],varres[[2]])
# colnames(SDdiff) = c("female", "male")
# SDdiff
# print("Looking only at sex, it seems like there is some small variation in measures.")
# 
# # using an F-test to compare sexes
# fres = c()
# for (i in outcome_vars){
#   f = formula(paste(i,"~sex"))
#   fres[i] = var.test(f,data=data)$p.val
# }
# SDdiff = SDdiff[1:length(outcome_vars),]
# SDdiff$SameVar = fres
# # F-test says that the following values have statistically the same variance
# SDdiff %>% filter(SameVar > 0.05)
# # and those variances differ:
# SDdiff %>% filter(SameVar < 0.05)
# 
# # using an F-test to compare variability across sexes and time points (not actual time)
# fres = c()
# for (i in outcome_vars){
#   f = formula(paste(i,"~TP"))
#   fres[i] = var.test(f,data=df01)$p.val
# }
# female_tdiff = data.frame(SDdiff$female,fres)
# # F-test says that the following values have statistically the same variance
# female_tdiff %>% filter(fres > 0.05)
# # and those variances differ:
# female_tdiff %>% filter(fres < 0.05)
# 
# ## >> nearly all of the diffusion scores for females do not seem to change over time!
# 
# fres = c()
# for (i in outcome_vars){
#   f = formula(paste(i,"~TP"))
#   fres[i] = var.test(f,data=df02)$p.val
# }
# male_tdiff = fres
# 
# male_tdiff = data.frame(SDdiff$male,fres)
# # F-test says that the following values have statistically the same variance
# male_tdiff %>% filter(fres > 0.05)
# # and those variances differ:
# male_tdiff %>% filter(fres < 0.05)
# # also for males, there seems to be no variability for sexes across time
# 
# # 1c) repeat for [SCANNER SITE & TIME]
# df01.0 = filter(data, site_t3=="Cheadle" & TP == 0)
# df01.1 = filter(data, site_t3=="Cheadle" & TP == 1)
# df02.0 = filter(data, site_t3=="Newcastle" & TP == 0)
# df02.1 = filter(data, site_t3=="Newcastle" & TP == 1)
# df03.0 = filter(data, site_t3=="Reading" & TP == 0)
# df03.1 = filter(data, site_t3=="Reading" & TP == 1)
# varlist = list(df01.0, df01.1, df02.0, df02.1, df03.0, df03.1)
# varres = list()
# for (i in 1:length(varlist)){
#   varres[[i]] = data.frame(apply(varlist[[i]],2,sd))
#   colnames(varres[[i]]) = ("SD")
# }
# SDdiff = data.frame(varres[[1]]-varres[[2]], varres[[3]]-varres[[4]], varres[[5]]-varres[[6]])
# colnames(SDdiff) = c("Cheadle", "Newcastle", "Reading")
# a = rownames(SDdiff %>% filter(Cheadle > 0 & Newcastle > 0 & Reading > 0))
# paste("The variability in",a, "score DEcreases over time and across sites")
# a = rownames(SDdiff %>% filter(Cheadle < 0 & Newcastle < 0 & Reading < 0))
# paste("The variability in",a, "score INcreases over time and across sites")
# print("Ergo, the variability in diffusion scores increases over time and across sites")
# 
# # looking only at scanner site
# df01 = filter(data, site_t3=="Cheadle")
# df02 = filter(data, site_t3=="Newcastle")
# df03 = filter(data, site_t3=="Reading")
# varlist = list(df01, df02, df03)
# varres = list()
# for (i in 1:length(varlist)){
#   varres[[i]] = data.frame(apply(varlist[[i]],2,sd))
#   colnames(varres[[i]]) = ("SD")
# }
# SDdiff = data.frame(varres[[1]],varres[[2]], varres[[3]])
# colnames(SDdiff) = c("Cheadle", "Newcastle", "Reading")
# SDdiff = SDdiff[1:27,]
# print("Looking only at scanner site, despite some exceptions, it seems like there is some small variation in measures.")
# 
# # we check Heteroscedasticity across scanner sites this using F-tests
# fres = c()
# fres1 = c()
# fres2 = c()
# for (i in outcome_vars){
#   f = formula(paste(i,"~TP"))
#   #Cheadle
#   fres[i] = var.test(f,data=df01)$p.val
#   #Newcastle
#   fres1[i] = var.test(f,data=df02)$p.val
#   #Reading
#   fres2[i] = var.test(f,data=df03)$p.val
# }
# Cheadle_tdiff = data.frame(SDdiff$Cheadle,fres)
# # F-test says that the following values have statistically the same variance over time in Cheadle
# Cheadle_tdiff %>% filter(fres > 0.05)
# # and those variances differ over time in Cheadle:
# Cheadle_tdiff %>% filter(fres < 0.05)
# 
# Newcastle_tdiff = data.frame(SDdiff$Newcastle,fres1)
# # F-test says that the following values have statistically the same variance over time in Newcastle
# Newcastle_tdiff %>% filter(fres1 > 0.05)
# # and those variances differ over time in Newcastle:
# Newcastle_tdiff %>% filter(fres1 < 0.05)
# 
# Reading_tdiff = data.frame(SDdiff$Reading,fres2)
# # F-test says that the following values have statistically the same variance over time in Reading
# Reading_tdiff %>% filter(fres2 > 0.05)
# # and those variances differ over time in Reading:
# Reading_tdiff %>% filter(fres2 < 0.05)
# # this shows that the time-dependent variability is different across sites!
# 
# # 2) test heteroskedasticity using Breush Pagan and NCV Test
# # H0 is that the variance of the residuals is constant for 
# # Breush Pagan p > 0.05 is homoscedastic
# # NCV p > .05 is homoscedastic
# heterosced = data.frame(matrix(ncol = 3, nrow = length(outcome_vars)))
# heterosced[,1] = outcome_vars
# for (i in 1:length(outcome_vars)){
#   heterosced[i,2] = lmtest::bptest(LM[[i]])$p.val
#   heterosced[i,3] = car::ncvTest(LM[[i]])$p
# }
# colnames(heterosced) = c("Var", "BP_Test", "NCV_Test")
# 
# # now we can display which metrics are homoscedastic based on both tests' results
# heterosced %>% filter(BP_Test > .05 & NCV_Test > .05)
# # and which ones are heteroscedastic based on both tests  results
# heterosced %>% filter(BP_Test < .05 & NCV_Test < .05)
# # >> most metrics are heteroscedastic
# 
# 
# #
# #
# #
# #
# print("Conclusion variance and heteroscedasticity tests: there is heteroscedasticity across time, site and sex, which needs to be accounted for.")

###### 3.3.2 MODEL FIT ####

# compute AIC, BIC, R^2
AIC = c()
BIC = c()
conditionalR2 = c()
marginalR2 = c()
for (i in 1:length(outcome_vars)){
  AIC[i] = AIC(RI[[i]])
  BIC[i] = BIC(RI[[i]])
  conditionalR2[i] = as.numeric(performance::r2(RI[[i]])[1])
  marginalR2[i] = as.numeric(performance::r2(RI[[i]])[2])
}
conditionalR2 = unlist(conditionalR2)
marginalR2 = unlist(marginalR2)
diagnostics=data.frame(var_labels,AIC,BIC,(conditionalR2),(marginalR2))
colnames(diagnostics) = c("Metric", "AIC", "BIC", "Conditional R2", "Marginal R2")
write.csv(diagnostics,"/home/max/Documents/Projects/Diffusion/UKBlong/Results/MLM_diagnostic.csv")


### 4) PLOT BETAS ####
# this is to plot the effect of age and time point on diffusion metrics

# create a plot function for standardized betas
## NOTE:
## ylim needs to be adapted for optimal visualisation
## A name vector for the variable names should be passed into the function for better Outcome variable names
## In that case, the ylab() can also be dropped
bplot = function(Beta_Table){
  ggplot(Beta_Table, aes(x=Outcome, y=Std.Beta, color=diffusion_models,stat="identity")) + 
    geom_pointrange(aes(ymin=Std.Beta-Std.Err, ymax=Std.Beta+Std.Err), stat="identity") +scale_colour_brewer(palette="Dark2") + theme_bw() + theme(legend.position = "none") +
    ylab("Standardized Beta") + xlab("") + scale_x_discrete(limits = rev(factor(var_labels)))+ coord_flip()
}
age_betas=bplot(betas[[1]])
timepoint_betas=bplot(betas[[3]])
agesex_betas=bplot(betas[[4]])
sex_betas=bplot(betas[[2]])

# merge & save plots
figure = ggarrange(age_betas, timepoint_betas, sex_betas, agesex_betas,labels = c("a","b","c","d"))
ggsave("/home/max/Documents/Projects/Diffusion/UKBlong/Results/DiffusionAssocitation.pdf", width = 10, height = 9, figure)
figure
