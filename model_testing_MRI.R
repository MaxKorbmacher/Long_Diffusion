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

## 1) PREP ENV ####

# load packages and install if not already installed
remove.packages("sjPlot")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(lme4, nlme, ggplot2, tidyverse, lm.beta, remotes, ggpubr, grid, lmtest, car, lmtest)
#install_version("sjstats","0.17.7")
library(sjstats)

## 2) PREP DATA ####
# set working directory where the data files are situated
setwd("/home/max/Documents/Projects/Diffusion/UKBlong/data_export/")
# load data (without DRAD extra outliers removed as there were many outliers)
T1 = read.csv("T1_noDRAD.csv")
T2 = read.csv("T2_noDRAD.csv")
# hence, we have to exclude DRAD extra from our analysis, and focus on global average scores here
T1 = T1 %>% dplyr::select(contains("mean"), age, sex, site_t3, eid) %>% dplyr::select(!Drad_extra_Mean)
T2 = T2 %>% dplyr::select(contains("mean"), age, sex, site_t3, eid) %>% dplyr::select(!Drad_extra_Mean)

# create time point dummy TP
T1$TP = replicate(nrow(T1),0)
T2$TP = replicate(nrow(T1),1)
T2$time = T2$age - T1$age
T1$time = 0

# make single df
data = rbind(T1, T2)
data$TP = factor(data$TP)
levels(data$TP) = c("baseline","repeat")

## 3) ANALYSE ####
#################### #
#### 3.1) run models #### 
#################### #

# make name list of outcome variables, and empty lists (to be filled with fitted models) 
outcome_vars = colnames(data[1:27])
LM = list()
MM = list()
RI = list()
# loop over outcome_vars and run the same LM, RI and MM models for each outcome var.
for (i in outcome_vars){
  f = formula(paste(i,"~age*sex + site_t3 + TP"))
  f2 = formula(paste(i,"~age*sex + TP+(1|site_t3)+(1|eid)"))
  f3 = formula(paste(i,"~age*sex+site_t3+TP"))
  LM[[i]] = lm(f,data=data)
  RI[[i]] = lmer(f2, data = data)
  MM[[i]] = gls(model = f3, data=data, correlation = corSymm(form = ~1|eid), weights = varIdent(form = ~1|factor(TP) * factor(site_t3) * factor(sex)))
}

# for the marginal model, we run the analysis on standardized coefficients 
# standardise each diffusion metrics (Z-scoring) for later comparison
for (i in 1:length(outcome_vars)){
  T1[,i] = (T1[,i]-(mean(T1[,i])))/(sd(T1[,i]))
  T2[,i] = (T2[,i]-(mean(T2[,i])))/(sd(T2[,i]))
}

# create time point dummy TP
T1$TP = replicate(nrow(T1),0)
T2$TP = replicate(nrow(T1),1)
T2$time = T2$age - T1$age
T1$time = 0
# make single df
data2 = rbind(T1, T2)

# run MM on standardized data
LMstd = list()
RIstd = list()
MMstd = list()
for (i in outcome_vars){
  f = formula(paste(i,"~age*sex + site_t3 + TP"))
  f2 = formula(paste(i,"~age*sex + TP+(1|site_t3)+(1|eid)"))
  f3 = formula(paste(i,"~age*sex+site_t3+TP"))
  LMstd[[i]] = lm(f,data=data2)
  RIstd[[i]] = lmer(f2, data = data2)
  MMstd[[i]] = gls(model = f3, data=data2, correlation = corSymm(form = ~1|eid), weights = varIdent(form = ~1|factor(TP) * factor(site_t3) * factor(sex)))
}

#################### #
#### 3.2) get betas  #### 
#################### #

### STANDARDIZED BETAS

# compile standardized beta values for a variable specified in ROW_NUMBER, e.g. for age ROW_NUMBER = 1
ROW_NUMBER = 1
# change the ROW_NUMBER value accordingly!

# estimate betas, standard errors, t-vals, p-vals
LMb = c()
LMstde = c()
LMt = c()
LMp = c()
RIb = c()
RIstde = c()
MMb = c()
MMstde = c()
MMt = c()
MMp = c()
for (i in 1:length(outcome_vars)){
  LMb[i] = as.numeric(LMstd[[i]]$coefficients[ROW_NUMBER+1])
  LMstde[i] = as.numeric(sqrt(diag(vcov(LMstd[[i]])))[ROW_NUMBER+1])
  LMt[i] = as.numeric(coef(summary(LMstd[[i]]))[, "t value"][ROW_NUMBER+1])
  LMp[i] = as.numeric(coef(summary(LMstd[[i]]))[, "Pr(>|t|)"][ROW_NUMBER+1])
  RIb[i] = (summary(RIstd[[i]])$coefficients[ROW_NUMBER+1,1])
  RIstde[i] = summary(RIstd[[i]])$coefficients[ROW_NUMBER+1,2]
  MMb[i] = as.numeric(MMstd[[i]]$coefficients[ROW_NUMBER+1])
  MMstde[i] = coef(summary(MMstd[[i]]))[ROW_NUMBER+1,2]
  MMt[i] = as.numeric(coef(summary(MMstd[[i]]))[, "t-value"][ROW_NUMBER+1])
  MMp[i] = as.numeric(coef(summary(MMstd[[i]]))[, "p-value"][ROW_NUMBER+1])
}
# make data frames
var_labels = c("BRIA - V intra", "BRIA - V extra", "BRIA - V csf", "BRIA - micro RD","BRIA - micro FA",
                       "BRIA - micro AX", "BRIA - micro ADC", "BRIA - DAX intra", "BRIA - DAX extra",
                       
                       "DKI - RK", "DKI - AK", "DKI -MK",
                       
                       "DTI - FA", "DTI - MD", "DTI - RD", "DTI - AD",
                       
                       "SMT - FA", "SMT - long", "SMT - MD", "SMT - trans",
                       
                       "SMTmc - intra","SMTmc - extra MD", "SMTmc - extra trans",  "SMTmc - diff",
                       
                       "WMTI - axEAD", "WMTI - AWF", "WMTI - radEAD")


LM_betas = data.frame(LMb, LMstde,LMt,LMp, var_labels)
RI_betas = data.frame(RIb, RIstde, var_labels)
MM_betas = data.frame(MMb, MMstde, MMt, MMp, var_labels)
colnames(LM_betas) = c("Std.Beta", "Std.Err", "T.val", "p.val", "Outcome")
colnames(RI_betas) =  c("Std.Beta", "Std.Err", "Outcome")
colnames(MM_betas) =  c("Std.Beta", "Std.Err", "T.val", "p.val", "Outcome")

# Holm correction
MM_betas$p.adj = p.adjust(MM_betas$p.val,method = "holm")
# Write the table with p-vals for the effect of time point on cognitive measures
p_table = MM_betas %>% select(Outcome,Std.Beta,Std.Err, T.val, p.val, p.adj)
write.csv(p_table, "P_Table_Age_DIFF.csv")

### UNSTANDARDIZED BETAS
# estimate betas and standard errors
LMb = c()
LMstde = c()
RIb = c()
RIstde = c()
MMb = c()
MMstde = c()

for (i in 1:length(outcome_vars)){
  LMb[i] = as.numeric(LM[[i]]$coefficients[ROW_NUMBER+1])
  LMstde[i] = as.numeric(sqrt(diag(vcov(LM[[i]])))[ROW_NUMBER+1])
  RIb[i] = as.numeric(fixef(RI[[i]])[ROW_NUMBER+1])
  RIstde[i] = summary(RI[[i]])$coefficients[ROW_NUMBER+1,2]
  MMb[i] = as.numeric(MM[[i]]$coefficients[ROW_NUMBER+1])
  MMstde[i] = coef(summary(MM[[i]]))[ROW_NUMBER+1,2]
}

# make data frames
LM_unstd_betas = data.frame(LMb, LMstde, outcome_vars)
RI_unstd_betas = data.frame(RIb, RIstde, outcome_vars)
MM_unstd_betas = data.frame(MMb, MMstde, outcome_vars)
colnames(LM_unstd_betas) = c("UnStd.Beta", "Std.Err", "Outcome")
colnames(RI_unstd_betas) = c("UnStd.Beta", "Std.Err", "Outcome")
colnames(MM_unstd_betas) = c("UnStd.Beta", "Std.Err", "Outcome")

##################### #
#### 3.3) diagnostics #### 
##################### #

###### 3.3.1 HETEROSCEDASTICITY ####

# 1a) inspect variability across metrics for the two measures [TIME]
df01 = filter(data, TP==0)
df02 = filter(data, TP==1)
df01 = data.frame(apply(df01,2,sd))
colnames(df01) = c("SD")
df02 = data.frame(apply(df02,2,sd))
colnames(df02) = c("SD")
df03 = df01-df02
rownames(df03) = rownames(df02)
# display difference: variability increases
df03
print("Ergo, the variability in diffusion scores increases over time")

# using an F-test to compare variability across time points (not actual time)
fres = c()
for (i in outcome_vars){
  f = formula(paste(i,"~TP"))
  fres[i] = var.test(f,data=data)$p.val
}
SDdiff = cbind(df01, df02)
colnames(SDdiff) = c("TP1", "TP2")
SDdiff = SDdiff[1:length(outcome_vars),]
SDdiff$SameVar = fres
# F-test says that the following values have statistically the same variance
SDdiff %>% filter(SameVar > 0.05)
# and those variances differ:
SDdiff %>% filter(SameVar < 0.05)

# 1b) repeat for [SEX & TIME]
df01.0 = filter(data, sex==0 & TP == 0)
df01.1 = filter(data, sex==0 & TP == 1)
df02.0 = filter(data, sex==1 & TP == 0)
df02.1 = filter(data, sex==1 & TP == 1)
varlist = list(df01.0, df01.1, df02.0, df02.1)
varres = list()
for (i in 1:length(varlist)){
  varres[[i]] = data.frame(apply(varlist[[i]],2,sd))
  colnames(varres[[i]]) = ("SD")
}
SDdiff = data.frame(varres[[1]]-varres[[2]], varres[[3]]-varres[[4]])
colnames(SDdiff) = c("female", "male")
a = rownames(SDdiff %>% filter(female > 0 & male > 0))
paste("The variability in",a, "score DEcreases over time and across sex")
a = rownames(SDdiff %>% filter(female < 0 & male < 0))
paste("The variability in",a, "score INcreases over time and across sex")
print("Ergo, the variability in diffusion scores increases over time and across sexes")

# looking only at sex
df01 = filter(data, sex==0)
df02 = filter(data, sex==1)
varlist = list(df01, df02)
varres = list()
for (i in 1:length(varlist)){
  varres[[i]] = data.frame(apply(varlist[[i]],2,sd))
  colnames(varres[[i]]) = ("SD")
}
SDdiff = data.frame(varres[[1]],varres[[2]])
colnames(SDdiff) = c("female", "male")
SDdiff
print("Looking only at sex, it seems like there is some small variation in measures.")

# using an F-test to compare sexes
fres = c()
for (i in outcome_vars){
  f = formula(paste(i,"~sex"))
  fres[i] = var.test(f,data=data)$p.val
}
SDdiff = SDdiff[1:length(outcome_vars),]
SDdiff$SameVar = fres
# F-test says that the following values have statistically the same variance
SDdiff %>% filter(SameVar > 0.05)
# and those variances differ:
SDdiff %>% filter(SameVar < 0.05)

# using an F-test to compare variability across sexes and time points (not actual time)
fres = c()
for (i in outcome_vars){
  f = formula(paste(i,"~TP"))
  fres[i] = var.test(f,data=df01)$p.val
}
female_tdiff = data.frame(SDdiff$female,fres)
# F-test says that the following values have statistically the same variance
female_tdiff %>% filter(fres > 0.05)
# and those variances differ:
female_tdiff %>% filter(fres < 0.05)

## >> nearly all of the diffusion scores for females do not seem to change over time!

fres = c()
for (i in outcome_vars){
  f = formula(paste(i,"~TP"))
  fres[i] = var.test(f,data=df02)$p.val
}
male_tdiff = fres

male_tdiff = data.frame(SDdiff$male,fres)
# F-test says that the following values have statistically the same variance
male_tdiff %>% filter(fres > 0.05)
# and those variances differ:
male_tdiff %>% filter(fres < 0.05)
# also for males, there seems to be no variability for sexes across time

# 1c) repeat for [SCANNER SITE & TIME]
df01.0 = filter(data, site_t3=="Cheadle" & TP == 0)
df01.1 = filter(data, site_t3=="Cheadle" & TP == 1)
df02.0 = filter(data, site_t3=="Newcastle" & TP == 0)
df02.1 = filter(data, site_t3=="Newcastle" & TP == 1)
df03.0 = filter(data, site_t3=="Reading" & TP == 0)
df03.1 = filter(data, site_t3=="Reading" & TP == 1)
varlist = list(df01.0, df01.1, df02.0, df02.1, df03.0, df03.1)
varres = list()
for (i in 1:length(varlist)){
  varres[[i]] = data.frame(apply(varlist[[i]],2,sd))
  colnames(varres[[i]]) = ("SD")
}
SDdiff = data.frame(varres[[1]]-varres[[2]], varres[[3]]-varres[[4]], varres[[5]]-varres[[6]])
colnames(SDdiff) = c("Cheadle", "Newcastle", "Reading")
a = rownames(SDdiff %>% filter(Cheadle > 0 & Newcastle > 0 & Reading > 0))
paste("The variability in",a, "score DEcreases over time and across sites")
a = rownames(SDdiff %>% filter(Cheadle < 0 & Newcastle < 0 & Reading < 0))
paste("The variability in",a, "score INcreases over time and across sites")
print("Ergo, the variability in diffusion scores increases over time and across sites")

# looking only at scanner site
df01 = filter(data, site_t3=="Cheadle")
df02 = filter(data, site_t3=="Newcastle")
df03 = filter(data, site_t3=="Reading")
varlist = list(df01, df02, df03)
varres = list()
for (i in 1:length(varlist)){
  varres[[i]] = data.frame(apply(varlist[[i]],2,sd))
  colnames(varres[[i]]) = ("SD")
}
SDdiff = data.frame(varres[[1]],varres[[2]], varres[[3]])
colnames(SDdiff) = c("Cheadle", "Newcastle", "Reading")
SDdiff = SDdiff[1:27,]
print("Looking only at scanner site, despite some exceptions, it seems like there is some small variation in measures.")

# we check Heteroscedasticity across scanner sites this using F-tests
fres = c()
fres1 = c()
fres2 = c()
for (i in outcome_vars){
  f = formula(paste(i,"~TP"))
  #Cheadle
  fres[i] = var.test(f,data=df01)$p.val
  #Newcastle
  fres1[i] = var.test(f,data=df02)$p.val
  #Reading
  fres2[i] = var.test(f,data=df03)$p.val
}
Cheadle_tdiff = data.frame(SDdiff$Cheadle,fres)
# F-test says that the following values have statistically the same variance over time in Cheadle
Cheadle_tdiff %>% filter(fres > 0.05)
# and those variances differ over time in Cheadle:
Cheadle_tdiff %>% filter(fres < 0.05)

Newcastle_tdiff = data.frame(SDdiff$Newcastle,fres1)
# F-test says that the following values have statistically the same variance over time in Newcastle
Newcastle_tdiff %>% filter(fres1 > 0.05)
# and those variances differ over time in Newcastle:
Newcastle_tdiff %>% filter(fres1 < 0.05)

Reading_tdiff = data.frame(SDdiff$Reading,fres2)
# F-test says that the following values have statistically the same variance over time in Reading
Reading_tdiff %>% filter(fres2 > 0.05)
# and those variances differ over time in Reading:
Reading_tdiff %>% filter(fres2 < 0.05)
# this shows that the time-dependent variability is different across sites!

# 2) test heteroskedasticity using Breush Pagan and NCV Test
# H0 is that the variance of the residuals is constant for 
# Breush Pagan p > 0.05 is homoscedastic
# NCV p > .05 is homoscedastic
heterosced = data.frame(matrix(ncol = 3, nrow = length(outcome_vars)))
heterosced[,1] = outcome_vars
for (i in 1:length(outcome_vars)){
  heterosced[i,2] = lmtest::bptest(LM[[i]])$p.val
  heterosced[i,3] = car::ncvTest(LM[[i]])$p
}
colnames(heterosced) = c("Var", "BP_Test", "NCV_Test")

# now we can display which metrics are homoscedastic based on both tests' results
heterosced %>% filter(BP_Test > .05 & NCV_Test > .05)
# and which ones are heteroscedastic based on both tests  results
heterosced %>% filter(BP_Test < .05 & NCV_Test < .05)
# >> most metrics are heteroscedastic


#
#
#
#
print("Conclusion variance and heteroscedasticity tests: there is heteroscedasticity across time, site and sex, which needs to be accounted for.")

###### 3.3.2 MODEL FIT ####

# compute AIC, BIC, R^2
diagnostics=list()
for (i in 1:length(outcome_vars)){
  aic = c((AIC(LM[[i]])), (AIC(RI[[i]])),AIC(MM[[i]]))
  bic = c(BIC(LM[[i]]), BIC(RI[[i]]),BIC(MM[[i]]))
  model = c("LM","RI", "MM")
  diagnostics[[i]]=data.frame(model,aic,bic)
}
# find smallest AIC & BIC and highest R squared vals to evaluate model fit
for (i in 1:length(outcome_vars)){
  print(paste("Next Metric:", outcome_vars[i]))
  print(paste("Lowest AIC value in ",diagnostics[[i]][which.min(diagnostics[[i]]$aic),]$model))
  print(paste("Lowest BIC value in ",diagnostics[[i]][which.min(diagnostics[[i]]$bic),]$model))
}
# these analyses show that random intercept models and marginal models outperform vanilla linear models
# but their performance is similar
# ...
# hence, we test the RI and MM models against each other using likelihood ratio tests
chis = c()
ps = c()
for (i in 1:length(outcome_vars)){
  test = lrtest(MM[[i]], RI[[i]])
  chis[i] = test$Chisq[2]
  ps[i] = test$`Pr(>Chisq)`[2]
}
LR_tests = data.frame(outcome_vars, chis, ps)
LR_tests$sig = ifelse(LR_tests$ps < 0.05, "*", "ns")

# inspect model comparisons:
LR_tests
# For SMT trans and smt FA it doesn't matter whether we pick random intercept or marginal models. (Those were the models with higher AIC and BIC for the randon intercepts model)
# For all other metrics marginal models provide better fit.

###### 3.4.3 FINE-TUNING THE MODEL ####
# run marginal models without different assumptions about the variance
MM1 = list()
MM2 = list()
MM3 = list()
for (i in outcome_vars){
  f = formula(paste(i,"~age*sex+site_t3+TP"))
  # no sex var assumption
  MM1[[i]] = gls(model = f3, data=data, correlation = corSymm(form = ~1|eid), weights = varIdent(form = ~1|factor(TP) * factor(site_t3)))
  # no site var assumption
  MM2[[i]] = gls(model = f3, data=data, correlation = corSymm(form = ~1|eid), weights = varIdent(form = ~1|factor(TP) *  factor(sex)))
  # no time point var assumption
  MM3[[i]] = gls(model = f3, data=data, correlation = corSymm(form = ~1|eid), weights = varIdent(form = ~1|factor(site_t3) * factor(sex)))
}

# Is any of the new models anticipating less variability in our confounder factors better in terms of AIC or BIC scores?
for (i in 1:length(outcome_vars)){
  if (AIC(MM[[i]])[1] > AIC(MM1[[i]])[1] || AIC(MM[[i]])[1] > AIC(MM2[[i]])[1] || AIC(MM[[i]])[1] > AIC(MM3[[i]])[1]){
    print(i)
    print(AIC(MM[[i]], MM1[[i]], MM2[[i]], MM3[[i]]))
  }
}
print("These are the diffusion metrics for which simpler models do better:")
print(outcome_vars[6])
print(outcome_vars[8])
print(outcome_vars[10])
print(outcome_vars[25])

# Likelihood Ratio tests for the different models
MMcomp = data.frame(matrix(ncol=6,nrow=length(outcome_vars)))
for (i in outcome_vars){
  test = lrtest(MM[[i]], MM1[[i]])
  test2 = lrtest(MM[[i]], MM2[[i]])
  test3 = lrtest(MM[[i]], MM3[[i]])
  MMcomp[i,1] = test$Chisq[2]
  MMcomp[i,2] = test$`Pr(>Chisq)`[2]
  MMcomp[i,3] = test2$Chisq[2]
  MMcomp[i,4] = test2$`Pr(>Chisq)`[2]
  MMcomp[i,5] = test3$Chisq[2]
  MMcomp[i,6] = test3$`Pr(>Chisq)`[2]
}
LR_tests = data.frame(outcome_vars, na.omit(MMcomp))
# inspect model comparisons:
LR_tests[6,]
LR_tests[8,]
LR_tests[10,]
LR_tests[25,]
# those are clear differences

# we could now go into specifics and check how coefficients compare across models:
# (the effect of TP misrepresented for 6, 8, and 25, as the slope's direction is inconsistent with the age slope's direction)
summary(MM[[6]])$coefficients
summary(MM1[[6]])$coefficients
summary(MM2[[6]])$coefficients
summary(MM3[[6]])$coefficients

summary(MM[[8]])$coefficients
summary(MM1[[8]])$coefficients
summary(MM2[[8]])$coefficients
summary(MM3[[8]])$coefficients

summary(MM[[10]])$coefficients # here, this effect makes more sense, as the metric decreases over time
summary(MM1[[10]])$coefficients
summary(MM2[[10]])$coefficients
summary(MM3[[10]])$coefficients

summary(MM[[25]])$coefficients
summary(MM1[[25]])$coefficients
summary(MM2[[25]])$coefficients
summary(MM3[[25]])$coefficients

# However, for simplicity, we can go with the standard model established in the beginning of the script.


### 4) PLOT ####
# this is to plot the effect of age and time point on diffusion metrics


# create a plot function for standardized betas
## NOTE:
## ylim needs to be adapted for optimal visualisation
## A name vector for the variable names should be passed into the function for better Outcome variable names
## In that case, the ylab() can also be dropped
bplot = function(Beta_Table){
  ggplot(Beta_Table, aes(x=Outcome, y=Std.Beta, color=Outcome)) + 
    geom_pointrange(aes(ymin=Std.Beta-Std.Err, ymax=Std.Beta+Std.Err)) + theme_bw() + theme(legend.position = "none") +
    ylab("Standardized Beta of Age") + xlab("") + ylim(-.075,.075) + scale_x_discrete(limits = rev(factor(var_labels)))+ coord_flip()
}
plot1=bplot(MM_betas)
# # unstandardized betas
# bplot = function(Beta_Table){
#   ggplot(Beta_Table, aes(x=outcome_vars, y=UnStd.Beta, color=outcome_vars)) + 
#     geom_pointrange(aes(ymin=UnStd.Beta-Std.Err, ymax=UnStd.Beta+Std.Err)) + theme_bw() + theme(legend.position = "none") +
#     ylab("") + xlab("Standardized Beta") + ylim(-.004,.004)+ coord_flip() 
# }
# bplot(MM_unstd_betas)
# 
# # Alternatively, we can plot all different models' standardized age coefficients:
# # This shows that random intercept and marginal models produce slighlty smaller effects than the vanilla linear models.
# # However, overall the predictor structures are very similar.
# bplot = function(Beta_Table){
#   ggplot(Beta_Table, aes(x=Outcome, y=Std.Beta, color=Outcome)) +
#     geom_pointrange(aes(ymin=Std.Beta-Std.Err, ymax=Std.Beta+Std.Err)) + theme_bw() + theme(legend.position = "none") +
#     ylab("") + xlab("Standardized Beta") + ylim(-.3,.3) + scale_x_discrete(limits = rev(factor(var_labels)))+ coord_flip()
# }
# # create one plot for each model (linear regression, random intercept, marginal model)
# LM_plot = bplot(LM_betas)
# RI_plot = bplot(RI_betas)
# MM_plot = bplot(MM_betas)
# 
# # merge plots
# figure = ggarrange(LM_plot, RI_plot, MM_plot, ncol = 3)
# # anotate figure
# annotate_figure(figure, #left = textGrob("Variables", rot = 90, vjust = 1, gp = gpar(cex = 1.3)),
#                 bottom = textGrob("Standardized Beta", gp = gpar(cex = 1.3)))
# 


# compile standardized beta values for a variable specified in ROW_NUMBER, e.g. for age ROW_NUMBER = 1
ROW_NUMBER = 5

# change the ROW_NUMBER value accordingly!

# estimate betas, standard errors, T-values and p-values
LMb = c()
LMstde = c()
LMt = c()
LMp = c()
RIb = c()
RIstde = c()
MMb = c()
MMstde = c()
MMt = c()
MMp = c()
for (i in 1:length(outcome_vars)){
  LMb[i] = as.numeric(LMstd[[i]]$coefficients[ROW_NUMBER+1])
  LMstde[i] = as.numeric(sqrt(diag(vcov(LMstd[[i]])))[ROW_NUMBER+1])
  LMt[i] = as.numeric(coef(summary(LMstd[[i]]))[, "t value"][ROW_NUMBER+1])
  LMp[i] = as.numeric(coef(summary(LMstd[[i]]))[, "Pr(>|t|)"][ROW_NUMBER+1])
  RIb[i] = as.numeric(fixef(RIstd[[i]])[ROW_NUMBER-1])
  RIstde[i] = summary(RIstd[[i]])$coefficients[ROW_NUMBER-1,2]
  MMb[i] = as.numeric(MMstd[[i]]$coefficients[ROW_NUMBER+1])
  MMstde[i] = coef(summary(MMstd[[i]]))[ROW_NUMBER+1,2]
  MMt[i] = as.numeric(coef(summary(MMstd[[i]]))[, "t-value"][ROW_NUMBER+1])
  MMp[i] = as.numeric(coef(summary(MMstd[[i]]))[, "p-value"][ROW_NUMBER+1])
}

# make data frames
LM_betas = data.frame(LMb, LMstde,LMt,LMp, var_labels)
RI_betas = data.frame(RIb, RIstde, var_labels)
MM_betas = data.frame(MMb, MMstde, MMt, MMp, var_labels)
colnames(LM_betas) = c("Std.Beta", "Std.Err", "T.val", "p.val", "Outcome")
colnames(RI_betas) =  c("Std.Beta", "Std.Err", "Outcome")
colnames(MM_betas) =  c("Std.Beta", "Std.Err", "T.val", "p.val", "Outcome")
bplotstd = function(Beta_Table){
  ggplot(Beta_Table, aes(x=Outcome, y=Std.Beta, color=outcome_vars)) + 
    geom_pointrange(aes(ymin=Std.Beta-Std.Err, ymax=Std.Beta+Std.Err)) + theme_bw() + theme(legend.position = "none") +
    ylab("Standardized Beta of Time Point") + xlab("") + ylim(-0.2,0.2) + scale_x_discrete(limits = rev(factor(var_labels)))+ coord_flip()
}
#bplotstd(LM_betas)
#bplotstd(RI_betas)
plot2 = bplotstd(MM_betas)
AgeCog = ggarrange(plot1, plot2, ncol = 2)
ggsave("DiffusionAgeAssocitation.pdf", width = 12, height = 9, AgeCog)

# Holm correction
MM_betas$p.adj = p.adjust(MM_betas$p.val,method = "holm")
# Write the table with p-vals for the effect of time point on cognitive measures
p_table = MM_betas %>% select(Outcome,Std.Beta,Std.Err, T.val, p.val, p.adj)
write.csv(p_table, "P_Table_Time_point_DIFF.csv")
