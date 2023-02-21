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
#### 4.1) Age-cognition associations              #
#### 4.2) Ageing & cognitive changes              #
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
T1 = read.csv("TP1_cog.csv")
T2 = read.csv("TP2_cog.csv")

# read in MRI data to get info on sex, site and scanner site
MRI1 = read.csv("T1_noDRAD.csv")
MRI2 = read.csv("T2_noDRAD.csv")
MRI1 = MRI1 %>% dplyr::select(eid, age, sex, site_t3)
MRI2 = MRI2 %>% dplyr::select(eid, age, sex, site_t3)

# create time point dummy TP
T1$TP = replicate(nrow(T1),0)
T2$TP = replicate(nrow(T1),1)

# add sex, age, site to data frames
T1 = merge(T1,MRI1, by = "eid")
T2 = merge(T2,MRI2, by = "eid")
T1$sex = factor(T1$sex)
T1$site_t3 = factor(T1$site_t3)
T1$eid = factor(T1$eid)
T1$TP = factor(T1$TP)
T2$sex = factor(T2$sex)
T2$site_t3 = factor(T2$site_t3)
T2$eid = factor(T2$eid)
T2$TP = factor(T2$TP)


# make single df
data = rbind(T1, T2)
data$AgeDiff = c(replicate(nrow(T1),0), T2$age - T1$age)
data$TP = factor(data$TP)
levels(data$TP) = c("baseline","repeat")
data$TP = factor(data$TP,levels=rev(levels(data$TP)))

############### #
### MISSINGNESS #
############### #

# show missingness by column and time point
for (i in 1:ncol(data)){
  print(colnames(data)[i])
  print(paste("Percentage of total N missing:", sum(is.na(data[i])/nrow(data)*100)))
  print(paste("Time point 1 percentage of total N missing: ", sum(is.na(T1[i]))/nrow(T1)*100))
  print(paste("Time point 2 percentage of total N missing: ", sum(is.na(T2[i]))/nrow(T2)*100))
}
# we should remove variables with more than 25% missingness in order to being able to exectute multiple imputation.
# This would leave us with the following variables:
colnames(data %>% dplyr::select(-inc_pair_matches_r3,-health_selfr,-contains("Matrix_RT"), 
                       -cor_Matrix_puzzles,-viewed_Matrix_puzzles, -cor_Tower_puzzles, -cor_SymDig_matches))
# This is clearly not enough, and we will simply go with row-wise exclusions for each outcome variable.

## 3) ANALYSE ####
#################### #
#### 3.1) run models #### 
#################### #

# make name list of outcome variables, and empty lists (to be filled with fitted models) 
outcome_vars = colnames(data[3:28])
LM = list()
MM = list()
RI = list()
# loop over outcome_vars and run the same LM, RI and MM models for each outcome var.
for (i in outcome_vars){
  f = formula(paste(i,"~age*sex + site_t3 + TP"))
  f2 = formula(paste(i,"~age*sex + TP+(1|eid)"))
  f3 = formula(paste(i,"~age*sex+site_t3+TP"))
  LM[[i]] = lm(f,data=data)
  RI[[i]] = lmer(f2, data = data)
  MM[[i]] = gls(model = f3, data=data, correlation = corSymm(form = ~1|eid), weights = varIdent(form = ~1|factor(TP) * 1|factor(sex)*1|factor(site_t3)), na.action = na.omit)
}

# for the marginal model, we run the analysis on standardized coefficients 
# standardise each diffusion metrics (Z-scoring) for later comparison
T1n = T1
T2n = T2
for (i in 2:28){
  T1n[,i] = (T1[,i]-(mean(T1[,i], na.rm = T)))/(sd(T1[,i], na.rm = T))
  T2n[,i] = (T2[,i]-(mean(T2[,i], na.rm = T)))/(sd(T2[,i], na.rm = T))
}

# create time point dummy TP
T1n$TP = replicate(nrow(T1),0)
T2n$TP = replicate(nrow(T1),1)

# make single df
data2 = rbind(T1n, T2n)
data2$eid = data$eid
data2$AgeDiff = c(replicate(nrow(T1),0), T2$age - T1$age)

# run MM on standardized data
LMstd = list()
RIstd = list()
MMstd = list()
for (i in outcome_vars){
  f = formula(paste(i,"~age*sex + site_t3 + TP"))
  f2 = formula(paste(i,"~age*sex + TP+(1|eid)"))
  f3 = formula(paste(i,"~age*sex+site_t3+TP"))
  LMstd[[i]] = lm(f,data=data2)
  RIstd[[i]] = lmer(f2, data = data2)
  MMstd[[i]] = gls(model = f3, data=data2, correlation = corSymm(form = ~1|eid), weights = varIdent(form = ~1|factor(TP) * 1|factor(sex)*1|factor(site_t3)), na.action = na.omit)
}

#################### #
#### 3.2) get betas  #### 
#################### #

### STANDARDIZED BETAS

# compile standardized beta values for a variable specified in ROW_NUMBER, e.g. for age ROW_NUMBER = 1
ROW_NUMBER = 1

# change the ROW_NUMBER value accordingly!

# estimate betas and standard errors
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
var_labels = c("Incorrect Pair Matches Round 1","Incorrect Pair Matches Round 2",
               "Incorrect Pair Matches Round 3","Maximum Number of Digits Remembered",
               "Fluid Intelligence", "Prospective Memory", "Self-Rated Health", 
               "Response Time Matrix Puzzle 1", "Response Time Matrix Puzzle 2", 
               "Response Time Matrix Puzzle 3", "Response Time Matrix Puzzle 4", 
               "Response Time Matrix Puzzle 5", "Response Time Matrix Puzzle 6", 
               "Response Time Matrix Puzzle 7", "Response Time Matrix Puzzle 8", 
               "Response Time Matrix Puzzle 9", "Response Time Matrix Puzzle 10", 
               "Response Time Matrix Puzzle 11", "Response Time Matrix Puzzle 12", 
               "Response Time Matrix Puzzle 13", "Response Time Matrix Puzzle 14", 
               "Response Time Matrix Puzzle 15", "Solved Matrix Puzzles", "Viewed Matrix Puzzles",
               "Solved Tower Puzzles", "Solved Symbol Digit Matches")
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
write.csv(p_table, "P_Table_Age.csv")

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
  RIb[i] = as.numeric(fixef(RI[[i]])[ROW_NUMBER-1])
  RIstde[i] = summary(RI[[i]])$coefficients[ROW_NUMBER-1,2]
  MMb[i] = as.numeric(MM[[i]]$coefficients[ROW_NUMBER+1])
  MMstde[i] = coef(summary(MM[[i]]))[ROW_NUMBER+1,2]
}

# make data frames
LM_unstd_betas = data.frame(LMb, LMstde, var_labels)
RI_unstd_betas = data.frame(RIb, RIstde, var_labels)
MM_unstd_betas = data.frame(MMb, MMstde, var_labels)
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
df01 = data.frame(apply(df01,2,sd, na.rm = T))
colnames(df01) = c("SD")
df02 = data.frame(apply(df02,2,sd, na.rm = T))
colnames(df02) = c("SD")
df03 = df01-df02
rownames(df03) = rownames(df02)
# display difference: variability increases
df03
print("Ergo, the variability in cognitive scores mostly increases over time")

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
  varres[[i]] = data.frame(apply(varlist[[i]],2,sd, na.rm = T))
  colnames(varres[[i]]) = ("SD")
}
SDdiff = data.frame(varres[[1]]-varres[[2]], varres[[3]]-varres[[4]])
colnames(SDdiff) = c("female", "male")
a = rownames(SDdiff %>% filter(female > 0 & male > 0))
paste("The variability in",a, "score DEcreases over time and across sex")
a = rownames(SDdiff %>% filter(female < 0 & male < 0))
paste("The variability in",a, "score INcreases over time and across sex")
print("Ergo, the variability in cognitive scores changes over time and across sexes")

# looking only at sex
df01 = filter(data, sex==0)
df02 = filter(data, sex==1)
varlist = list(df01, df02)
varres = list()
for (i in 1:length(varlist)){
  varres[[i]] = data.frame(apply(varlist[[i]],2,sd, na.rm = T))
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
## Many of the females' cognitive scores' variance change over time!

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
# For males, variability seems to be more stable than for females across time

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
  varres[[i]] = data.frame(apply(varlist[[i]],2,sd, na.rm = T))
  colnames(varres[[i]]) = ("SD")
}
SDdiff = data.frame(varres[[1]]-varres[[2]], varres[[3]]-varres[[4]], varres[[5]]-varres[[6]])
colnames(SDdiff) = c("Cheadle", "Newcastle", "Reading")
a = rownames(SDdiff %>% filter(Cheadle > 0 & Newcastle > 0 & Reading > 0))
paste("The variability in",a, "score DEcreases over time and across sites")
a = rownames(SDdiff %>% filter(Cheadle < 0 & Newcastle < 0 & Reading < 0))
paste("The variability in",a, "score INcreases over time and across sites")
print("Ergo, the variability in only a few cognitive scores changes over time and across sites")

# looking only at scanner site
df01 = filter(data, site_t3=="Cheadle")
df02 = filter(data, site_t3=="Newcastle")
df03 = filter(data, site_t3=="Reading")
varlist = list(df01, df02, df03)
varres = list()
for (i in 1:length(varlist)){
  varres[[i]] = data.frame(apply(varlist[[i]],2,sd, na.rm = T))
  colnames(varres[[i]]) = ("SD")
}
SDdiff = data.frame(varres[[1]],varres[[2]], varres[[3]])
colnames(SDdiff) = c("Cheadle", "Newcastle", "Reading")
SDdiff = SDdiff[3:28,]
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
  #fres2[i] = var.test(f,data=df03)$p.val
}
# WE CANNOT TEST THE VARIANCE BETWEEN TIME POINTS FOR ALL READING VARIABLES, AS THERE IS TOO LITTLE DATA
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

# again, Reading is not on the list due to missingness

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
# For a few cognitive scores it doesn't matter whether we pick random intercept or marginal models. (Those were the models with higher AIC and BIC for the random intercepts model)
# For all other metrics marginal models provide better fit.

###### 3.4.3 FINE-TUNING THE MODEL ####
# run marginal models without different assumptions about the variance
MM1 = list()
MM2 = list()
MM3 = list()
for (i in outcome_vars){
  f = formula(paste(i,"~age*sex+site_t3+TP"))
  # no sex var assumption
  MM1[[i]] = gls(model = f, data=data, correlation = corSymm(form = ~1|eid), weights = varIdent(form = ~1|factor(TP) * factor(site_t3)), na.action = na.omit)
  # no site var assumption
  MM2[[i]] = gls(model = f, data=data, correlation = corSymm(form = ~1|eid), weights = varIdent(form = ~1|factor(TP) *  factor(sex)), na.action = na.omit)
  # no time point var assumption
  MM3[[i]] = gls(model = f, data=data, correlation = corSymm(form = ~1|eid), weights = varIdent(form = ~1|factor(site_t3) * factor(sex)), na.action = na.omit)
}

# Is any of the new models which are all anticipating more  variability in our confounder factors better in terms of AIC or BIC scores?
for (i in 1:length(outcome_vars)){
  if (AIC(MM[[i]])[1] > AIC(MM1[[i]])[1] || AIC(MM[[i]])[1] > AIC(MM2[[i]])[1] || AIC(MM[[i]])[1] > AIC(MM3[[i]])[1]){
    print(i)
    print(AIC(MM[[i]], MM1[[i]], MM2[[i]], MM3[[i]]))
  }
}
print("Simpler models provide slightly better fit, e.g. only anticipating the variability of sex and time points seems a good choice")
print("This makes sense, as the factor scanner site is most affected of the factors.")

# Likelihood Ratio tests for the different models
MMcomp = data.frame(matrix(ncol=6,nrow=length(outcome_vars)))
for (i in 1:length(outcome_vars)){
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
LR_tests
# those are clear differences

# However, for simplicity, we can go with the standard model established in the beginning of the script.


### 4) PLOT ####
##### 4.1) Age-cognition associations ####
# create a plot function
## NOTE:
## ylim needs to be adapted for optimal visualisation
## A name vector for the variable names should be passed into the function for better Outcome variable names
## In that case, the ylab() can also be dropped

# standardized betas
bplotstd = function(Beta_Table){
  ggplot(Beta_Table, aes(x=Outcome, y=Std.Beta, color=outcome_vars)) + 
    geom_pointrange(aes(ymin=Std.Beta-Std.Err, ymax=Std.Beta+Std.Err)) + theme_bw() + theme(legend.position = "none") +
    ylab("Standardized Beta of Age") + xlab("") + ylim(-0.075,0.04) + scale_x_discrete(limits = rev(factor(var_labels)))+ coord_flip()
}
plot1 = bplotstd(MM_betas)
ggsave("CogAgeAssocitation.pdf", width = 7, height = 9, plot1)

# 
# # unstandardized betas, if of interest
# bplot = function(Beta_Table){
#   ggplot(Beta_Table, aes(x=outcome_vars, y=UnStd.Beta, color=outcome_vars)) + 
#     geom_pointrange(aes(ymin=UnStd.Beta-Std.Err, ymax=UnStd.Beta+Std.Err)) + theme_bw() + theme(legend.position = "none") +
#     ylab("Unstandardized Beta") + xlab("") + ylim(-0.75,2.75)+ coord_flip() 
# }
# bplot(MM_unstd_betas)
# 
# Alternatively, we can plot all different models' standardized coefficients: e.g. of age
# This shows that random intercept and marginal models produce slighlty smaller effects than the vanilla linear models.
# However, overall the predictor structures are very similar.
# bplot = function(Beta_Table){
#   ggplot(Beta_Table, aes(x=outcome_vars, y=Std.Beta, color=outcome_vars)) + 
#     geom_pointrange(aes(ymin=Std.Beta-Std.Err, ymax=Std.Beta+Std.Err)) + theme_bw() + theme(legend.position = "none") +
#     ylab("") + xlab("Standardized Beta") + ylim(-.1,.1)+ coord_flip() 
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


# # check p-vals between time points, if wanted.
# for (i in 3:28){
#   print(outcome_vars[i-2])
#   print((t.test(T1[i], T2[i], use="complete.obs"))$estimate)
#   print((t.test(T1[i], T2[i], use="complete.obs"))$p.value)
# }
