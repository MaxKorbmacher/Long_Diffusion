# regional WM ageing
# Max Korbmacher, 13.02.2022

# 1) PREP ####
## 1.1) PREP PACKAGES ####
# load packages and install if not already installed
remove.packages("sjPlot")
if (!require("pacman")) install.packages("pacman")
pacman::p_load(lme4, nlme, ggplot2, tidyverse, lm.beta, remotes, ggpubr, grid, lmtest, car, lmtest,lmeInfo,lmerTest,sjstats,effsize,Rmpfr,ggrepel,PASWR2)
#install_version("sjstats","0.17.7")
library(sjstats)

## 1.2) PREP DATA ####
# set working directory where the data files are situated
setwd("/home/max/Documents/Projects/Diffusion/UKBlong/data_export/")
# load data (without DRAD extra outliers removed as there were many outliers)
T1 = read.csv("T1_noDRAD.csv")
T2 = read.csv("T2_noDRAD.csv")
# hence, we have to exclude DRAD extra from our analysis, and focus on global average scores here
T1 = T1 %>% dplyr::select(-contains("mean")) %>% dplyr::select(-contains("Drad_extra_")) %>% dplyr::select(-starts_with("X"))
T2 = T2 %>% dplyr::select(-contains("mean")) %>% dplyr::select(-contains("Drad_extra_")) %>% dplyr::select(-starts_with("X"))

# make unstandardized data frame prior to standardization
T1_unstd = T1
T2_unstd = T2
data_unstd = rbind(T1, T2)

# standardise each diffusion metrics (Z-scoring) for later comparison
T1 = T1%>%mutate_if(is.numeric,scale)
T2 = T2%>%mutate_if(is.numeric,scale)

# create time point dummy TP
T1$TP = replicate(nrow(T1),0)
T2$TP = replicate(nrow(T1),1)
T2$time = T2$age - T1$age
T1$time = 0

# make single df and add factor order for time point 
data = rbind(T1, T2)
data$TP = factor(data$TP)
levels(data$TP) = c("baseline","repeat")
data_unstd$TP = factor(data$TP)
levels(data_unstd$TP) = c("baseline","repeat")

# 2) ANALYSE REGIONAL AGEING ####
# make eid last column, where also age, sex, site and time point are located
data = data%>% relocate(eid, .after = last_col())
data_unstd = data_unstd%>% relocate(eid, .after = last_col())
T1_unstd = T1_unstd%>% relocate(eid, .after = last_col())
T2_unstd = T2_unstd%>% relocate(eid, .after = last_col())

# columns 1-1836 are diffusion metrics
#colnames(data[1836:1843])
# these columns will be our outcome variables
outcome_vars = colnames(data[1:1836])

#################### #
#### 2.0) paired samples t.tests for time point differences #### 
#################### #
p.val = c()
t.val = c()
z = c()
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
  # some of the p values will be rounded to 0. We need to use a more accurate approach (see below)
  t.val[i] = t.test(f, paired = TRUE)$statistic
}
# we use Rmpfr to estimate accurate p-vals
.N <- function(.) mpfr(., precBits = 200)
p.val = 2 * pnorm(-abs(.N(t.val)))

for (i in 1:length(outcome_vars)){
  mean_TP1[i] = mean(T1_unstd[,i])
  sd_TP1[i] = sd(T1_unstd[,i])
  mean_TP2[i] = mean(T2_unstd[,i])
  sd_TP2[i] = sd(T2_unstd[,i])
  d[i] = cohen.d(T2_unstd[,i], T1_unstd[,i], paired=TRUE)$estimate
  d_lower[i] = as.numeric(cohen.d(T2_unstd[,i], T1_unstd[,i], paired=TRUE)$conf.int[1])
  d_upper[i] = as.numeric(cohen.d(T2_unstd[,i], T1_unstd[,i], paired=TRUE)$conf.int[2])
}
p.adj = mpfr(p.val*length(outcome_vars),100)
paired_t = data.frame(outcome_vars, mean_TP1, sd_TP1, mean_TP2, sd_TP2, t.val, formatMpfr(p.val),formatMpfr(p.adj), d,d_lower, d_upper)
#paired_t_out = data.frame(lapply(paired_t_out,function(x) if(is.numeric(x)) round(x, 3) else x))
paired_t_out = paired_t
colnames(paired_t_out) = c("Metric", "Mean TP1","SD TP1",  "Mean TP2","SD TP1", "T", "p", "adjusted p", "Cohens's d", "Lower 95% CI", "Upper 95% CI")
# write the table
write.csv(paired_t_out, "/home/max/Documents/Projects/Diffusion/UKBlong/Results/REGIONAL_paired_t_results.csv")
# make a figure summarizing the findings

## add for this -log10 adjusted p-values and labels for colour
# add log transformed p-vals
logp = mpfr(-log10(p.adj),100)
paired_t$logp = formatMpfr(logp)
# add a column of directionality of effect
paired_t$Effect <- "No time point difference"
paired_t$Effect[paired_t$d > 0 & as.numeric(paired_t$logp) > -log10(0.05)] <- "Increase"
paired_t$Effect[paired_t$d < 0 & (paired_t$logp) > -log10(0.05)] <- "Decrease"
# add cleaner text for labels which will be printed
paired_t$outcome_vars = gsub("radEAD_Fornix_StriaterminalisR", "WMTI - radEAD Fornix-StTer(R)",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("radEAD_Fornix", "WMTI - radEAD Fornix",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("FA_Fornix_StriaterminalisR", "DTI - FA Fornix-StTer(R)",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("AD_Fornix", "DTI - AD Fornix",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("AD_Middlecerebellarpeduncle", "DTI - AD MCP",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("MD_Middlecerebellarpeduncle", "DTI - MD MCP",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("FA_Fornix", "DTI - FA Fornix",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("FA_FMIN", "DTI - FA FMIN",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("FA_BodyCC", "DTI - FA BodyCC",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("RD_BodyCC", "DTI - RD BodyCC",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("MD_Fornix", "DTI - MD Fornix",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("RD_Fornix", "DTI - RD Fornix",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("RD_Fornix_StriaterminalisR", "DTI - RD Fornix-StTer(R)",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("rk_Fornix", "DKI - RK Fornix",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("mk_Fornix", "DKI - MK Fornix",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("ak_Fornix", "DKI - AK Fornix",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("micro_FA_Fornix", "BRIA - micro FA Fornix",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("v_intra_Fornix", "BRIA - v intra Fornix",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("v_extra_Fornix", "BRIA - v extra Fornix",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("micro_ADC_Fornix", "BRIA - micro ADC Fornix",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("v_csf_Fornix", "BRIA - v csf Fornix",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("v_csf_Middlecerebellarpeduncle", "BRIA - v csf MCP",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("micro_Rd_Middlecerebellarpeduncle", "BRIA - micro RD MCP",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("micro_Rd_Fornix", "BRIA - micro RD Fornix",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("micro_Ax_Fornix", "BRIA - micro AX Fornix",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("micro_Ax_SLTFR", "BRIA - micro AX SLTF(R)",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("micro_Ax_Middlecerebellarpeduncle", "BRIA - micro AX MCP",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("Dax_intra_Middlecerebellarpeduncle", "BRIA - DAX intra MCP",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("Dax_extra_Middlecerebellarpeduncle", "BRIA - DAX extra MCP",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("micro_ADC_Middlecerebellarpeduncle", "BRIA - micro ADC MCP",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("smt_long_Middlecerebellarpeduncle", "SMT - long MCP",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("smt_trans_SuperiorcerebellarpeduncleR", "SMT - trans SCP(R)",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("smt_fa_SuperiorcerebellarpeduncleR", "SMT - FA SCP(R)",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("smt_md_InferiorcerebellarpeduncleR", "SMT - MD ICP(R)",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("smt_md_Middlecerebellarpeduncle", "SMT - MD MCP",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("smt_long_InferiorcerebellarpeduncleR", "SMT - long ICP(R)",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("smt_mc_intra_Fornix", "SMTmc - intra Fornix",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("smt_mc_diff_Middlecerebellarpeduncle", "SMTmc - MD MCP",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("smt_mc_extramd_Fornix", "SMTmc - extra MD Fornix",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("smt_mc_extramd_Middlecerebellarpeduncle", "SMTmc - extra MD MCP",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("smt_mc_extratrans_Fornix", "SMTmc - extra trans Fornix",paired_t$outcome_vars)
paired_t$outcome_vars = gsub("awf_Fornix", "WMTI - AWF Fornix",paired_t$outcome_vars)
# plot
plot = ggplot(data=paired_t, aes(x=d, y=as.numeric(logp), label=outcome_vars, col = Effect)) +
  geom_point() + 
  theme_minimal() +
  # make labels based on Cohen's |d| > 0.5
  geom_text_repel(data = subset(paired_t, d < -0.5),colour='black', nudge_x = -1.7, direction = "y", segment.size = 0.1, xlim = c(-0.3,-1.6))+
  geom_text_repel(data = subset(paired_t, d > 0.5),colour='black', nudge_x = 1, direction = "y",segment.size = 0.1, xlim = c(0.2,1.3))+
  # make labels based on extreme p-values
  geom_text_repel(data = subset(paired_t, as.numeric(logp) > 500 & d < 0),colour='black', nudge_x = -1.7, direction = "y", segment.size = 0.1, xlim = c(-0.3,-1.6))+
  geom_text_repel(data = subset(paired_t, as.numeric(logp) > 500 & d > 0),colour='black', nudge_x = 1, direction = "y",segment.size = 0.1, xlim = c(0.2,1.3))+
  # geom_text_repel(
  #   force        = 0.25,
  #   #nudge_x      = -0.25,
  #   #direction    = "y",
  #   #hjust        = 1,
  #   #segment.size = 0.1,
  #   max.overlaps = Inf,
  #   colour='black')+
  #scale_fill_brewer(palette="Dark2")
  #scale_color_manual(values=c("#56B4E9","#D55E00", "#999999", "#0072B2", "#E69F00"))+
  #scale_color_manual(values=c("#0072B2","#D55E00", "#999999", "#56B4E9", "#E69F00"))+
  scale_color_manual(values=c("#56B4E9","#D55E00", "#999999"))+
  #scale_color_manual(values=c("blue","red", "black", "orange", "purple"))+
  geom_vline(xintercept = c(-0.5, 0.5), color = "red", linetype = "dashed", cex = 1, alpha = 0.2) +
  geom_hline(yintercept=c(-log10(0.05), 500), color = "red", linetype = "dashed", cex = 1, alpha = 0.2) +
  #geom_hline(yintercept= 500, col="black") + 
  #geom_hline(yintercept=(600), col="black") + 
  xlab("")+
  ylab("-log10(Bonferroni-corrected p)")+
  xlim(-1.7,1.5)+
  theme(legend.position="bottom")#+labs(title = "Time point changes in WM microstructure")
ggsave("/home/max/Documents/Projects/Diffusion/UKBlong/Results/volcano.pdf",plot, height = 7, width = 15)

# count how many were significant
prop.table(table(paired_t$Effect))

#################### #
#### 2.1) run linear random effects models #### 
#################### #
# make sure that TP2 -TP1 to see the right direction of developments 
data$TP = ordered(data$TP, levels =  c("repeat","baseline"))
data$sex = factor(data$sex)
levels(data$sex) = c("female", "male")
# make empty list to be filled
RI = list()
# loop over outcome_vars and run the same LM, RI and MM models for each outcome var.
for (i in outcome_vars){
  f = formula(paste(i,"~age*sex + TP+(1|site_t3/eid)"))
  RI[[i]] = lmer(f, data = data)
}

#################### ###
#### 2.2) get betas  #### 
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
  RI_betas = data.frame(outcome_vars, RIb, RIstde, RIp)
  # correct labelling of column names
  colnames(RI_betas) =  c("Outcome","Std.Beta", "Std.Err", "p")
  # Holm correction
  RI_betas$p.adj = p.adjust(RI_betas$p,method = "bonferroni")
  RI_betas$Model = diffusion_models
  betas[[ROW_NUMBER]] = RI_betas
  #betas[[ROW_NUMBER]] = data.frame(lapply(RI_betas,function(x) if(is.numeric(x)) round(x, 3) else x)) # this line can be used for rounding
}

# save  
write.csv(betas[[1]], "/home/max/Documents/Projects/Diffusion/UKBlong/Results/REGIONAL_effect_of_age.csv")
write.csv(betas[[2]], "/home/max/Documents/Projects/Diffusion/UKBlong/Results/REGIONAL_effect_of_sex.csv")
write.csv(betas[[3]], "/home/max/Documents/Projects/Diffusion/UKBlong/Results/REGIONAL_effect_of_timepoint.csv")
write.csv(betas[[4]], "/home/max/Documents/Projects/Diffusion/UKBlong/Results/REGIONAL_effect_of_age_sex_interaction.csv")

# we plot AGE effects (betas and p-vals)
dmri = betas[[1]]
# we use log-transformed p-values for visualisation purposes
dmri$logp = -log(dmri$p.adj)
# first, add a column of directionality for beta values / slopes
dmri$Slope <- "No age effect"
dmri$Slope[dmri$Std.Beta > 0 & dmri$logp > -log10(0.05)] <- "Positive association"
dmri$Slope[dmri$Std.Beta < 0 & dmri$logp > -log10(0.05)] <- "Negative association"

# fix labels
dmri$Outcome = gsub("RD_Fornix_StriaterminalisR", "DTI - RD Fornix-StTer(R)",dmri$Outcome)
dmri$Outcome = gsub("radEAD_Fornix_StriaterminalisR", "WMTI - radEAD Fornix-StTer(R)",dmri$Outcome)
dmri$Outcome = gsub("FA_Fornix_StriaterminalisR", "DTI - FA Fornix-StTer(R)",dmri$Outcome)
dmri$Outcome = gsub("micro_FA_Fornix", "BRIA - micro FA Fornix",dmri$Outcome)
dmri$Outcome = gsub("micro_Rd_Fornix", "BRIA - micro RD Fornix",dmri$Outcome)
dmri$Outcome = gsub("AD_Fornix", "DTI - AD Fornix",dmri$Outcome)
dmri$Outcome = gsub("FA_Fornix", "DTI - FA Fornix",dmri$Outcome)
dmri$Outcome = gsub("FA_FMIN", "DTI - FA FMIN",dmri$Outcome)
dmri$Outcome = gsub("MD_Fornix", "DTI - MD Fornix",dmri$Outcome)
dmri$Outcome = gsub("RD_Fornix", "DTI - RD Fornix",dmri$Outcome)
dmri$Outcome = gsub("rk_Fornix", "DKI - RK Fornix",dmri$Outcome)
dmri$Outcome = gsub("mk_Fornix", "DKI - MK Fornix",dmri$Outcome)
dmri$Outcome = gsub("ak_Fornix", "DKI - AK Fornix",dmri$Outcome)
dmri$Outcome = gsub("v_intra_Fornix", "BRIA - v intra Fornix",dmri$Outcome)
dmri$Outcome = gsub("v_extra_Fornix", "BRIA - v extra Fornix",dmri$Outcome)
dmri$Outcome = gsub("micro_ADC_Fornix", "BRIA - micro ADC Fornix",dmri$Outcome)
dmri$Outcome = gsub("v_csf_Fornix", "BRIA - v csf Fornix",dmri$Outcome)
dmri$Outcome = gsub("micro_Ax_Fornix", "BRIA - micro AX Fornix",dmri$Outcome)
dmri$Outcome = gsub("smt_mc_intra_Fornix", "SMTmc - intra Fornix",dmri$Outcome)
dmri$Outcome = gsub("smt_mc_extramd_Fornix", "SMTmc - extra MD Fornix",dmri$Outcome)
dmri$Outcome = gsub("smt_mc_extratrans_Fornix", "SMTmc - extra trans Fornix",dmri$Outcome)
dmri$Outcome = gsub("awf_Fornix", "WMTI - AWF Fornix",dmri$Outcome)
dmri$Outcome = gsub("radEAD_Fornix", "WMTI - radEAD Fornix",dmri$Outcome)

# then plot (volcano plot)
plot2 = ggplot(data=dmri, aes(x=Std.Beta, y=logp,col = Slope, label=Outcome)) +
  geom_point() + 
  theme_minimal() +
  geom_vline(xintercept = c(-0.5, 0.5), color = "red", linetype = "dashed", cex = 1, alpha = 0.2) +
  geom_hline(yintercept=-log10(0.05), color = "red", linetype = "dashed", cex = 1, alpha = 0.2) +
  geom_text_repel(data = subset(dmri, Std.Beta < -0.5),colour='black', nudge_x = -0.7, direction = "y", segment.size = 0.1, xlim = c(-0.3,-1.2))+
  geom_text_repel(data = subset(dmri, Std.Beta > 0.5),colour='black', nudge_x = 0.7, direction = "y",segment.size = 0.1, xlim = c(0.3,1.2))+
  # geom_text_repel(
  #   force        = 0.25,
  #   #nudge_x      = -0.25,
  #   #direction    = "y",
  #   #hjust        = 1,
  #   #segment.size = 0.1,
  #   max.overlaps = Inf,
  #   colour='black')+
  #scale_fill_brewer(palette="Dark2")
  #scale_color_manual(values=c("#56B4E9","#D55E00", "#999999", "#0072B2", "#E69F00"))+
  scale_color_manual(values=c("#56B4E9","#999999","#D55E00"))+
  #scale_color_manual(values=c("blue","red", "black", "orange", "purple"))+
  #geom_vline(xintercept=c(-0.075, 0.075), col="black") +
  #geom_hline(yintercept=-log10(0.05), col="black") + 
  xlab("")+ylab("-log10(Bonferroni-corrected p)")+
  xlim(-1.2,1.2)+
  theme(legend.position="bottom")#+labs(title = "Associations of MRI features' lateralisation and age")
ggsave("/home/max/Documents/Projects/Diffusion/UKBlong/Results/REGIONAL_volcano_AGE.pdf",plot2, height = 7, width = 11)


# we plot only time point differences (betas and p-vals)
dmri = betas[[2]]
# we use log-transformed p-values for visualisation purposes
dmri$logp = -log(dmri$p.adj)
# first, add a column of directionality for beta values / slopes
dmri$Slope <- "No sex difference"
dmri$Slope[dmri$Std.Beta > 0 & dmri$logp > -log10(0.05)] <- "Increase"
dmri$Slope[dmri$Std.Beta < 0 & dmri$logp > -log10(0.05)] <- "Decrease"

# then plot (volcano plot)
plot2 = ggplot(data=dmri, aes(x=Std.Beta, y=logp,col = Slope, label=Outcome)) +
  geom_point() + 
  theme_minimal() +
  geom_vline(xintercept = c(-0.5, 0.5), color = "red", linetype = "dashed", cex = 1, alpha = 0.2) +
  geom_hline(yintercept=-log10(0.05), color = "red", linetype = "dashed", cex = 1, alpha = 0.2) +
  geom_text_repel(data = subset(dmri, Std.Beta < -0.5),colour='black', nudge_x = -0.7, direction = "y", segment.size = 0.1, xlim = c(-0.3,-1.2))+
  geom_text_repel(data = subset(dmri, Std.Beta > 0.5),colour='black', nudge_x = 0.7, direction = "y",segment.size = 0.1, xlim = c(0.3,1.2))+
  # geom_text_repel(
  #   force        = 0.25,
  #   #nudge_x      = -0.25,
  #   #direction    = "y",
  #   #hjust        = 1,
  #   #segment.size = 0.1,
  #   max.overlaps = Inf,
  #   colour='black')+
  #scale_fill_brewer(palette="Dark2")
  #scale_color_manual(values=c("#56B4E9","#D55E00", "#999999", "#0072B2", "#E69F00"))+
  scale_color_manual(values=c("#56B4E9","#D55E00","#999999"))+
  #scale_color_manual(values=c("blue","red", "black", "orange", "purple"))+
  #geom_vline(xintercept=c(-0.075, 0.075), col="black") +
  #geom_hline(yintercept=-log10(0.05), col="black") + 
  xlab("")+ylab("-log10(Bonferroni-corrected p)")+
  xlim(-1.2,1.2)+
  theme(legend.position="bottom")#+labs(title = "Associations of MRI features' lateralisation and age")
ggsave("/home/max/Documents/Projects/Diffusion/UKBlong/Results/REGIONAL_volcano_SEX.pdf",plot2, height = 7, width = 11)


# we plot only time point differences (betas and p-vals)
dmri = betas[[3]]
# we use log-transformed p-values for visualisation purposes
dmri$logp = -log(dmri$p.adj)
# first, add a column of directionality for beta values / slopes
dmri$Slope <- "No time point difference"
dmri$Slope[dmri$Std.Beta > 0 & dmri$logp > -log10(0.05)] <- "Increase"
dmri$Slope[dmri$Std.Beta < 0 & dmri$logp > -log10(0.05)] <- "Decrease"

# then plot (volcano plot)
plot2 = ggplot(data=dmri, aes(x=Std.Beta, y=logp,col = Slope, label=Outcome)) +
  geom_point() + 
  theme_minimal() +
  geom_vline(xintercept = c(-0.25, 0.25), color = "red", linetype = "dashed", cex = 1, alpha = 0.2) +
  geom_hline(yintercept=-log10(0.05), color = "red", linetype = "dashed", cex = 1, alpha = 0.2) +
  geom_text_repel(data = subset(dmri, Std.Beta < -0.25),colour='black', nudge_x = -0.5, direction = "y", segment.size = 0.1, xlim = c(-0.3,-0.9))+
  geom_text_repel(data = subset(dmri, Std.Beta > 0.25),colour='black', nudge_x = 0.5, direction = "y",segment.size = 0.1, xlim = c(0.3,0.9))+
  # geom_text_repel(
  #   force        = 0.25,
  #   #nudge_x      = -0.25,
  #   #direction    = "y",
  #   #hjust        = 1,
  #   #segment.size = 0.1,
  #   max.overlaps = Inf,
  #   colour='black')+
  #scale_fill_brewer(palette="Dark2")
  #scale_color_manual(values=c("#56B4E9","#D55E00", "#999999", "#0072B2", "#E69F00"))+
  scale_color_manual(values=c("#56B4E9","#D55E00", "#999999"))+
  #scale_color_manual(values=c("blue","red", "black", "orange", "purple"))+
  #geom_vline(xintercept=c(-0.075, 0.075), col="black") +
  #geom_hline(yintercept=-log10(0.05), col="black") + 
  xlab("")+ylab("-log10(Bonferroni-corrected p)")+
  xlim(-0.9,0.9)+
  theme(legend.position="bottom")#+labs(title = "Associations of MRI features' lateralisation and age")
ggsave("/home/max/Documents/Projects/Diffusion/UKBlong/Results/REGIONAL_volcano_TIMEPOINT.pdf",plot2, height = 7, width = 15)

hist(dmri$Std.Beta)

# we plot only time point differences (betas and p-vals)
dmri = betas[[4]]
# we use log-transformed p-values for visualisation purposes
dmri$logp = -log(dmri$p.adj)
# first, add a column of directionality for beta values / slopes
dmri$Slope <- "No sex difference"
dmri$Slope[dmri$Std.Beta > 0 & dmri$logp > -log10(0.05)] <- "Increase"
dmri$Slope[dmri$Std.Beta < 0 & dmri$logp > -log10(0.05)] <- "Decrease"

# then plot (volcano plot)
plot2 = ggplot(data=dmri, aes(x=Std.Beta, y=logp,col = Slope, label=Outcome)) +
  geom_point() + 
  theme_minimal() +
  geom_vline(xintercept = c(-0.25, 0.25), color = "red", linetype = "dashed", cex = 1, alpha = 0.2) +
  geom_hline(yintercept=-log10(0.05), color = "red", linetype = "dashed", cex = 1, alpha = 0.2) +
  geom_text_repel(data = subset(dmri, Std.Beta < -0.25),colour='black', nudge_x = -0.5, direction = "y", segment.size = 0.1, xlim = c(-0.3,-0.9))+
  geom_text_repel(data = subset(dmri, Std.Beta > 0.25),colour='black', nudge_x = 0.5, direction = "y",segment.size = 0.1, xlim = c(0.3,0.9))+
  # geom_text_repel(
  #   force        = 0.25,
  #   #nudge_x      = -0.25,
  #   #direction    = "y",
  #   #hjust        = 1,
  #   #segment.size = 0.1,
  #   max.overlaps = Inf,
  #   colour='black')+
  #scale_fill_brewer(palette="Dark2")
  #scale_color_manual(values=c("#56B4E9","#D55E00", "#999999", "#0072B2", "#E69F00"))+
  scale_color_manual(values=c("#56B4E9","#D55E00", "#999999"))+
  #scale_color_manual(values=c("blue","red", "black", "orange", "purple"))+
  #geom_vline(xintercept=c(-0.075, 0.075), col="black") +
  #geom_hline(yintercept=-log10(0.05), col="black") + 
  xlab("")+ylab("-log10(Bonferroni-corrected p)")+
  xlim(-0.9,0.9)+
  theme(legend.position="bottom")#+labs(title = "Associations of MRI features' lateralisation and age")
ggsave("/home/max/Documents/Projects/Diffusion/UKBlong/Results/REGIONAL_volcano_SEX_AGE_INTERACTION.pdf",plot2, height = 7, width = 15)

