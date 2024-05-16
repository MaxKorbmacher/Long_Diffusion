# visualise relationships between longitudinal changes in white matter and PGRS
# 15.08.2022; Max Korbmacher (max.korbmacher@gmail.com)
# due to model convergence issues, mixed effects were abandoned. Instead, simple linear models were applied for the PGRS-WMM relationships.

# PREP ####

# load packages
if (!require("pacman")) install.packages("pacman")
pacman::p_load(lme4, nlme, ggplot2, tidyverse, lm.beta, remotes, ggpubr, grid, lmtest, car, lmtest,lmeInfo,lmerTest,sjstats,effsize,Rmpfr,ggrepel,PASWR2, reshape2)

# load data
PGRS = read.csv("/cluster/projects/p33/users/maxk/UKB/genetics/PRS.csv")
T1 = read.csv("/cluster/projects/p33/users/maxk/UKB/longitudinal/data/noICD10_T1.csv")
T2 = read.csv("/cluster/projects/p33/users/maxk/UKB/longitudinal/data/T2_noDRAD.csv")
cross_sectional = read.csv("/cluster/projects/p33/users/maxk/UKB/batch2023/final_dMRI_data.csv")
outliers = read.csv("/cluster/projects/p33/users/maxk/UKB/longitudinal/data/outliers.csv")
PGRS_copy = PGRS
# # make copies
# T1c = T1
# T2c = T2
# T1 = T1c
# T2 = T2c

# exclude outliers / impossible values
T1 = T1[!T1$eid %in% outliers$x,]
T2 = T2[!T2$eid %in% outliers$x,]
T2 = T2[T2$eid %in% T1$eid,]
# additionally we have to exclude DRAD extra from our analysis, as this metric did poor in the QC on whole skeleton average values
# the same QC concern is true for axEAD estimated in the Fornix. We hence exclude the metric.
# we also extract only mean / whole skeleton average values in this step of the analysis
T1 = T1 %>% dplyr::select(!(contains("Drad_extra")))%>% dplyr::select(!(contains("axEAD")))
T2 = T2 %>% dplyr::select(!(contains("Drad_extra")))%>% dplyr::select(!(contains("axEAD")))

# merge eids for PGRS and WM / diffusion metrics
T1 = T1[T1$eid %in% PGRS$FID,]
T2 = T2[T2$eid %in% PGRS$FID,]
PGRS = PGRS[PGRS$FID %in% T1$eid,]

# make copy of MRI data frames for later
allT1 = T1
allT2 = T2

# select mean/global averages only
T1 = T1 %>% dplyr::select(contains("mean"), age, sex, site_t3, eid)
T2 = T2 %>% dplyr::select(contains("mean"), age, sex, site_t3, eid)

# now, select only actual PGRS scores for data frame
PGRS1 = PGRS[3:ncol(PGRS)]

## Global Annual Change - PGRS associations ####
# now, we can subtract one time point from the other, the last column is age, for which we estimate the delta
ROC = T2[1:27]-T1[1:27]
lm_list = list()
age = T1$age
sex = T1$sex
site = T1$site_t3
betas = c()
outdat = data.frame(matrix(ncol = ncol(PGRS1), nrow = ncol(ROC)-1))
ps = outdat
for (o in 1:ncol(PGRS1)){
  tmpPGRS = scale(PGRS1[[o]])
  for (i in 1:26){
    tmpROC = scale(ROC[i]/abs(ROC[27]), center = FALSE)
    tmp_df = data.frame(age, sex, site, tmpROC, tmpPGRS)
    tmp_model = lm(tmpROC ~ tmpPGRS + age*sex + site, data = tmp_df)
    outdat[i,o] = summary(tmp_model)$coefficients[2]
    ps[i,o] = summary(tmp_model)$coefficients[2,4]
  }
}

#  rename columns
## for beta values
colnames(outdat) = colnames(PGRS1)
var_labels = c("BRIA - V intra", "BRIA - V extra", "BRIA - V csf", "BRIA - micro RD","BRIA - micro FA",
               "BRIA - micro AX", "BRIA - micro ADC", "BRIA - DAX intra", "BRIA - DAX extra",
               
               "DKI - RK", "DKI - AK", "DKI -MK",
               
               "DTI - FA", "DTI - MD", "DTI - RD", "DTI - AD",
               
               "SMT - FA", "SMT - long", "SMT - MD", "SMT - trans",
               
               "SMTmc - intra","SMTmc - extra MD", "SMTmc - extra trans",  "SMTmc - diff",
               
               "WMTI - AWF", "WMTI - radEAD")
outdat$names = var_labels
## for p values
colnames(ps) = colnames(PGRS1)
ps$names = var_labels

# melt data frame into long format for plotting

plot_dat = reshape2::melt(outdat)
p_frame = melt(ps)


# this is what happens when adjusting the p-values. (No significant associations left after correcting for multiple comparison)
min(p.adjust(p_frame$value, method = "fdr"))

# we could also select findings based on p-values before correction... (this is currently used in the plot)
p_frame %>% select(names, variable, value) %>% dplyr::filter(value < .05)
plot_dat$p = p_frame$value

# to frame only significant associations define the colour of the frames
colors = c(ifelse(plot_dat$p < .05,"black", "white"))

### TILE PLOT for Global Change ####
plot = ggplot(plot_dat, aes(variable, names,fill = value)) +
  geom_tile(lwd = 0.75, linetype = 1, color = colors, width = 0.9, height = 0.9) +  
  geom_text(aes(label=round(value,2)),size=4)+
  scale_fill_gradient(low = "#0072B2", high = "#F0E442") +
  ylab("") + xlab("") + theme_bw() + guides(fill=guide_legend(title="Association Strength"))
ggsave("/tsd/p33/home/p33-maxk/export/PGRS_annual_change_plots.pdf",plot, width = 10, height = 12)
ggsave("/cluster/projects/p33/users/maxk/UKB/export/LONG_WM_change_PGRS.pdf",plot, height = 12, width = 10)

## STATS for GLOBAL ARoC-PGRS associations ####
plot_dat %>% summarise(M = mean(abs(value)), SD = sd(abs(value)), Md= median(abs(value)),MAD = mad(abs(value)))
plot_beta %>% group_by(variable) %>% summarise(M = mean(abs(value)), SD = sd(abs(value)), Md= median(abs(value)),MAD = mad(abs(value)))

## Global Metrics - PGRS associations ####
# scale mean metrics
T1[1:26] = T1[1:26]%>%mutate_if(is.numeric,scale)
T2[1:26] = T2[1:26]%>%mutate_if(is.numeric,scale)

# get PGRS-WM associations at time point 1
outdat1 = data.frame(matrix(ncol = ncol(PGRS1), nrow = (26)))
ps1 = outdat1
for (o in 1:ncol(PGRS1)){
  tmpPGRS = scale(PGRS1[[o]])
  for (i in 1:26){
    tmpROC = T1[[i]]
    age = T1$age
    tmp_df = data.frame(age, sex, site, tmpROC, tmpPGRS)
    tmp_model = lm(tmpROC ~ tmpPGRS + age*sex + (site), data = tmp_df)
    outdat1[i,o] = summary(tmp_model)$coefficients[2]
    ps1[i,o] = summary(tmp_model)$coefficients[2,4]
  }
}

# get PGRS-WM associations at time point 2
outdat2 = data.frame(matrix(ncol = ncol(PGRS1), nrow = (26)))
ps2 = outdat2
for (o in 1:ncol(PGRS1)){
  tmpPGRS = scale(PGRS1[[o]])
  for (i in 1:26){
    tmpROC = T2[[i]]
    age = T2$age
    tmp_df = data.frame(age, sex, site, tmpROC, tmpPGRS)
    tmp_model = lm(tmpROC ~ tmpPGRS + age*sex + (site), data = tmp_df)
    outdat2[i,o] = summary(tmp_model)$coefficients[2]
    ps2[i,o] = summary(tmp_model)$coefficients[2,4]
  }
}

#  rename columns
## for beta values
colnames(outdat1) = colnames(PGRS1)
colnames(outdat2) = colnames(PGRS1)
## add labels
outdat1$names = var_labels
outdat2$names = var_labels

## for p values
colnames(ps1) = colnames(PGRS1)
colnames(ps2) = colnames(PGRS1)
ps1$names = var_labels
ps2$names = var_labels

# melt data frame into long format for plotting
plot_dat1 = melt(outdat1)
p_frame1 = melt(ps1)
plot_dat2 = melt(outdat2)
p_frame2 = melt(ps2)

# this is what happens when adjusting the p-values. (No significant associations left after correcting for multiple comparison)
min(p.adjust(p_frame1$value, method = "fdr"))
min(p.adjust(p_frame2$value, method = "fdr"))
# we could also select findings based on p-values before correction... (this is currently used in the plots)
p_frame1 %>% select(names, variable, value) %>% dplyr::filter(value < .05)
p_frame2 %>% select(names, variable, value) %>% dplyr::filter(value < .05)

# use UNCORRECTED p-values for visualisation
plot_dat1$p = p_frame1$value
plot_dat2$p = p_frame2$value

# to frame only significant associations define the colour of the frames
colors1 = c(ifelse(plot_dat1$p < .05,"black", "white"))
colors2 = c(ifelse(plot_dat2$p < .05,"black", "white"))

### TILE PLOT FOR GLOBAL METRICS AT EACH TIME POINT ####
plot1 = ggplot(plot_dat1, aes(variable, names,fill = value)) +
  geom_tile(lwd = 0.75, linetype = 1, color = colors1, width = 0.9, height = 0.9) +  
  geom_text(aes(label=round(value,2)),size=4)+
  scale_fill_gradient(low = "#0072B2", high = "#F0E442") +
  ylab("") + xlab("") + theme_bw() + guides(fill=guide_legend(title="Association Strength"))
plot2 = ggplot(plot_dat2, aes(variable, names,fill = value)) +
  geom_tile(lwd = 0.75, linetype = 1, color = colors2, width = 0.9, height = 0.9) +  
  geom_text(aes(label=round(value,2)),size=4)+
  scale_fill_gradient(low = "#0072B2", high = "#F0E442") +
  ylab("") + xlab("") + theme_bw() + guides(fill=guide_legend(title="Association Strength"))
plot3g = ggarrange(plot1,plot2, labels = c("a","b"),common.legend = TRUE)
ggsave("/tsd/p33/home/p33-maxk/export/PGRS_T1T2_plots.pdf",plot3g, width = 17, height = 12)
ggsave("/cluster/projects/p33/users/maxk/UKB/export/LONG_WM_T1_T2_PGRS_plots.pdf",plot3g, width = 17, height = 12)

#############################
## Regional Annual Change - PGRS associations ####
# select mean/global averages only
allT1 = allT1 %>% dplyr::select(-c(X, X.1,X.2,site_t4)) %>% dplyr::select(-sex, -site_t3, -eid, -age, everything())
allT2 = allT2 %>% dplyr::select(-c(X,site_t4)) %>% dplyr::select(-sex, -site_t3, -eid, -age, everything())
# number of columns with WM metrics
Ncol = ncol(allT1)-4
ageDiff = allT2$age-allT1$age
# now, estimate the annual regional change and how it relates to PGRS
ROC = allT2[1:Ncol]-allT1[1:Ncol]
lm_list = list()
age = allT1$age
sex = allT1$sex
site = allT1$site_t3
outdat = data.frame(matrix(ncol = ncol(PGRS1), nrow = ncol(ROC)))
ps = outdat
for (o in 1:ncol(PGRS1)){
  tmpPGRS = scale(PGRS1[[o]])
  for (i in 1:ncol(ROC)){
    tmpROC = scale(ROC[i]/abs(ageDiff), center = FALSE)
    tmp_df = data.frame(age, sex, site, tmpROC, tmpPGRS)
    tmp_model = lm(tmpROC ~ tmpPGRS + age*sex + (site), data = tmp_df)
    outdat[i,o] = summary(tmp_model)$coefficients[2]
    ps[i,o] = summary(tmp_model)$coefficients[2,4]
  }
}
### PLOT for Regional Change ####
colnames(outdat)= colnames(PGRS1)
colnames(ps)= colnames(PGRS1)

# adjust p-values(?)
correctedP = ps*(Ncol*ncol(PGRS1))
for (i in colnames(correctedP)){
  print(correctedP %>% filter(i < .05))
}
# this is presented as a question, as there are no p-values surviving the multiple comparison.

# However, we can look at the smallest p-values and largest associations in a volcano plot
plot_beta = melt(outdat)
plot_p = melt(ps)
plot_beta$p = plot_p$value
plot_beta$Outcome = colnames(ROC)

# write the plot_beta data frame for potential later purposes
write.csv(plot_beta, "/cluster/projects/p33/users/maxk/UKB/longitudinal/results/Longitudinal/regional_WMM_ARoC_PGRS_associations.csv")
## One could add some colour to the direction of the effect, but that would require some consideration as there are many different observed groups
# dmri$Slope <- "No age effect"
# dmri$Slope[dmri$Std.Beta > 0 & dmri$logp > -log10(0.05)] <- "Positive association"
# dmri$Slope[dmri$Std.Beta < 0 & dmri$logp > -log10(0.05)] <- "Negative association"

# rename the outcome variables to become nice labels
plot_beta$Outcome = gsub("rk_InferiorcerebellarpeduncleR", "DKI - RK rICP",plot_beta$Outcome)
plot_beta$Outcome = gsub("smt_long_Fornix", "SMT - long Fornix",plot_beta$Outcome)
plot_beta$Outcome = gsub("ak_SuperiorfrontooccipitalfasciculusR", "DKI - AK rSFOF",plot_beta$Outcome)
plot_beta$Outcome = gsub("smt_fa_MediallemniscusR", "SMT - FA rMediallemniscus",plot_beta$Outcome)
plot_beta$Outcome = gsub("Dax_extra_CerebralpeduncleR", "BRIA - extra rCereb.Ped.",plot_beta$Outcome)
plot_beta$Outcome = gsub("ak_UFR", "DKI - AK rUF",plot_beta$Outcome)
plot_beta$Outcome = gsub("smt_trans_CingulumcingulategyrusR", "SMT - trans rCingulum",plot_beta$Outcome)
plot_beta$Outcome = gsub("ak_RetrolenticularpartofinternalcapsuleR", "DKI - AK rRIC",plot_beta$Outcome)
plot_beta$Outcome = gsub("v_intra_InferiorcerebellarpeduncleR", "BRIA - Vintra rInf.Cer.Ped.",plot_beta$Outcome)
plot_beta$Outcome = gsub("smt_mc_intra_InferiorcerebellarpeduncleR", "SMTmc - intra rInf.Cer.Ped.",plot_beta$Outcome)
plot_beta$Outcome = gsub("smt_fa_InferiorcerebellarpeduncleR", "SMT - FA rInf.Cer.Ped.",plot_beta$Outcome)
plot_beta$Outcome = gsub("smt_fa_InferiorcerebellarpeduncleL", "SMT - FA lInf.Cer.Ped.",plot_beta$Outcome)
plot_beta$Outcome = gsub("smt_mc_diff_InferiorcerebellarpeduncleR", "SMTmc - diff rInf.Cer.Ped.",plot_beta$Outcome)
plot_beta$Outcome = gsub("micro_Ax_InferiorcerebellarpeduncleR", "BRIA - microAX rInf.Cer.Ped.",plot_beta$Outcome)
plot_beta$Outcome = gsub("AD_PosteriorthalamicradiationL", "DTI - AD lPost.Thal.Rad.",plot_beta$Outcome)
plot_beta$Outcome = gsub("v_extra_InferiorcerebellarpeduncleR", "BRIA - Vextra rInf.Cer.Ped.",plot_beta$Outcome)
plot_beta$Outcome = gsub("micro_Ax_CSTL", "BRIA - microAX lCST",plot_beta$Outcome)
plot_beta$Outcome = gsub("AD_CSTL", "DTI - AD lCST",plot_beta$Outcome)
plot_beta$Outcome = gsub("v_extra_AnteriorcoronaradiataL", "BRIA - Vextra lAnt.Cor.Rad.",plot_beta$Outcome)
plot_beta$Outcome = gsub("micro_ADC_UncinatefasciculusL", "BRIA - miscroADC lUni.Fasc.",plot_beta$Outcome)
plot_beta$Outcome = gsub("awf_CerebralpeduncleL", "WMTI - AWF lCereb.Ped.",plot_beta$Outcome)
plot_beta$Outcome = gsub("v_csf_Middlecerebellarpeduncle", "BRIA - V CSF Mid.Cereb.Ped.",plot_beta$Outcome)
plot_beta$Outcome = gsub("micro_FA_InferiorcerebellarpeduncleR", "BRIA - microFA rInf.Cereb.Ped.",plot_beta$Outcome)
plot_beta$Outcome = gsub("v_csf_CorticospinaltractL", "BRIA - V CSF Corticospinal Tr.",plot_beta$Outcome)
plot_beta$Outcome = gsub("smt_mc_diff_AnteriorcoronaradiataL", "SMTmc - diff lAnt.Cor.Rad.",plot_beta$Outcome)
plot_beta$Outcome = gsub("smt_mc_extramd_AnteriorlimbofinternalcapsuleL", "SMTmc - extraMD lAnt.Limbinf.Ext.Caps.",plot_beta$Outcome)


plot_beta = rename(plot_beta, PGRS = variable)
region_plot = ggplot(data=plot_beta, aes(x=value, y=-log10(p), group = PGRS, label=Outcome)) +
  geom_point(aes(shape = PGRS, color = PGRS)) + 
  theme_minimal() +
  geom_vline(xintercept = c(-0.5, 0.5), color = "red", linetype = "dashed", cex = 1, alpha = 0.2) +
  geom_hline(yintercept=-log10(0.005), color = "red", linetype = "dashed", cex = 1, alpha = 0.2) +
  geom_text_repel(data = subset(plot_beta, value < -0.05 & p < 0.005),colour='black', nudge_x = -0.1, direction = "y", segment.size = 0.1, xlim = c(-0.1,-0.3))+
  geom_text_repel(data = subset(plot_beta, value > 0.05 & p < 0.005),colour='black', nudge_x = 0.1, direction = "y",segment.size = 0.1, xlim = c(0.1,0.3))+
  #scale_color_manual(values=c("#56B4E9","#999999","#D55E00"))+
  #scale_color_manual(values=c("blue","red", "black", "orange", "purple"))+
  #geom_vline(xintercept=c(-0.075, 0.075), col="black") +
  #geom_hline(yintercept=-log10(0.05), col="black") + 
  xlab("Association Strength")+ylab("Un-corrected -log10(p)")+
  xlim(-0.3,0.3)+ 
  scale_shape_manual(values = 1:nlevels(factor(plot_beta$PGRS)))+
  theme(legend.position="bottom")+scale_colour_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))

Rplot = ggarrange(plot, region_plot, ncol = 2, labels = c("a","b"))
ggsave("/tsd/p33/home/p33-maxk/export/PGRS_region_plot.pdf",Rplot, width = 18, height = 8)
ggsave("/cluster/projects/p33/users/maxk/UKB/export/LONG_WMM_change_PGRS.pdf",Rplot, width = 18, height = 8)

### STATS for ARoC-PGRS associations ####
# overall
plot_beta %>% summarise(M = mean(abs(value)), SD = sd(abs(value)), Md= median(abs(value)),MAD = mad(abs(value)))
# by disorder
plot_beta %>% group_by(PGRS) %>% summarise(M = mean(abs(value)), SD = sd(abs(value)), Md= median(abs(value)),MAD = mad(abs(value)))
# by strongest associations
## inferior cerebellar peduncle
### overall
plot_beta %>% filter(str_detect(Outcome, "Inferiorcerebellarpeduncle")) %>% summarise(M = mean(abs(value)), SD = sd(abs(value)), Md= median(abs(value)),MAD = mad(abs(value)), MAX = max(abs(value)))
### by disorder PGRS
plot_beta %>% filter(str_detect(Outcome, "Inferiorcerebellarpeduncle")) %>% group_by(PGRS) %>% summarise(M = mean(abs(value)), SD = sd(abs(value)), Md= median(abs(value)),MAD = mad(abs(value)), MAX = max(abs(value)))
## middle cerebellar peduncle
### overall
plot_beta %>% filter(str_detect(Outcome, "Middlecerebellarpeduncle")) %>% summarise(M = mean(abs(value)), SD = sd(abs(value)), Md= median(abs(value)),MAD = mad(abs(value)), MAX = max(abs(value)))
### by disorder PGRS
plot_beta %>% filter(str_detect(Outcome, "Middlecerebellarpeduncle")) %>% group_by(variable) %>% summarise(M = mean(abs(value)), SD = sd(abs(value)), Md= median(abs(value)),MAD = mad(abs(value)), MAX = max(abs(value)))
## superior cerebellar peduncle
### overall
plot_beta %>% filter(str_detect(Outcome, "Superiorcerebellarpeduncle")) %>% summarise(M = mean(abs(value)), SD = sd(abs(value)), Md= median(abs(value)),MAD = mad(abs(value)), MAX = max(abs(value)))
### by disorder PGRS
plot_beta %>% filter(str_detect(Outcome, "Superiorcerebellarpeduncle")) %>% group_by(variable) %>% summarise(M = mean(abs(value)), SD = sd(abs(value)), Md= median(abs(value)),MAD = mad(abs(value)), MAX = max(abs(value)))

# # Unicate fasciculus
# ## by disease PGRS
# plot_beta %>% filter(str_detect(Outcome, "Uncinatefasciculus")) %>% group_by(variable) %>% summarise(M = mean(abs(value)), SD = sd(abs(value)), Md= median(abs(value)),MAD = mad(abs(value)), MAX = max(abs(value)))
# ## overall
# plot_beta %>% filter(str_detect(Outcome, "Uncinatefasciculus")) %>% summarise(M = mean(abs(value)), SD = sd(abs(value)), Md= median(abs(value)),MAD = mad(abs(value)), MAX = max(abs(value)))



## fornix
### overall
plot_beta %>% filter(str_detect(Outcome, "ornix"))%>% filter(!str_detect(Outcome, "ornix_")) %>%  summarise(M = mean(abs(value)), SD = sd(abs(value)), Md= median(abs(value)),MAD = mad(abs(value)), MAX = max(abs(value)))
### by disorder
plot_beta %>% filter(str_detect(Outcome, "ornix"))%>% filter(!str_detect(Outcome, "ornix_")) %>% group_by(PGRS) %>% summarise(M = mean(abs(value)), SD = sd(abs(value)), Md= median(abs(value)),MAD = mad(abs(value)), MAX = max(abs(value)))

############################
# REGIONAL PGRS ASSOCIATIONS FOR EACH TIME POINT ####
# now, estimate the annual regional change and how it relates to PGRS
ROC1 = allT1[1:Ncol]
lm_list = list()
age = allT1$age
sex = allT1$sex
site = allT1$site_t3
outdat01 = data.frame(matrix(ncol = ncol(PGRS1), nrow = ncol(ROC1)))
ps01 = outdat
for (o in 1:ncol(PGRS1)){
  tmpPGRS = scale(PGRS1[[o]])
  for (i in 1:ncol(ROC1)){
    tmpROC = scale(ROC1[i], center = FALSE)
    tmp_df = data.frame(age, sex, site, tmpROC, tmpPGRS)
    tmp_model = lm(tmpROC ~ tmpPGRS + age*sex + (site), data = tmp_df)
    outdat01[i,o] = summary(tmp_model)$coefficients[2]
    ps01[i,o] = summary(tmp_model)$coefficients[2,4]
  }
}
ROC2 = allT2[1:Ncol]
lm_list = list()
age = allT2$age
sex = allT2$sex
site = allT2$site_t3
outdat02 = data.frame(matrix(ncol = ncol(PGRS1), nrow = ncol(ROC2)))
ps02 = outdat
for (o in 1:ncol(PGRS1)){
  tmpPGRS = scale(PGRS1[[o]])
  for (i in 1:ncol(ROC2)){
    tmpROC = scale(ROC2[i], center = FALSE)
    tmp_df = data.frame(age, sex, site, tmpROC, tmpPGRS)
    tmp_model = lm(tmpROC ~ tmpPGRS + age*sex + (site), data = tmp_df)
    outdat02[i,o] = summary(tmp_model)$coefficients[2]
    ps02[i,o] = summary(tmp_model)$coefficients[2,4]
  }
}

## Plot regional associations for each time point ####
colnames(outdat01)= colnames(PGRS1)
colnames(ps01)= colnames(PGRS1)
colnames(outdat02)= colnames(PGRS1)
colnames(ps02)= colnames(PGRS1)

# adjust p-values(?)
correctedP = ps01*(Ncol*ncol(PGRS1))
for (i in colnames(correctedP)){
  print(correctedP %>% filter(i < .05))
}
correctedP = ps02*(Ncol*ncol(PGRS1))
for (i in colnames(correctedP)){
  print(correctedP %>% filter(i < .05))
}
# this is presented as a question, as there are no p-values surviving the multiple comparison.

# However, we can look at the smallest p-values and largest associations in a volcano plot
## TP 1
plot_beta01 = melt(outdat01)
plot_p01 = melt(ps01)
plot_beta01$p = plot_p01$value
plot_beta01$Outcome = colnames(ROC)
## TP 2
plot_beta02 = melt(outdat02)
plot_p02 = melt(ps02)
plot_beta02$p = plot_p02$value
plot_beta02$Outcome = colnames(ROC)
# write the plot_beta data frame for potential later purposes
write.csv(plot_beta01, "/cluster/projects/p33/users/maxk/UKB/longitudinal/results/Longitudinal/regional_WMM_TP1_PGRS_associations.csv")
write.csv(plot_beta02, "/cluster/projects/p33/users/maxk/UKB/longitudinal/results/Longitudinal/regional_WMM_TP2_PGRS_associations.csv")

# now make the plot
plot_beta01 = rename(plot_beta01, PGRS = variable)
plot_beta02 = rename(plot_beta02, PGRS = variable)
# rename features
plot_beta01$Outcome = gsub("smt_fa_CINGR", "SMT - FA rCING",plot_beta01$Outcome)
plot_beta01$Outcome = gsub("v_csf_.PosteriorlimbofinternalcapsuleR", "BRIA - vCSF rPost.Limb.Int.Caps.",plot_beta01$Outcome)
plot_beta01$Outcome = gsub("micro_Ax_.PosteriorlimbofinternalcapsuleR", "WMTI - microAX rPost.Limb.Int.Caps.",plot_beta01$Outcome)
plot_beta01$Outcome = gsub("v_csf_UFR", "BRIA - vCSF UFR",plot_beta01$Outcome)
plot_beta01$Outcome = gsub("radEAD_UFR", "WMTI - radEAD UFR",plot_beta01$Outcome)
plot_beta01$Outcome = gsub("micro_ADC_SuperiorlongitudinalfasciculusR", "BRIA - microADC rSup.Long.Fasc.",plot_beta01$Outcome)
plot_beta01$Outcome = gsub("smt_mc_extramd_UFR", "SMTmc - extraMD rUF", plot_beta01$Outcome)
plot_beta01$Outcome = gsub("smt_long_SLFTL", "SMT - long lSLTF", plot_beta01$Outcome)
plot_beta01$Outcome = gsub("smt_mc_extramd_SLFTL", "SMTmc - extraMD lSLTF", plot_beta01$Outcome)
plot_beta01$Outcome = gsub("radEAD_SLFTL", "WMTI - radEAD lSLTF", plot_beta01$Outcome)
plot_beta01$Outcome = gsub("Dax_extra_SLFTL", "BRIA - DAX extra lSLTF", plot_beta01$Outcome)
plot_beta01$Outcome = gsub("Dax_intra_SLFTL", "BRIA - DAX intra lSLTF", plot_beta01$Outcome)
plot_beta01$Outcome = gsub("smt_md_SLFTL", "SMT - MD lSLTF", plot_beta01$Outcome)
plot_beta01$Outcome = gsub("micro_Ax_SLFTL", "BRIA - microAX lSLTF", plot_beta01$Outcome)
plot_beta01$Outcome = gsub("micro_ADC_SLFTL", "BRIA - microADC lSLTF", plot_beta01$Outcome)
plot_beta01$Outcome = gsub("smt_mc_diff_ExternalcapsuleR", "SMTmc - diff rExt.Caps.", plot_beta01$Outcome)
plot_beta01$Outcome = gsub("micro_Ax_AnteriorcoronaradiataL", "BRIA - microAX lAnt.Cor.Rad.", plot_beta01$Outcome)
plot_beta01$Outcome = gsub("MD_SLFTL", "DTI - MD lSLTF", plot_beta01$Outcome)
plot_beta01$Outcome = gsub("RD_UFR", "DTI - RD UFR", plot_beta01$Outcome)
plot_beta01$Outcome = gsub("v_csf_SLFTL", "BRIA - V CSF lSLTF", plot_beta01$Outcome)
plot_beta01$Outcome = gsub("AD_CerebralpeduncleR", "DTI - AD rCereb.Ped.", plot_beta01$Outcome)
plot_beta01$Outcome = gsub("smt_mc_diff_SuperiorlongitudinalfasciculusL", "SMTmc - diff lSup.Long.Fasc.", plot_beta01$Outcome)
plot_beta01$Outcome = gsub("awf_UncinatefasciculusL", "WMTI - AWF lUF", plot_beta01$Outcome)

plot_beta01$PGRS = plot_beta01$variable


# plot
region_plot01 = ggplot(data=plot_beta01, aes(x=value, y=-log10(p), group = PGRS, label=Outcome)) +
  geom_point(aes(shape = PGRS, color = PGRS)) + 
  theme_minimal() +
  geom_vline(xintercept = c(-0.001, 0.001), color = "red", linetype = "dashed", cex = 1, alpha = 0.2) +
  geom_hline(yintercept=-log10(0.001), color = "red", linetype = "dashed", cex = 1, alpha = 0.2) +
  geom_text_repel(data = subset(plot_beta01, value < -0.001 & p < 0.001),colour='black', nudge_x = -0.025, direction = "y", segment.size = 0.1, xlim = c(-0.01,-0.023))+
  geom_text_repel(data = subset(plot_beta01, value > 0.001 & p < 0.001),colour='black', nudge_x = 0.025, direction = "y",segment.size = 0.1, xlim = c(0.01,0.023))+
  #scale_color_manual(values=c("#56B4E9","#999999","#D55E00"))+
  #scale_color_manual(values=c("blue","red", "black", "orange", "purple"))+
  #geom_vline(xintercept=c(-0.075, 0.075), col="black") +
  #geom_hline(yintercept=-log10(0.05), col="black") + 
  xlab("Association Strength")+ylab("Un-corrected -log10(p)")+
  xlim(-0.022,0.022)+ 
  scale_shape_manual(values = 1:nlevels(factor(plot_beta01$PGRS)))+
  theme(legend.position="bottom")+scale_colour_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))

# rename features
plot_beta02$Outcome = gsub("ak_FMAJ", "DKI - AK FMAJ",plot_beta02$Outcome)
plot_beta02$Outcome = gsub("v_extra_ExternalcapsuleR", "BRIA - Vextra rExt. Caps.",plot_beta02$Outcome)
plot_beta02$Outcome = gsub("v_extra_ExternalcapsuleL", "BRIA - Vextra lExt. Caps.",plot_beta02$Outcome)
plot_beta02$Outcome = gsub("smt_mc_diff_ExternalcapsuleR", "SMTmc - diff rExt. Caps.",plot_beta02$Outcome)
plot_beta02$Outcome = gsub("smt_mc_intra_ExternalcapsuleL", "SMTmc - intra lExt. Caps.",plot_beta02$Outcome)
plot_beta02$Outcome = gsub("v_extra_TapetumR", "BRIA - Vextra rTapetum",plot_beta02$Outcome)
plot_beta02$PGRS=plot_beta02$value

# then plot
region_plot02 = ggplot(data=plot_beta02, aes(x=value, y=-log10(p), group = PGRS, label=Outcome)) +
  geom_point(aes(shape = PGRS, color = PGRS)) + 
  theme_minimal() +
  geom_vline(xintercept = c(-0.001, 0.001), color = "red", linetype = "dashed", cex = 1, alpha = 0.2) +
  geom_hline(yintercept=-log10(0.001), color = "red", linetype = "dashed", cex = 1, alpha = 0.2) +
  geom_text_repel(data = subset(plot_beta02, value < -0.001 & p < 0.001),colour='black', nudge_x = -0.004, direction = "y", segment.size = 0.05, xlim = c(-0.005,-0.11))+
  geom_text_repel(data = subset(plot_beta02, value > 0.001 & p < 0.001),colour='black', nudge_x = 0.004, direction = "y",segment.size = 0.05, xlim = c(0.005,0.11))+
  xlab("Association Strength")+ylab("Un-corrected -log10(p)")+
  xlim(-0.012,0.012)+ 
  scale_shape_manual(values = 1:nlevels(factor(plot_beta02$PGRS)))+
  theme(legend.position="bottom")+scale_colour_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))
T1T2plot = ggarrange(region_plot01,region_plot02, ncol = 2, labels = c("a","b"))
ggsave("/tsd/p33/home/p33-maxk/export/LONG_WMM_T1T2_REGIONAL_PGRS.pdf",T1T2plot, width = 18, height = 10)
ggsave("/cluster/projects/p33/users/maxk/UKB/export/LONG_WMM_T1T2_REGIONAL_PGRS.pdf",T1T2plot, width = 18, height = 10)


##################################
## CROSS-SECTIONAL WMM-PGRS Associations ####
PGRS = PGRS_copy

# exclude outliers / impossible values
cross_sectional = cross_sectional[!cross_sectional$eid %in% outliers$x,]
# make a copy of the data frame
cross_sectional1 = cross_sectional
# exclude the longitudinal data from the cross-sectional sample
cross_sectional = cross_sectional[!cross_sectional$eid %in% T1$eid,]
# merge eids for PGRS and WM / diffusion metrics
cross_sectional = cross_sectional[cross_sectional$eid %in% PGRS$FID,]
PGRS = PGRS[PGRS$FID %in% cross_sectional$eid,]
#### GLOBAL ASSOCIATIONS ####
# prep PGRS and WM frames
# select mean/global averages only
cross_sectional = cross_sectional %>% dplyr::select(contains("mean"), age, sex, site, eid)
# unselect the excluded metrics: BRIA DRADextra and WMTI axEAD
cross_sectional = cross_sectional %>% dplyr::select(-contains("DRAD"))%>% dplyr::select(-contains("axEAD"))

# now, select only actual PGRS scores for data frame
PGRS1 = PGRS[3:ncol(PGRS)]

## Global Metrics - PGRS associations ####
# scale mean metrics
cross_sectional[1:26] = cross_sectional[1:26]%>%mutate_if(is.numeric,scale)

# get PGRS-WM associations at time point 1 BUT FOR THE CROSS-SECTIONAL VALIDATION SAMPLE
outdat1 = data.frame(matrix(ncol = ncol(PGRS1), nrow = (26)))
ps1 = outdat1
age = cross_sectional$age
sex = cross_sectional$sex
site = cross_sectional$site
for (o in 1:ncol(PGRS1)){
  tmpPGRS = scale(PGRS1[[o]])
  for (i in 1:26){
    tmpROC = cross_sectional[[i]]
    age = cross_sectional$age
    tmp_df = data.frame(age, sex, site, tmpROC, tmpPGRS)
    tmp_model = lm(tmpROC ~ tmpPGRS + age*sex + (site), data = tmp_df)
    outdat1[i,o] = summary(tmp_model)$coefficients[2]
    ps1[i,o] = summary(tmp_model)$coefficients[2,4]
  }
}
#  rename columns
## for beta values
colnames(outdat1) = colnames(PGRS1)
## add labels
outdat1$names = var_labels

## for p values
colnames(ps1) = colnames(PGRS1)
ps1$names = var_labels

# melt data frame into long format for plotting
plot_dat1 = melt(outdat1)
p_frame1 = melt(ps1)

# this is what happens when adjusting the p-values. (No significant associations left after correcting for multiple comparison)
min(p.adjust(p_frame1$value, method = "fdr"))
# we could also select findings based on p-values before correction... (this is currently used in the plots)
p_frame1 %>% select(names, variable, value) %>% dplyr::filter(value < .05)

# use UNCORRECTED p-values for visualisation
plot_dat1$p = p_frame1$value

# to frame only significant associations define the colour of the frames
colors1 = c(ifelse(plot_dat1$p < .05,"black", "white"))

### Global Associations PLOTS ####
plot1gg = ggplot(plot_dat1, aes(variable, names,fill = value)) +
  geom_tile(lwd = 0.75, linetype = 1, color = colors1, width = 0.9, height = 0.9) +  
  geom_text(aes(label=round(value,2)),size=4)+
  scale_fill_gradient(low = "#0072B2", high = "#F0E442") +
  ylab("") + xlab("") + theme_bw() + guides(fill=guide_legend(title="Association Strength"))
# perhaps combine plots?
plot3gg = ggarrange(plot1,plot2,plot1gg, labels = c("a","b", "c"),ncol = 3, common.legend = TRUE)
ggsave("/tsd/p33/home/p33-maxk/export/LONG_crossWMM_PGRS.pdf",plot3gg, width = 17, height = 12)
ggsave("/cluster/projects/p33/users/maxk/UKB/export/LONG_crossWMM_PGRS.pdf",plot3gg, width = 17, height = 12)


########
#### REGIONAL ASSOCIATIONS ####
# exclude the longitudinal data from the cross-sectional sample
cross_sectional1 = cross_sectional1[!cross_sectional1$eid %in% T1$eid,]
# merge eids for PGRS and WM / diffusion metrics
cross_sectional1 = cross_sectional1[cross_sectional1$eid %in% PGRS$FID,]
PGRS = PGRS[PGRS$FID %in% cross_sectional1$eid,]
# now, select only actual PGRS scores for data frame
PGRS1 = PGRS[3:ncol(PGRS)]
# make some empty lists
ROCcross = cross_sectional1[1:Ncol]
lm_list = list()
age = cross_sectional1$age
sex = cross_sectional1$sex
site = cross_sectional1$site
outdat03 = data.frame(matrix(ncol = ncol(PGRS1), nrow = ncol(ROCcross)))
ps03 = outdat
# then start looping
for (o in 1:ncol(PGRS1)){
  tmpPGRS = scale(PGRS1[[o]])
  for (i in 1:ncol(ROCcross)){
    tmpROC = scale(ROCcross[i], center = FALSE)
    tmp_df = data.frame(age, sex, site, tmpROC, tmpPGRS)
    tmp_model = lm(tmpROC ~ tmpPGRS + age*sex + site, data = tmp_df)
    outdat03[i,o] = summary(tmp_model)$coefficients[2]
    ps03[i,o] = summary(tmp_model)$coefficients[2,4]
  }
}

## Plot regional associations for each time point ####
colnames(outdat03)= colnames(PGRS1)
colnames(ps03)= colnames(PGRS1)

# adjust p-values(?)
correctedP = ps03*(Ncol*ncol(PGRS1))
for (i in colnames(correctedP)){
  print(correctedP %>% filter(i < .05))
}
# this is presented as a question, as there are no p-values surviving the multiple comparison.

# However, we can look at the smallest p-values and largest associations in a volcano plot
plot_beta03 = melt(outdat03)
plot_p03 = melt(ps03)
plot_beta03$p = plot_p03$value
plot_beta03$Outcome = colnames(ROCcross)

# write the plot_beta data frame for potential later purposes
write.csv(plot_beta03, "/cluster/projects/p33/users/maxk/UKB/longitudinal/results/Longitudinal/regional_WMM_CROSS_PGRS_associations.csv")

# now make the plot
plot_beta03 = rename(plot_beta03, PGRS = variable)

# rename features
plot_beta03$Outcome = gsub("DKI.DKI.MK_Pontine", "DKI - MK Pontine",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("DKI.DKI.AK_CorticospinaltractR", "DKI - AK rCor.Sp.Tract",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("DKI.DKI.MK_CorticospinaltractR", "DKI - MK rCor.Sp.Tract",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("DKI.DKI.AK_Pontine", "DKI - AK Pontine",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("DKI.DKI.MK_CorticospinaltractR", "DKI - MK rCor.Sp.Tract",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("DKI.DKI.MK_CerebralpeduncleR", "DKI - MK rCereb.Peduncle",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("DKI.DKI.MK_CorticospinaltractL", "DKI - MK rCor.Sp.Tract",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("SMT_mc.SMT_mc.smt_mc_intra_CorticospinaltractR", "SMTmc - intra rCor.Sp.Tract",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("BRIA.BRIA.v_intra_CorticospinaltractR", "BRIA - Vintra rCor.Sp.Tract",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("BRIA.BRIA.v_intra_Pontine", "BRIA - Vintra Pontine",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("BRIA.BRIA.v_extra_CingulumhippocampusR", "BRIA - Vextra rCing.Hippocamp.",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("SMT_mc.SMT_mc.smt_mc_intra_Pontine", "SMTmc - extratrans Pontine",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("SMT_mc.SMT_mc.smt_mc_extratrans_CorticospinaltractR", "SMTmc - extratrans rCor.Sp.Tract",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("SMT_mc.SMT_mc.smt_mc_extratrans_CorticospinaltractL", "SMTmc - extratrans lCorticosp.Tract",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("SMT_mc.SMT_mc.smt_mc_extratrans_Pontine", "SMTmc - extratrans Pontine",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("DTI.DTI.MD_CorticospinaltractR", "DTI - MD rCor.Sp.Tract",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("DTI.DTI.RD_CorticospinaltractR", "DTI - RD rCor.Sp.Tract",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("BRIA.BRIA.micro_FA_CorticospinaltractR", "BRIA - Vextra rCing.Hippocamp.",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("SMT_mc.SMT_mc.smt_mc_extramd_CorticospinaltractR", "SMTmc - extraMD rCor.Sp.Tract",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("DTI.DTI.FA_UncinatefasciculusR", "DTI - FA rUnci.Fasc.",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("SMT.SMT.smt_trans_CorticospinaltractR", "SMT - trans rCor.Sp.Tract",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("DTI.DTI.FA_CINGL", "DTI - FA lCING",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("DTI.DTI.FA_CingulumhippocampusL", "DTI - FA lCing.Hippocamp.",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("BRIA.BRIA.micro_Rd_CorticospinaltractR", "BRIA - microRD rCor.Sp.Tract",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("BRIA.BRIA.v_extra_UncinatefasciculusL", "BRIA - vextra lUnci.Fasc.",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("SMT_mc.SMT_mc.smt_mc_extratrans_MediallemniscusR", "SMTmc - extratrans rMed.Lemn.",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("DTI.DTI.FA_SuperiorfrontooccipitalfasciculusR", "DTI - FA lSup.Frontoocc.Fasc.",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("BRIA.BRIA.micro_Rd_UncinatefasciculusL", "BRIA - microRD lUF",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("BRIA.BRIA.v_extra_UncinatefasciculusR", "BRIA - vextra rUnci.Fasc.",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("SMT_mc.SMT_mc.smt_mc_extratrans_ILFR", "SMTmc - extra trans rILF",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("DTI.DTI.RD_.PosteriorlimbofinternalcapsuleL", "DTI - RD lPost.Limb Int.Caps",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("DTI.DTI.RD_.PosteriorlimbofinternalcapsuleR", "DTI - RD rPost.Limb Int.Caps",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("WMTI.WMTI.axEAD_FMIN", "WMTI - axEAD FMIN",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("DKI.DKI.AK_SLTFR", "DKI - AK rSLTF",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("BRIA.BRIA.micro_ADC_AnteriorlimbofinternalcapsuleL", "BRIA - microADC lAnt.Limb",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("BRIA.BRIA.micro_ADC_AnteriorlimbofinternalcapsuleR", "BRIA - microADC rAnt.Limb",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("DTI.DTI.RD_CSTR", "DTI - RD rCST",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("DTI.DTI.RD_AnteriorlimbofinternalcapsuleR", "DTI - RD rAnt.Limb Ext.Caps",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("DTI.DTI.RD_AnteriorlimbofinternalcapsuleL", "DTI - RD lAnt.Limb Ext.Caps",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("DTI.DTI.RD_SuperiorfrontooccipitalfasciculusR", "DTI - RD rSup.Frontoocc. Fasc.",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("DTI.DTI.RD_SuperiorfrontooccipitalfasciculusL", "DTI - RD lSup.Frontoocc. Fasc.",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("DTI.DTI.RD_SLTFR", "DTI - RD rSLTF",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("DKI.DKI.AK_CINGR", "DKI - AK rCING",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("WMTI.WMTI.axEAD_AnteriorcoronaradiataL", "WMTI - axEAD lAnt.Cor.Rad",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("WMTI.WMTI.axEAD_AnteriorcoronaradiataR", "WMTI - axEAD rAnt.Cor.Rad",plot_beta03$Outcome)



# plot
cross_plot = ggplot(data=plot_beta03, aes(x=value, y=-log10(p), group = PGRS, label=Outcome)) +
  geom_point(aes(shape = PGRS, color = PGRS)) + 
  theme_minimal() +
  geom_vline(xintercept = c(-0.001, 0.001), color = "red", linetype = "dashed", cex = 1, alpha = 0.2) +
  geom_hline(yintercept=-log10(0.001), color = "red", linetype = "dashed", cex = 1, alpha = 0.2) +
  geom_text_repel(data = subset(plot_beta03, value < -0.001 & p < 0.001),colour='black', nudge_x = -0.05, direction = "y", segment.size = 0.1, xlim = c(-0.01,-0.02))+
  geom_text_repel(data = subset(plot_beta03, value > 0.001 & p < 0.001),colour='black', nudge_x = 0.05, direction = "y",segment.size = 0.1, xlim = c(0.01,0.02))+
  xlab("Association Strength")+ylab("Un-corrected -log10(p)")+
  xlim(-0.02,0.02)+ 
  scale_shape_manual(values = 1:nlevels(factor(plot_beta03$PGRS)))+
  theme(legend.position="bottom")+scale_colour_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"))



regional_plots = ggarrange(region_plot01,region_plot02,cross_plot, ncol = 3, labels = c("a","b", "c"))
ggsave("/tsd/p33/home/p33-maxk/export/LONG_WMM_ALL_REGIONAL_PGRS.pdf",regional_plots, width = 19, height = 8)
ggsave("/cluster/projects/p33/users/maxk/UKB/export/LONG_WMM_ALL_REGIONAL_PGRS.pdf",regional_plots, width = 18, height = 8)
# cross_only
ggsave("/tsd/p33/home/p33-maxk/export/LONG_WMM_CROSS_REGIONAL_PGRS.pdf",cross_plot, width = 9, height = 8)
ggsave("/cluster/projects/p33/users/maxk/UKB/export/LONG_WMM_CROSS_REGIONAL_PGRS.pdf",cross_plot, width = 9, height = 8)

### COMBINE ALL GLOBAL AND REGIONAL CROSS SECTIONAL PLOTS ####
# require(gridExtra)
all_tile = ggarrange(plot1,plot2,plot1gg,region_plot01,region_plot02,cross_plot, ncol = 3, nrow = 1, labels = c("a","b", "c"),common.legend = T)
all_volcano = ggarrange(region_plot01,region_plot02,cross_plot, ncol = 3, nrow = 1, labels = c("d", "e","f"),common.legend = T)
all_plot_combi = ggarrange((plot1),plot2,(plot1gg),region_plot01,region_plot02,cross_plot, ncol = 3, nrow = 2, labels = c("a","b", "c", "d", "e","f"),common.legend = T)
ggsave("/tsd/p33/home/p33-maxk/export/LONG_ALL_crossWMM_PGRS.pdf",all_plot_combi, width = 27, height = 20)
ggsave("/cluster/projects/p33/users/maxk/UKB/export/LONG_ALL_crossWMM_PGRS.pdf",all_plot_combi, width = 27, height = 20)
ggsave("/tsd/p33/home/p33-maxk/export/LONG_TILE_crossWMM_PGRS.pdf",all_tile, width = 27, height = 10)
ggsave("/cluster/projects/p33/users/maxk/UKB/export/LONG_TILE_crossWMM_PGRS.pdf",all_tile, width = 27, height = 10)
ggsave("/tsd/p33/home/p33-maxk/export/LONG_VOLCANO_crossWMM_PGRS.pdf",all_volcano, width = 27, height = 10)
ggsave("/cluster/projects/p33/users/maxk/UKB/export/LONG_VOLCANO_crossWMM_PGRS.pdf",all_volcano, width = 27, height = 10)



### PLOTTING PGRS ASSOCIATIONS WITH FORNIX MICROSTRUCTURE ####
# make 4 plots for WMM-PGRS associations:
## 1) annual rate of change
## 2) time point 1
## 3) time point 2
## 4) time point 1 validation

### STATS for regional cross-sectional ARoC-PGRS associations ####
# overall
plot_beta01 %>% summarise(M = mean(abs(value)), SD = sd(abs(value)), Md= median(abs(value)),MAD = mad(abs(value)))
plot_beta02 %>% summarise(M = mean(abs(value)), SD = sd(abs(value)), Md= median(abs(value)),MAD = mad(abs(value)))
plot_beta03 %>% summarise(M = mean(abs(value)), SD = sd(abs(value)), Md= median(abs(value)),MAD = mad(abs(value)))

#### FORNIX ####
plot_beta$Outcome = colnames(ROC)
plot_beta01$Outcome = colnames(ROC)
plot_beta02$Outcome = colnames(ROC)
plot_beta03$Outcome =colnames(ROCcross)
plot_beta03$Outcome = gsub("WMTI.WMTI.","",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("BRIA.BRIA.","",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("SMT.SMT.","",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("DKI.DKI.","",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("DTI.DTI.","",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("SMT_mc.SMT_mc.","",plot_beta03$Outcome)
fornix1$Outcome
fornixlist = list(plot_beta, plot_beta01, plot_beta02, plot_beta03)
fornixplots = list()
for (i in 1:length(fornixlist)){
  fornix1 = (dplyr::filter(fornixlist[[i]], !grepl("_Fornix_", Outcome)))
  fornix1 = (dplyr::filter(fornix1, grepl("Fornix", Outcome)))
  fornix1$Outcome = gsub("_Fornix","",fornix1$Outcome)
  fornix1$Outcome = var_labels
  colors = c(ifelse(fornix1$p < .05,"black", "white"))
  fornixplots[[i]] = ggplot(fornix1, aes(PGRS, Outcome,fill = value)) +
    geom_tile(lwd = 0.75, linetype = 1, color = colors, width = 0.9, height = 0.9) +  
    geom_text(aes(label=round(value,2)),size=4)+
    scale_fill_gradient(low = "#0072B2", high = "#F0E442") +
    ylab("") + xlab("") + theme_bw() + guides(fill=guide_legend(title="Association Strength"))
}
figure = ggarrange(fornixplots[[1]], fornixplots[[2]] , fornixplots[[3]], fornixplots[[4]], ncol = 4, common.legend = T, labels = c("a", "b", "c", "d"))

# possible to remove text with jrremove("y.text")
# perhaps change font size:
# font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top"))
# 
# test_fig = cowplot::plot_grid(fornixplots[[1]], 
#                    fornixplots[[2]] + 
#                      theme(axis.text.y = element_blank(),
#                            axis.ticks.y = element_blank(),
#                            axis.title.y = element_blank(),
#                            plot.margin = margin(r = 1) ),
#                    fornixplots[[3]] + 
#                      theme(axis.text.y = element_blank(),
#                            axis.ticks.y = element_blank(),
#                            axis.title.y = element_blank(),
#                            plot.margin = margin(r = 1) ),
#                    fornixplots[[4]] + 
#                      theme(axis.text.y = element_blank(),
#                            axis.ticks.y = element_blank(),
#                            axis.title.y = element_blank(),
#                            plot.margin = margin(r = 1) ), 
#                    nrow = 1,
#                    labels = "auto",
#                    align = "v")

# save the fornix figure with all panels for appendix
ggsave("/cluster/projects/p33/users/maxk/UKB/export/LONG_ALL_Fornix.pdf",figure, width = 30, height = 14)
ggsave("/tsd/p33/home/p33-maxk/export/LONG_ALL_FORNIX.pdf",figure, width = 30, height = 14)


#### CEREBRAL PREDUNCLE ####
# make plots for Cerebral Peduncle annual regional WMM change and cross-sectional associations
pedplots = list()
for (i in 1:length(fornixlist)){
  fornix1 = (dplyr::filter(fornixlist[[i]], grepl('Middlecerebellarpeduncle', Outcome)))
  fornix1$Outcome = gsub("_Middlecerebellarpeduncle","",fornix1$Outcome)
  fornix1$Outcome = var_labels
  colors = c(ifelse(fornix1$p < .05,"black", "white"))
  pedplots[[i]] = ggplot(fornix1, aes(PGRS, Outcome,fill = value)) +
    geom_tile(lwd = 0.75, linetype = 1, color = colors, width = 0.9, height = 0.9) +  
    geom_text(aes(label=round(value,3)),size=4)+
    scale_fill_gradient(low = "#0072B2", high = "#F0E442") +
    ylab("") + xlab("") + theme_bw() + guides(fill=guide_legend(title="Association Strength"))
}
pedfig = ggarrange(pedplots[[1]], pedplots[[2]] , pedplots[[3]], pedplots[[4]], ncol = 4, common.legend = T, labels = c("a", "b", "c", "d"))
ggsave("/tsd/p33/home/p33-maxk/export/LONG_regional_cross_sec_Peduncle.pdf",pedfig, width = 30, height = 14)
ggsave("/cluster/projects/p33/users/maxk/UKB/export/LONG_regional_cross_sec_Peduncle.pdf",pedfig, width = 30, height = 14)

# make a special plot including a panel for peduncle in addition to mean values and all regional values
plot99 = plot+ theme(legend.position="bottom")
plot100 = pedplots[[1]] + theme(legend.position="bottom")
newRplot = ggarrange(plot99, region_plot,plot100, ncol = 3, labels = c("a","b", "c"))
ggsave("/tsd/p33/home/p33-maxk/export/LONG_regional_PEDUNCLE.pdf",newRplot, width = 30, height = 14)
ggsave("/cluster/projects/p33/users/maxk/UKB/export/LONG_regional_PEDUNCLE.pdf",newRplot, width = 30, height = 14)


# make a plot for Cerebellar Peduncle annual regional WMM change
peddat = dplyr::filter(plot_beta03, grepl('eduncle', Outcome))
pedcolors = c(ifelse(peddat$p < .05,"black", "white"))
pedplot = ggplot(peddat, aes(PGRS, Outcome,fill = value)) +
  geom_tile(lwd = 0.75, linetype = 1, color = pedcolors, width = 0.9, height = 0.9) +  
  geom_text(aes(label=round(value,3)),size=4)+
  scale_fill_gradient(low = "#0072B2", high = "#F0E442") +
  ylab("") + xlab("") + theme_bw() + guides(fill=guide_legend(title="Association Strength"))
ggsave("/tsd/p33/home/p33-maxk/export/LONG_regional_Peduncle.pdf",pedplot, width = 13, height = 30)
ggsave("/cluster/projects/p33/users/maxk/UKB/export/LONG_regional_Peduncle.pdf",pedplot, width = 13, height = 30)