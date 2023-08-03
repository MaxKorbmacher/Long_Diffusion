# test-retest reliability in UKB regional average white matter microstructure data
# Max Korbmacher, 25.07.2023
# load packages
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)

######## MEAN SKELETON VALUES #####

# set working directory where the data files are situated
setwd("/home/max/Documents/Projects/Diffusion/UKBlong/data_export/")
# load data (without DRAD extra outliers removed as there were many outliers)
T1 = read.csv("T1_noDRAD.csv")
T2 = read.csv("T2_noDRAD.csv")
# hence, we have to exclude DRAD extra from our analysis
# also, we start focusing on global average scores here
T1 = T1 %>% dplyr::select(contains("mean"), age, sex, site_t3, eid) %>% dplyr::select(!Drad_extra_Mean)
T2 = T2 %>% dplyr::select(contains("mean"), age, sex, site_t3, eid) %>% dplyr::select(!Drad_extra_Mean)
for (i in 1:27){
  T1[,i] = (T1[,i]-(mean(T1[,i])))/(sd(T1[,i]))
  T2[,i] = (T2[,i]-(mean(T2[,i])))/(sd(T2[,i]))
}
T1$TP = replicate(nrow(T1),0)
T2$TP = replicate(nrow(T1),1)
df = rbind(T1,T2)
age_delta = T2$age-T1$age
# run linear models correcting
model=list()
model2=list()
names1 = colnames(T1[1:27])
betas = c()
betas_uncor = c()
se = c()
se_uncor = c()
for (i in 1:27){
  f = formula(paste("T1$",names1[i],"~T2$",names1[i],"+age_delta+T1$sex+T1$age+T1$site_t3"))
  f2 = formula(paste("T1$",names1[i],"~T2$",names1[i]))
  model[[i]] = lm(f)
  model2[[i]] = lm(f2)
  betas[i] = summary(model[[i]])$coefficients[2,1]
  se[i] = summary(model[[i]])$coefficients[2,2]
  betas_uncor[i] = summary(model2[[i]])$coefficients[2,1]
  se_uncor[i] = summary(model2[[i]])$coefficients[2,2]
}
mean = data.frame(names1,betas,se,betas_uncor, se_uncor)
write.csv(mean,"/home/max/Documents/Projects/Diffusion/test-retest/mean_skeleton.csv")
mean(mean$betas)
sd(mean$betas)
#
########## ALL OTHER VALUES #######
# load data (without DRAD extra outliers removed as there were many outliers)
T1 = read.csv("T1_noDRAD.csv")
T2 = read.csv("T2_noDRAD.csv")
# hence, we have to exclude DRAD extra from our analysis
# also, we start focusing on global average scores here
T1d = T1 %>% dplyr::select(-c(Drad_extra_Mean,X,X.1,X.2,sex, age, site_t3,eid))
T2d = T2 %>% dplyr::select(-c(Drad_extra_Mean, X,eid)) %>%select(everything(), sex, age, site_t3)
# normalize data
T1d = T1d%>%mutate_if(is.numeric,scale)
T2d = T2d%>%mutate_if(is.numeric,scale)
names1 = colnames(T1d[1:1931])
for (i in length(names1)){
  f = formula(paste("T1d$",names1[i],"~T2d$",names1[i],"+age_delta+T2d$sex+T2d$age+T2d$site_t3"))
  f2 = formula(paste("T1d$",names1[i],"~T2d$",names1[i]))
  model[[i]] = lm(f)
  model2[[i]] = lm(f2)
  betas[i] = summary(model[[i]])$coefficients[2,1]
  se[i] = summary(model[[i]])$coefficients[2,2]
  betas_uncor[i] = summary(model2[[i]])$coefficients[2,1]
  se_uncor[i] = summary(model2[[i]])$coefficients[2,2]
}
all = data.frame(names1,betas,se,betas_uncor, se_uncor)
write.csv(all,"/home/max/Documents/Projects/Diffusion/test-retest/regional_skeleton.csv")
mean(all$betas)
sd(all$betas)

# plot skeleton averages
plot1 = ggplot(mean) +
  geom_bar( aes(x=names1, y=betas), stat="identity", fill="skyblue", alpha=0.7) +
  geom_pointrange( aes(x=names1, y=betas, ymin=betas-se, ymax=betas+se), colour="orange", alpha=0.9, size=0.7)+
  coord_flip() + theme_bw() + xlab("Diffusion Metric") + ylab("Adjusted test-retest correlation")
ggsave("/home/max/Documents/Projects/Diffusion/test-retest/test_retest_mean_plot.pdf", plot1, height = 5, width = 9)
# plot region and tract averages (rough density plot throwing all data together)
plot2 = ggplot(all, aes(x = betas)) +
  geom_density(color = 4,fill = 4, alpha = 0.25) + theme_bw() + 
  ylab("Density") + xlab("Adjusted test-retest correlation")
ggsave("/home/max/Documents/Projects/Diffusion/test-retest/test_retest_regions_plot1.pdf", plot2, height = 6, width = 6)
plot3 = ggplot(all, aes(x = betas_uncor)) +
  geom_density(color = 4,fill = 4, alpha = 0.25) + theme_bw() + 
  ylab("Density") + xlab("Un-adjusted test-retest correlation")
ggsave("/home/max/Documents/Projects/Diffusion/test-retest/Unadjusted_test_retest_regions_plot2.pdf", plot3, height = 6, width = 6)

# plot each metric
## create name list and empty list
names2 = c("v_intra","v_extra","v_csf",  "micro_Rd","micro_FA","micro_Ax","micro_ADC","Dax_intra","Dax_extra","rk", "ak", "mk",
           "FA_", "MD_", "RD_", "AD_", "smt_fa", "smt_long","smt_md", "smt_trans","smt_mc_intra","smt_mc_extramd","smt_mc_extratrans", 
           "smt_mc_diff","axEAD","awf","radEAD")
out = list()
## loop over name list for adjusted test-retest associations
for (i in 1:length(names2)){
  tmp.dat = all %>% filter(grepl(paste("^",names2[i], sep = ""), names1))
  tmp.plot = ggplot(tmp.dat) +
    geom_bar( aes(x=names1, y=betas), stat="identity", fill="skyblue", alpha=0.7) +
    geom_pointrange( aes(x=names1, y=betas, ymin=betas-se, ymax=betas+se), colour="orange", alpha=0.9, size=0.7)+
    coord_flip() + theme_bw() + xlab("Diffusion Metric") + ylab("Adjusted test-retest correlation")
  out[[i]] = tmp.plot
}
plot4 = do.call(ggarrange, c(out, widths = c(2, 1)))
ggsave("/home/max/Documents/Projects/Diffusion/test-retest/Adjusted_test_retest_regions_plot.pdf", plot4, height = 40, width = 30)
## loop over name list for UNadjusted test-retest associations
out = list()
for (i in 1:length(names2)){
  tmp.dat = all %>% filter(grepl(paste("^",names2[i], sep = ""), names1))
  tmp.plot = ggplot(tmp.dat) +
    geom_bar( aes(x=names1, y=betas_uncor), stat="identity", fill="skyblue", alpha=0.7) +
    geom_pointrange( aes(x=names1, y=betas_uncor, ymin=betas_uncor-se_uncor, ymax=betas_uncor+se_uncor), colour="orange", alpha=0.9, size=0.7)+
    coord_flip() + theme_bw() + xlab("Diffusion Metric") + ylab("Un-adjusted test-retest correlation")
  out[[i]] = tmp.plot
}
plot5 = do.call(ggarrange, c(out, widths = c(2, 1)))
ggsave("/home/max/Documents/Projects/Diffusion/test-retest/UnAdjusted_test_retest_regions_plot.pdf", plot5, height = 40, width = 30)

# NEXT STEP: CLOUD OF TEXT ALONG A GRADIENT
out1 = list()
for (i in 1:length(names2)){
  tmp.dat = all %>% filter(grepl(paste("^",names2[i], sep = ""), names1))
  tmp.plot1 <- ggplot(tmp.dat, aes(x = betas, y=se, label = names1)) +
    geom_point(aes(color = betas)) +
    geom_text_repel(aes(label = names1, color = betas)) + theme_bw() + ylab("Standard Error") + xlab("Adjusted test-retest correlation") + scale_color_gradient(low = "blue", high = "red")
  out1[[i]] = tmp.plot1
  }
# we want each page to include 2 plots 
plot6 = ggarrange(plotlist = out1, ncol = 2)
# then we can save them all separetely
for (i in 1:length(plot6)){
  ggsave(paste("/home/max/Documents/Projects/Diffusion/test-retest/Cloud_adjusted_test_retest_regions_plot",i,".pdf", sep = ""),plot6[[i]],height = 12, width = 17)
}
