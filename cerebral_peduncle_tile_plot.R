# remake of the cerebral peduncle plot

# libs&data
library(ggplot2)
library(ggpubr)
plot_beta = read.csv("/cluster/projects/p33/users/maxk/UKB/longitudinal/results/Longitudinal/regional_WMM_ARoC_PGRS_associations.csv")
plot_beta01 = read.csv("/cluster/projects/p33/users/maxk/UKB/longitudinal/results/Longitudinal/regional_WMM_TP1_PGRS_associations.csv")
plot_beta02 = read.csv("/cluster/projects/p33/users/maxk/UKB/longitudinal/results/Longitudinal/regional_WMM_TP2_PGRS_associations.csv")
plot_beta03 = read.csv("/cluster/projects/p33/users/maxk/UKB/longitudinal/results/Longitudinal/regional_WMM_CROSS_PGRS_associations.csv")

# prep (rename badbly labelled data)
plot_beta03$Outcome = gsub("WMTI.WMTI.","",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("BRIA.BRIA.","",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("SMT.SMT.","",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("DKI.DKI.","",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("DTI.DTI.","",plot_beta03$Outcome)
plot_beta03$Outcome = gsub("SMT_mc.SMT_mc.","",plot_beta03$Outcome)

#### CEREBRAL PREDUNCLE ####
# make plots for Cerebra Peduncel annual regional WMM change and cross-sectional associations
pedunclelist = list(plot_beta, plot_beta01, plot_beta02, plot_beta03)
pedplots = list()
for (i in 1:length(pedunclelist)){
  peduncle1 = (dplyr::filter(pedunclelist[[i]], grepl('eduncle', Outcome)))
  # rename Cerebralpeduncle
  # peduncle1$Outcome=gsub("Cerebralpeduncle","Cer.Ped.",peduncle1$Outcome)
  colors = c(ifelse(peduncle1$p < .05,"black", "white"))
  pedplots[[i]] = ggplot(peduncle1, aes(variable, Outcome,fill = value)) +
    geom_tile(lwd = 0.75, linetype = 1, color = colors, width = 0.9, height = 0.9) +  
    geom_text(aes(label=round(value,3)),size=4)+
    scale_fill_gradient(low = "#0072B2", high = "#F0E442") +
    ylab("") + xlab("") + theme_bw() + guides(fill=guide_legend(title="Association Strength"))
}
pedfig = ggarrange(pedplots[[1]], pedplots[[2]] , pedplots[[3]], pedplots[[4]], ncol = 4, common.legend = T, labels = c("a", "b", "c", "d"))
ggsave("/tsd/p33/home/p33-maxk/export/LONG_regional_cross_sec_Peduncle.pdf",pedfig, width = 40, height = 30)
ggsave("/cluster/projects/p33/users/maxk/UKB/export/LONG_regional_cross_sec_Peduncle.pdf",pedfig, width = 40, height = 30)

# 
# # make a plot for Cerebra Peduncle annual regional WMM change
# peddat = dplyr::filter(plot_beta03, grepl('eduncle', Outcome))
# pedcolors = c(ifelse(peddat$p < .05,"black", "white"))
# pedplot = ggplot(peddat, aes(PGRS, Outcome,fill = value)) +
#   geom_tile(lwd = 0.75, linetype = 1, color = pedcolors, width = 0.9, height = 0.9) +  
#   geom_text(aes(label=round(value,3)),size=4)+
#   scale_fill_gradient(low = "#0072B2", high = "#F0E442") +
#   ylab("") + xlab("") + theme_bw() + guides(fill=guide_legend(title="Association Strength"))
# ggsave("/tsd/p33/home/p33-maxk/export/LONG_regional_Peduncle.pdf",pedplot, width = 13, height = 30)
# ggsave("/cluster/projects/p33/users/maxk/UKB/export/LONG_regional_Peduncle.pdf",pedplot, width = 13, height = 30)