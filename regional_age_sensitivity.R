# HEMI SCRIPT!
# This script estimates AGE SENSITIVITY of ABSOLUTE brain asymmetries
# 30.10.2023 Max Korbmacher (max.korbmacher@gmail.com)
# 
# # PREP ####
# # load packages
# if (!require("pacman")) install.packages("pacman")
# pacman::p_load(dplyr, marginaleffects, ggpubr, ggplot2, Rmpfr, mgcv)
# # load data
# T1 = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/both/T1.csv")
# dMRI = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/both/dMRI.csv")
# multimodal = read.csv("/cluster/projects/p33/users/maxk/UKB/hemispheric_age/both/combi.csv")
# demo = read.delim("/cluster/projects/p33/users/maxk/UKB/data/50k_demogShort_020523.txt")
# 
# # split left from right
# LT1 = T1 %>% dplyr::select(eid, starts_with(c("lh", "Left")))
# RT1 = T1 %>% dplyr::select(eid, starts_with(c("rh","Right")))
# LdMRI = dMRI %>% dplyr::select(eid, ends_with(c("L")))
# RdMRI = dMRI %>% dplyr::select(eid, ends_with(c("R")))
# L = merge(LT1, LdMRI, by = "eid")
# R = merge(RT1,RdMRI, by = "eid")
# 
# # add demo to L for later regression
# L = merge(L,demo,by="eid")
# LdMRI = merge(LdMRI, demo, by = "eid")
# LT1 = merge(LT1,demo,by="eid")
# 
# ## for multimodal data
# L = L %>% dplyr::select(-eid) 
# R = R %>% dplyr::select(-eid) 
# 
# ## dMRI only
# LdMRI = LdMRI %>% dplyr::select(-eid) 
# RdMRI = RdMRI %>% dplyr::select(-eid) 
# ## T1 only
# LT1 = LT1 %>% dplyr::select(-eid) 
# RT1 = RT1 %>% dplyr::select(-eid) 
# 
# # define laterality index/asymmetry function (ABSOLUTE VALUES)
# LI = function(L,R){
#   abs((L-R)/(L+R))
# }
# # now, produce mean-centered (to zero), scaled laterality indexed regional metrics
# LIdat = data.frame(matrix(ncol = ncol(R),nrow = nrow(R)))
# LIdMRI = data.frame(matrix(ncol = ncol(RdMRI),nrow = nrow(RdMRI)))
# LIT1 = data.frame(matrix(ncol = ncol(RT1),nrow = nrow(RT1)))
# for (i in 1:ncol(R)){
#   LIdat[i] = as.numeric(scale(LI(L[[i]], R[[i]])))
# }
# for (i in 1:ncol(RdMRI)){
#   LIdMRI[i] = as.numeric(scale(LI(LdMRI[i], RdMRI[i])))
# }
# for (i in 1:ncol(RT1)){
#   LIT1[i] = as.numeric(scale(LI(LT1[i], RT1[i])))
# }
# names(LIdat) = names(R)
# names(LIdMRI) = names(RdMRI)
# names(LIT1) = names(RT1)
# 
# # remove mean average remainders
# remainders = function(data){
#   data %>% dplyr::select(!contains("mean"), !contains("total"))
# }
# LIdat = remainders(LIdat)
# LIdMRI = remainders(LIdMRI)
# LIT1 = remainders(LIT1)
# #
# # add sex, age, site
# LIdat$sex = L$Sex
# LIdat$age = L$Age
# LIdat$site = L$Scanner
# LIdMRI$sex = LdMRI$Sex
# LIdMRI$age = LdMRI$Age
# LIdMRI$site = LdMRI$Scanner
# LIT1$sex = LT1$Sex
# LIT1$age = LT1$Age
# LIT1$site = LT1$Scanner
# #
# # scale age for easier plotting of coefficients
# LIdat$age_scaled = scale(LIdat$age)
# LIdMRI$age_scaled = scale(LIdMRI$age)
# LIT1$age_scaled = scale(LIT1$age)
# # make list of LI data frames
# LIdfs = list(LIdat, LIdMRI, LIT1)
# #
# # TEST LINEAR EFFECTS ####
# # testing age associations
# # (due to convergence problems, we use simple linear and generalized additive models, instead of site as random effect)
# regional.lm = list()
# # we use Rmpfr to estimate accurate p-vals
# .N <- function(.) mpfr(., precBits = 200)
# # now, loop over list of data frames: again, using linear models!
# for (df in 1:length(LIdfs)){
#   coeff = c()
#   tvals = c()
#   p.val = c()
#   p.adj = c()
#   dat = LIdfs[[df]]
#   metrics = dat %>% dplyr::select(-site, -sex, -age, -age_scaled) %>% names()
#   for (i in metrics){
#     f = formula(paste(i,"~age_scaled+sex+site"))
#     tmpmodel = lm(f, data = LIdat)
#     coeff[i] = summary(tmpmodel)$coefficients[2]
#     tvals[i] = summary(tmpmodel)$coefficients[2,3]
#     print(paste("Done with estimating age-associations for", i))
#   }
#   p.val = 2 * pnorm(-abs(.N(tvals)))
#   p.adj = mpfr(p.val*length(metrics),100)
#   regional.lm[[df]] = data.frame(metrics,coeff,formatMpfr(p.val), formatMpfr(p.adj))
#   print(paste("Done with estimating age-associations in data frame number", df))
# }
# # save the three data frames
# # set work dir
# write.csv(regional.lm[[1]], "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_LI_LINEAR_age_associations_multimodal.csv")
# write.csv(regional.lm[[2]], "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_LI_LINEAR_age_associations_dMRI.csv")
# write.csv(regional.lm[[3]], "/cluster/projects/p33/users/maxk/UKB/export/Hemi_NEW_LI_LINEAR_age_associations_T1.csv")
