# clean demo data

# load packages
library(dplyr)

# set working directory to the data location
setwd("/home/max/Documents/Projects/Diffusion/UKBlong/data_export/")

# load data
demo = read.csv("demo.csv")
MRI = read.csv("T1_noDRAD.csv")

# select variables for the appropiate time points
demo = demo %>% dplyr::select("eid", contains(".2."), contains(".3."))

# drop the general vascular diagnosis variable, as there is only one time point available.
demo2 = demo %>% dplyr::select(!contains("6150")) %>% dplyr::select(!contains("6332"))

# name columns in a more understandable way
colnames(demo2) = c("eid","inc_pair_matches_r1_TP1","inc_pair_matches_r2_TP1","inc_pair_matches_r3_TP1","max_digits_rem_TP1", "fluid_intel_TP1", "prosp_mem_TP1","nb_meds", "health_selfr_TP1", 
  "Matrix_RT1_TP1", "Matrix_RT2_TP1","Matrix_RT3_TP1","Matrix_RT4_TP1", "Matrix_RT5_TP1","Matrix_RT6_TP1","Matrix_RT7_TP1", "Matrix_RT8_TP1","Matrix_RT9_TP1",
  "Matrix_RT10_TP1", "Matrix_RT11_TP1","Matrix_RT12_TP1","Matrix_RT13_TP1", "Matrix_RT14_TP1","Matrix_RT15_TP1", "cor_Matrix_puzzles_TP1","viewed_Matrix_puzzles_TP1",
  "cor_Tower_puzzles_TP1", "cor_SymDig_matches_TP1", "inc_pair_matches_R1_TP2", "inc_pair_matches_R2_TP2", "inc_pair_matches_R3_TP2","max_digits_rem_TP2",  "fluid_intel_TP2","prosp_mem_TP2","nb_meds", "health_selfr_TP2", 
  "Matrix_RT1_TP2", "Matrix_RT2_TP2","Matrix_RT3_TP2","Matrix_RT4_TP2", "Matrix_RT5_TP2","Matrix_RT6_TP2","Matrix_RT7_TP2", "Matrix_RT8_TP2","Matrix_RT9_TP2",
  "Matrix_RT10_TP2", "Matrix_RT11_TP2","Matrix_RT12_TP2","Matrix_RT13_TP2", "Matrix_RT14_TP2","Matrix_RT15_TP2", "cor_Matrix_puzzles_TP2","viewed_Matrix_puzzles_TP2",
  "cor_Tower_puzzles_TP2", "cor_SymDig_matches_TP2")

# match eids with the ones for which we have MRI data 
demo2 = demo2[demo2$eid %in% MRI$eid,]
# make time-point specific data frames
TP1 = demo2 %>% select(eid, ends_with("TP1"))
TP2 = demo2 %>% select(eid, ends_with("TP2"))
# remove suffix from dfs
TP1 = TP1 %>% 
  rename_at(.vars = vars(ends_with("_TP1")),
            .funs = funs(sub("_TP1$", "", .)))
colnames(TP2) = colnames(TP1)
# write the final file with the cognitive scores
write.csv(TP1, "TP1_cog.csv")
write.csv(TP2, "TP2_cog.csv")
