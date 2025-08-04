---
#title: "Post-hoc ANOVA test Peritumoral area"
#author: "Mona Zimmermann"
#date: ""
#output: html_document

#### Email: m.l.m.zimmermann@amsterdamumc.nl
#### Purpose of script: Statistical analysis for the post hoc test of clinical 
#### relevance in activity changes for the peritumoral area

#### Status:
#Under development/reviewed/final version

#### Review History
#Reviewed by: 
#Date: 

#### Load libraries
library(tidyverse) 
library(ggpubr)
library(rstatix)
library(MASS)
library(dplyr)

# https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/
# mixed anova: https://www.datanovia.com/en/lessons/mixed-anova-in-r/


##### ----- Data preparation ----- #####
##### IMPORT DATA
##### ----- Data preparation ----- #####
file_path_peri_progression <- 'M:\\path\\to\20250404_dataframe_peri_baseline_FU_long_format_cov_prog.csv'
file_path_peri_epilepsy <- 'M:\\path\\to\20250404_dataframe_peri_baseline_FU_long_format_cov_epi.csv'
file_path_peri_mol <- 'M:\\path\\to\20250404_dataframe_peri_baseline_FU_long_format_cov_mol.csv'
file_path_peri_rth_grade <- 'M:\\path\\to\20250701_dataframe_peri_baseline_FU_long_format_cov_rth_grade.csv'

df_peri_prog <- read_csv(file_path_peri_progression)
df_peri_epi <- read_csv(file_path_peri_epilepsy)
df_peri_mol <- read_csv(file_path_peri_mol)
df_peri_rth_grade <- read_csv(file_path_peri_rth_grade)

#turn necessary categorical columns into factors to use them as covariates in the analysis (returns same results as anova already interpreted the categorical vars as factors)

factorize <-function(df, var_name){
  df$MM_factor <- factor(df$MM)
  df[[paste0(var_name, "_factor")]] <- factor(df[[var_name]])
  return(df)
}

df_peri_prog <- factorize(df_peri_prog, "progression")
df_peri_epi <- factorize(df_peri_epi, "epilepsy_aed")
df_peri_mol <- factorize(df_peri_mol, "IDH_1p19q")
df_peri_rth_grade <- factorize(df_peri_rth_grade, "Syntax_dummy_RTH_all")
df_peri_rth_grade <- factorize(df_peri_rth_grade, "graad")



#### ----- Summary statistics ----- ####

# --- BB_welch_z --- #
#MM alone 
df_peri_prog %>%
  group_by(MM) %>% 
  get_summary_stats(BB_welch_z, type = "mean_sd")

#MM and progression 
df_peri_prog %>%
  group_by(progression, MM) %>% 
  get_summary_stats(BB_welch_z, type = 'mean_sd')

##MM and epilepsy
df_peri_epi %>%
  group_by(epilepsy_aed, MM) %>% 
  get_summary_stats(BB_welch_z, type = 'mean_sd')

##MM and grade 
df_peri_rth_grade %>%
  group_by(graad, MM) %>% 
  get_summary_stats(BB_welch_z, type = 'mean_sd')


# --- offset_z --- #
#MM alone
df_peri_mol %>%
  group_by(MM) %>% 
  get_summary_stats(offset_z, type = "mean_sd")

#MM and progression 
df_peri_prog %>%
  group_by(progression, MM) %>% 
  get_summary_stats(offset_z, type = 'mean_sd')

##MM and epilepsy
df_peri_epi %>%
  group_by(epilepsy_aed, MM) %>% 
  get_summary_stats(offset_z, type = 'mean_sd')


#### ----- Visualize ----- ####

# --- BB_welch_z --- #
#MM alone
bxp_BB <- ggboxplot(df, x = 'MM', y ='BB_welch_z', add = 'point')
bxp_BB

#MM and progression 
bxp_BB_both <- ggboxplot(df, x = 'MM', y = 'BB_welch_z', color = 'progression')
bxp_BB_both

#MM and epilepsy 
bxp_BB_both <- ggboxplot(df, x = 'MM', y = 'BB_welch_z', color = 'epilepsy_aed')
bxp_BB_both

#MM and grade 
bxp_BB_both <- ggboxplot(df_peri_rth_grade, x = 'MM', y = 'BB_welch_z', color = 'graad')
bxp_BB_both

#MM and rth 
bxp_BB_both <- ggboxplot(df_peri_rth_grade, x = 'MM', y = 'BB_welch_z', color = "Syntax_dummy_RTH_all")
bxp_BB_both

# --- offset_z --- #
#MM_alone
bxp_offset <- ggboxplot(df, x = 'MM', y ='offset_z', add = 'point')
bxp_offset

#MM and progression 
bxp_offset_both <- ggboxplot(df, x = 'MM', y = 'offset_z', color = 'progression')
bxp_offset_both

#MM and epilepsy 
bxp_offset_both <- ggboxplot(df, x = 'MM', y = 'offset_z', color = 'epilepsy_aed')
bxp_offset_both

#MM and grade
bxp_offset_both <- ggboxplot(df_peri_rth_grade, x = 'MM', y = 'offset_z', color = 'graad')
bxp_offset_both





#### ----- Two-way repeated measures ANOVA (without taking into consideration transformations) ----- ####
# --- BB_welch_z --- #
#MM
res_aov_BB <- anova_test(data = df_peri_prog, dv = 'BB_welch_z', wid = 'sub', within = 'MM_factor')
get_anova_table(res_aov_BB)

df_res_aov_BB <- data.frame(res_aov_BB)


#Progression
#simple implementation
res_aov_BB_progression <- anova_test(data = df_peri_prog, dv = 'BB_welch_z', wid = 'sub', within = 'MM_factor', between = 'progression_factor')
get_anova_table(res_aov_BB_progression)

df_res_aov_BB_progression <- data.frame(res_aov_BB_progression)


#Epilepsy
#simple implementation
res_aov_BB_epilepsy <- anova_test(data = df_peri_epi, dv = 'BB_welch_z', wid = 'sub', within = 'MM_factor', between = 'epilepsy_aed_factor')
get_anova_table(res_aov_BB_epilepsy)

df_res_aov_BB_epilepsy <- data.frame(res_aov_BB_epilepsy)


#Molecular marker 
#simple implementation
res_aov_BB_mol <- anova_test(data = df_peri_mol, dv = 'BB_welch_z', wid = 'sub', within = 'MM_factor', between = 'IDH_1p19q_factor')
get_anova_table(res_aov_BB_mol)

df_res_aov_BB_mol <- data.frame(res_aov_BB_mol)

### Post hoc test (grade & RTH as covariates) ###
#Tumor grade 
#simple implementation
res_aov_BB_grade <- anova_test(data = df_peri_rth_grade, dv = 'BB_welch_z', wid = 'sub', within = 'MM_factor', between = 'graad_factor')
get_anova_table(res_aov_BB_grade)

df_res_aov_BB_grade <- data.frame(res_aov_BB_grade)

#RTH 
#simple implementation
res_aov_BB_RTH <- anova_test(data = df_peri_rth_grade, dv = 'BB_welch_z', wid = 'sub', within = 'MM_factor', between = 'Syntax_dummy_RTH_all_factor')
get_anova_table(res_aov_BB_RTH)

df_res_aov_BB_RTH <- data.frame(res_aov_BB_RTH)



# --- offset_z --- #
#MM
res_aov_offset <- anova_test(data = df_peri_prog, dv = 'offset_z', wid = 'sub', within = 'MM_factor')
get_anova_table(res_aov_BB)

df_res_aov_offset <- data.frame(res_aov_offset)

#Progression
#simple implementation
res_aov_offset_progression <- anova_test(data = df_peri_prog, dv = 'offset_z', wid = 'sub', within = 'MM_factor', between = 'progression_factor')
get_anova_table(res_aov_offset_progression)

df_res_aov_offset_progression <- data.frame(res_aov_offset_progression)


#Epilepsy
#simple implementation
res_aov_offset_epilepsy <- anova_test(data = df_peri_epi, dv = 'offset_z', wid = 'sub', within = 'MM_factor', between = 'epilepsy_aed_factor')
get_anova_table(res_aov_offset_epilepsy)

df_res_aov_offset_epilepsy <- data.frame(res_aov_offset_epilepsy)


#Molecular marker 
#simple implementation
res_aov_offset_mol <- anova_test(data = df_mol, dv = 'offset_z', wid = 'sub', within = 'MM_factor', between = 'IDH_1p19q_factor')
get_anova_table(res_aov_offset_mol)

df_res_aov_offset_mol <- data.frame(res_aov_offset_mol)


### Post hoc test (grade & RTH as covariates) ###
#Tumor grade 
#simple implementation
res_aov_offset_grade <- anova_test(data = df_peri_rth_grade, dv = 'offset_z', wid = 'sub', within = 'MM_factor', between = 'graad_factor')
get_anova_table(res_aov_offset_grade)

df_res_aov_offset_grade <- data.frame(res_aov_offset_grade)

# --> Grade is a significant main effect. So do pairwise t-tests to understand which grades differ from each other.
pwc_graad <- df_peri_rth_grade %>%
  pairwise_wilcox_test(
    offset_z ~ graad_factor, 
    p.adjust.method = "fdr",
  )
pwc_graad

#Assumptions
#Normality
df_peri_rth_grade %>%  
  group_by(graad_factor) %>%  
  shapiro_test(offset_z)
#--> normality violated for grade 2, use pairwise Mann-Whitney U tests instead

#get the median and IQR of all groups
df_peri_rth_grade %>%  
  group_by(MM_factor, graad_factor) %>%  
  summarise(
    Median = median(offset_z),
    IQR = IQR(offset_z)
  )  

#plot the interaction even though not significant
ggline(df_peri_rth_grade, 
       x = "MM_factor", y = "offset_z", 
       color = "graad_factor",
       add = c("mean_se"),      
       palette = "jco",         
       ylab = "Offset dev", 
       xlab = "Time Point",
       title = "Offset dev by Tumor grade Across Time") +
  theme_minimal() +
  labs(color = "Tumor grade") +
  scale_x_discrete(
    labels = c("T1_Baseline" = "T1", "T2_FU" = "T2")) +
  scale_color_discrete(
    labels = c("II" = "2", "III" = "3", "IV" = "4")) 


#RTH 
#simple implementation
res_aov_offset_RTH <- anova_test(data = df_peri_rth_grade, dv = 'offset_z', wid = 'sub', within = 'MM_factor', between = 'Syntax_dummy_RTH_all_factor')
get_anova_table(res_aov_offset_RTH)

df_res_aov_offset_RTH <- data.frame(res_aov_offset_RTH)
write.csv2(df_res_aov_offset_RTH, 'M:\\MULTINET\\GOALS2\\06_projecten\\2023_activity_EF\\04_analysis\\01_results\\20250701_anova_offset_peritumor_non_transformed_RTH.csv', row.names = TRUE)

df_peri_rth_grade %>%  
  group_by(MM_factor, Syntax_dummy_RTH_all_factor) %>%  
  summarise(
    Mean = mean(offset_z),
    Std = sd(offset_z)
  )

#### ---- Plotting ---- #### 
## Raincloud plots 
#https://github.com/jorvlan/raincloudplots

#install.packages('ggrain')
if (!require(remotes)) {
  install.packages("remotes")
}
remotes::install_github('jorvlan/raincloudplots')

library(raincloudplots)

#Prepare the dataframe and test the scripts 
T1 <- df_peri_prog %>%
  filter(MM == 'T1_Baseline')%>%
  dplyr::select(BB_welch_z) %>% 
  pull() %>%
  as.array()

T2 <- df_peri_prog %>%
  filter(MM == 'T2_FU')%>%
  dplyr::select(BB_welch_z)%>%
  pull() %>%
  as.array()

df_1x1 <- data_1x1(
  array_1 = T1,
  array_2 = T2,
  jit_distance = .09,
  jit_seed = 321)

raincloud_1_v <- raincloud_1x1(
  data = df_1x1, 
  colors = (c('dodgerblue','darkorange')), 
  fills = (c('dodgerblue','darkorange')), 
  size = 1, 
  alpha = .6, 
  ort = 'v') +
  
  scale_x_continuous(breaks=c(1,2), labels=c("Baseline", "FU"), limits=c(0, 3)) +
  xlab("Groups") + 
  ylab("Score") +
  theme_classic()

raincloud_1_v

### Repeated measures ###
raincloud_2 <- raincloud_1x1_repmes(
  data = df_1x1,
  colors = (c('dodgerblue', 'darkorange')),
  fills = (c('dodgerblue', 'darkorange')),
  line_color = 'gray',
  line_alpha = .3,
  size = 1,
  alpha = .6,
  align_clouds = FALSE) +
  
  scale_x_continuous(breaks=c(1,2), labels=c("Pre", "Post"), limits=c(0, 3)) +
  xlab("Time") + 
  ylab("Score") +
  theme_classic()

raincloud_2


### --- With progression groups --- ### 
T1_no_prog <- df_peri_prog %>%
  filter((MM == 'T1_Baseline') & (progression == 'no_progression'))%>%
  dplyr::select(BB_welch_z) %>% 
  pull() %>%
  as.array()

T1_prog <-  df_peri_prog %>%
  filter((MM == 'T1_Baseline') & (progression == 'progression'))%>%
  dplyr::select(BB_welch_z) %>% 
  pull() %>%
  as.array()

T2_no_prog <- df_peri_prog %>%
  filter((MM == 'T2_FU') & (progression == 'no_progression'))%>%
  dplyr::select(BB_welch_z) %>% 
  pull() %>%
  as.array()

T2_prog <- df_peri_prog %>%
  filter((MM == 'T2_FU') & (progression == 'progression'))%>%
  dplyr::select(BB_welch_z) %>% 
  pull() %>%
  as.array()

df_2x2 <- data_2x2(
  array_1 = T1_no_prog,
  array_2 = T1_prog,
  array_3 = T2_no_prog,
  array_4 = T2_prog,
  labels = (c('no_progression','progression')),
  jit_distance = .09,
  jit_seed = 321,
  spread_x_ticks = FALSE)  

#RAINCLOUD
raincloud_2x2 <- raincloud_2x2_repmes(
  data = df_2x2,
  colors = (c('dodgerblue', 'darkorange', 'dodgerblue', 'darkorange')),
  fills = (c('dodgerblue', 'darkorange', 'dodgerblue', 'darkorange')),
  size = 1,
  alpha = .6,
  spread_x_ticks = FALSE) +
  
  scale_x_continuous(breaks=c(1,2), labels=c("Pre", "Post"), limits=c(0, 3)) +
  xlab("Time") + 
  ylab("Score") +
  theme_classic()

raincloud_2x2


### --- With epilepsy groups --- ### 
T1_epilepsy <- df %>%
  filter((MM == 'T1_Baseline') & (epilepsy_aed == 'epilepsy_aed'))%>%
  dplyr::select(BB_welch_z) %>% 
  pull() %>%
  as.array()

T1_no_epilepsy <-  df %>%
  filter((MM == 'T1_Baseline') & (epilepsy_aed == 'no_epilepsy'))%>%
  dplyr::select(BB_welch_z) %>% 
  pull() %>%
  as.array()

T2_epilepsy <- df %>%
  filter((MM == 'T2_FU') & (epilepsy_aed == 'epilepsy_aed'))%>%
  dplyr::select(BB_welch_z) %>% 
  pull() %>%
  as.array()

T2_no_epilepsy <- df %>%
  filter((MM == 'T2_FU') & (epilepsy_aed == 'no_epilepsy'))%>%
  dplyr::select(BB_welch_z) %>% 
  pull() %>%
  as.array()

df_2x2 <- data_2x2(
  array_1 = T1_epilepsy,
  array_2 = T1_no_epilepsy,
  array_3 = T2_epilepsy,
  array_4 = T2_no_epilepsy,
  labels = (c('epilepsy','no epilepsy')),
  jit_distance = .09,
  jit_seed = 321,
  spread_x_ticks = FALSE)  

#RAINCLOUD
raincloud_2x2 <- raincloud_2x2_repmes(
  data = df_2x2,
  colors = (c('dodgerblue', 'darkorange', 'dodgerblue', 'darkorange')),
  fills = (c('dodgerblue', 'darkorange', 'dodgerblue', 'darkorange')),
  size = 1,
  alpha = .6,
  spread_x_ticks = FALSE) +
  
  scale_x_continuous(breaks=c(1,2), labels=c("Pre", "Post"), limits=c(0, 3)) +
  xlab("Time") + 
  ylab("Score") +
  theme_classic()

raincloud_2x2