#title: "Post-hoc ANOVA test Peri/resection cavity area"
#author: "Mona Zimmermann"
#date: ""
#output: html_document

#### Email: m.l.m.zimmermann@amsterdamumc.nl
#### Purpose of script: Statistical analysis for the post hoc test of clinical 
#### relevance in activity changes for the peri-resection cavity

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

# https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/
# mixed anova: https://www.datanovia.com/en/lessons/mixed-anova-in-r/


##### ----- Data preparation ----- #####
file_path_cavity_progression <- 'M:\\path\\to\\20250404_dataframe_cavity_baseline_FU_long_format_cov_prog.csv'
file_path_cavity_epilepsy <- 'M:\\path\\to\20250404_dataframe_cavity_baseline_FU_long_format_cov_epi.csv'
file_path_cavity_mol <- 'M:\\path\\to\20250404_dataframe_cavity_baseline_FU_long_format_cov_mol.csv'
file_path_cavity_rth_grade <- 'M:\\path\\to\\20250701_dataframe_cavity_baseline_FU_long_format_cov_rth_grade.csv'

df_cavity_prog <- read_csv(file_path_cavity_progression)
df_cavity_epi <- read_csv(file_path_cavity_epilepsy)
df_cavity_mol <- read_csv(file_path_cavity_mol)
df_cavity_rth_grade <- read_csv(file_path_cavity_rth_grade)

#turn necessary categorical columns into factors to use them as covariates in the analysis (returns same results as anova already interpreted the categorical vars as factors)

factorize <-function(df, var_name){
  df$MM_factor <- factor(df$MM)
  df[[paste0(var_name, "_factor")]] <- factor(df[[var_name]])
  return(df)
}

df_cavity_prog <- factorize(df_cavity_prog, "progression")
df_cavity_epi <- factorize(df_cavity_epi, "epilepsy_aed")
df_cavity_mol <- factorize(df_cavity_mol, "IDH_1p19q")
df_cavity_rth_grade <- factorize(df_cavity_rth_grade, "Syntax_dummy_RTH_all")
df_cavity_rth_grade <- factorize(df_cavity_rth_grade, "graad")
  



#### ----- Summary statistics ----- ####

# --- BB_welch_z --- #
#MM alone 
df_cavity_prog %>%
  group_by(MM) %>% 
  get_summary_stats(BB_welch_z, type = "mean_sd")


#MM and progression 
df_cavity_prog %>%
  group_by(progression, MM) %>% 
  get_summary_stats(BB_welch_z, type = 'mean_sd')

##MM and epilepsy
df_cavity_epi %>%
  group_by(epilepsy_aed, MM) %>% 
  get_summary_stats(BB_welch_z, type = 'mean_sd')

##MM and molecular status
df_cavity_mol %>%
  group_by(IDH_1p19q, MM) %>%
  get_summary_stats(BB_welch_z, type = 'mean_sd')


# --- offset_z --- #
#MM alone
df_cavity_prog %>%
  group_by(MM) %>% 
  get_summary_stats(offset_z, type = "mean_sd")

#MM and progression 
df_cavity_prog %>%
  group_by(progression, MM) %>% 
  get_summary_stats(offset_z, type = 'mean_sd')

##MM and epilepsy
df_cavity_epi %>%
  group_by(epilepsy_aed, MM) %>% 
  get_summary_stats(offset_z, type = 'mean_sd')

##MM and molecular status
df_cavity_mol %>%
  group_by(IDH_1p19q, MM) %>%
  get_summary_stats(offset_z, type = 'mean_sd')



#### ----- Two-way repeated measures ANOVA (without taking into considertation transformations) ----- ####
# --- BB_welch_z --- #
#MM
res_aov_BB <- anova_test(data = df_cavity_prog, dv = 'BB_welch_z', wid = 'sub', within = 'MM_factor')
get_anova_table(res_aov_BB)

df_res_aov_BB <- data.frame(res_aov_BB)


#Progression
#simple implementation
res_aov_BB_progression <- anova_test(data = df_cavity_prog, dv = 'BB_welch_z', wid = 'sub', within = 'MM', between = 'progression_factor')
get_anova_table(res_aov_BB_progression)

df_res_aov_BB_progression <- data.frame(res_aov_BB_progression)

#Epilepsy
#simple implementation
res_aov_BB_epilepsy <- anova_test(data = df_cavity_epi, dv = 'BB_welch_z', wid = 'sub', within = 'MM', between = 'epilepsy_aed_factor')
get_anova_table(res_aov_BB_epilepsy)

df_res_aov_BB_epilepsy <- data.frame(res_aov_BB_epilepsy)

#Molecular status
#simple implementation

res_aov_BB_mol <- anova_test(data = df_cavity_mol, dv = 'BB_welch_z', wid = 'sub', within = 'MM', between = 'IDH_1p19q_factor')
get_anova_table(res_aov_BB_mol)

df_res_aov_BB_mol <- data.frame(res_aov_BB_mol)

### Post-hoc test: RTH and grade ###
#RTH
#simple implementation
res_aov_BB_rth <- anova_test(data = df_cavity_rth_grade, dv = 'BB_welch_z', wid = 'sub', within = 'MM', between = 'Syntax_dummy_RTH_all_factor')
get_anova_table(res_aov_BB_rth)

df_res_aov_BB_rth <- data.frame(res_aov_BB_rth)

df_cavity_rth_grade %>%  
  group_by(MM_factor, Syntax_dummy_RTH_all_factor) %>%  
  summarise(
    Mean = mean(BB_welch_z),
    Std = sd(BB_welch_z)
  )

#grade
#simple implementation
res_aov_BB_grade <- anova_test(data = df_cavity_rth_grade, dv = 'BB_welch_z', wid = 'sub', within = 'MM', between = 'graad_factor')
get_anova_table(res_aov_BB_grade)

df_res_aov_BB_grade <- data.frame(res_aov_BB_grade)



# --- Offset_z --- #
#MM
res_aov_offset <- anova_test(data = df_cavity_prog, dv = 'offset_z', wid = 'sub', within = 'MM_factor')
get_anova_table(res_aov_offset)

df_res_aov_offset <- data.frame(res_aov_offset)

#Progression
#simple implementation
res_aov_offset_progression <- anova_test(data = df_cavity_prog, dv = 'offset_z', wid = 'sub', within = 'MM', between = 'progression_factor')
get_anova_table(res_aov_offset_progression)

df_res_aov_offset_progression <- data.frame(res_aov_offset_progression)

#Epilepsy
#simple implementation
res_aov_offset_epilepsy <- anova_test(data = df_cavity_epi, dv = 'offset_z', wid = 'sub', within = 'MM', between = 'epilepsy_aed_factor')
get_anova_table(res_aov_offset_epilepsy)

df_res_aov_offset_epilepsy <- data.frame(res_aov_offset_epilepsy)

#Molecular status
#simple implementation
res_aov_offset_mol <- anova_test(data = df_cavity_mol, dv = 'offset_z', wid = 'sub', within = 'MM', between = 'IDH_1p19q_factor')
get_anova_table(res_aov_offset_mol)

df_res_aov_offset_mol <- data.frame(res_aov_offset_mol)

### Post-hoc test: RTH and grade ###
#RTH
#simple implementation
res_aov_offset_rth <- anova_test(data = df_cavity_rth_grade, dv = 'offset_z', wid = 'sub', within = 'MM', between = 'Syntax_dummy_RTH_all_factor')
get_anova_table(res_aov_offset_rth)

df_res_aov_offset_rth <- data.frame(res_aov_offset_rth)

df_cavity_rth_grade %>%  
  group_by(MM_factor, Syntax_dummy_RTH_all_factor) %>%  
  summarise(
    Mean = mean(offset_z),
    Std = sd(offset_z)
  )

#grade
#simple implementation
res_aov_offset_grade <- anova_test(data = df_cavity_rth_grade, dv = 'offset_z', wid = 'sub', within = 'MM', between = 'graad_factor')
get_anova_table(res_aov_offset_grade)

df_res_aov_offset_grade <- data.frame(res_aov_offset_grade)



