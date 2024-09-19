#title: "Post-hoc ANOVA test Enhancing tumor area"
#author: "Mona Zimmermann"
#date: ""
#output: html_document

#### Email: m.l.m.zimmermann@amsterdamumc.nl
#### Purpose of script: Statistical analysis for the post hoc test of clinical 
#### relevance in activity changes for the enhancing tumor area

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
file_path_enhancing_tumor <- 'M:\\MULTINET\\GOALS2\\06_projecten\\2023_activity_EF\\04_analysis\\03_dataframes\\20240221_dataframe_enhancing_tumor_baseline_FU_long_format_progression_epilepsy_mol.csv'
file_path_enhancing_tumor_epilepsy <- 'M:\\MULTINET\\GOALS2\\06_projecten\\2023_activity_EF\\04_analysis\\03_dataframes\\20240221_dataframe_enhancing_tumor_baseline_FU_long_format_epilepsy_filtered.csv'

df_enhancing_tumor <- read_csv(file_path_enhancing_tumor)
df_epilepsy <- read_csv(file_path_enhancing_tumor_epilepsy)

#turn necessary categorical columns into factors to use them as covariates in the analysis (returns same results as anova already interpreted the categorical vars as factors)

#df_enhancing_tumor$MM_factor <- factor(df_enhancing_tumor$MM)
#df_enhancing_tumor$progression_factor <- factor(df_enhancing_tumor$progression)
#df_enhancing_tumor$IDH_1p19q_factor <- factor(df_enhancing_tumor$IDH_1p19q)
#df_enhancing_tumor$epilepsy_aed_factor <- factor(df_enhancing_tumor$epilepsy_aed)


#### ----- Summary statistics ----- ####

# --- BB_welch_z --- #
#MM alone 
df_enhancing_tumor %>%
  group_by(MM) %>% 
  get_summary_stats(BB_welch_z, type = "mean_sd")

#MM and progression 
df_enhancing_tumor %>%
  group_by(progression, MM) %>% 
  get_summary_stats(BB_welch_z, type = 'mean_sd')

##MM and epilepsy
df_enhancing_tumor %>%
  group_by(epilepsy_aed, MM) %>% 
  get_summary_stats(BB_welch_z, type = 'mean_sd')

##MM and molecular status
df_enhancing_tumor %>%
  group_by(IDH_1p19q, MM) %>%
  get_summary_stats(BB_welch_z, type = 'mean_sd')


# --- offset_z --- #
#MM alone
df_enhancing_tumor %>%
  group_by(MM) %>% 
  get_summary_stats(offset_z, type = "mean_sd")

#MM and progression 
df_enhancing_tumor %>%
  group_by(progression, MM) %>% 
  get_summary_stats(offset_z, type = 'mean_sd')

##MM and epilepsy
df_enhancing_tumor %>%
  group_by(epilepsy_aed, MM) %>% 
  get_summary_stats(offset_z, type = 'mean_sd')

##MM and molecular status
df_enhancing_tumor %>%
  group_by(IDH_1p19q, MM) %>%
  get_summary_stats(offset_z, type = 'mean_sd')

# --> don't compute differences with molecular status and progression --> for both variables, only one subject has no_progression or an IDH-mutant glioma

#### ----- Two-way repeated measures ANOVA (without taking into considertation transformations) ----- ####
#BB_welch_z
#MM
res_aov_BB <- anova_test(data = df_enhancing_tumor, dv = 'BB_welch_z', wid = 'sub', within = 'MM')
get_anova_table(res_aov_BB)

df_res_aov_BB <- data.frame(res_aov_BB)
write.csv2(df_res_aov_BB, 'M:\\MULTINET\\GOALS2\\06_projecten\\2023_activity_EF\\03_analysis\\01_results\\20240221_anova_BB_enhancing_tumor_non_transformed.csv', row.names = TRUE)

#Epilepsy
#simple implementation
res_aov_BB_epilepsy <- anova_test(data = df_epilepsy, dv = 'BB_welch_z', wid = 'sub', within = 'MM', between = 'epilepsy_aed')
get_anova_table(res_aov_BB_epilepsy)

df_res_aov_BB_epilepsy <- data.frame(res_aov_BB_epilepsy)
write.csv2(df_res_aov_BB_epilepsy, 'M:\\MULTINET\\GOALS2\\06_projecten\\2023_activity_EF\\03_analysis\\01_results\\20240327_anova_BB_enhancing_tumor_non_transformed_epilepsy.csv', row.names = TRUE)


#Offset_z
#MM
res_aov_offset <- anova_test(data = df_enhancing_tumor, dv = 'offset_z', wid = 'sub', within = 'MM')
get_anova_table(res_aov_offset)

df_res_aov_offset <- data.frame(res_aov_offset)
write.csv2(df_res_aov_offset, 'M:\\MULTINET\\GOALS2\\06_projecten\\2023_activity_EF\\03_analysis\\01_results\\20240221_anova_offset_enhancing_tumor_non_transformed.csv', row.names = TRUE)

#Epilepsy
#simple implementation
res_aov_offset_epilepsy <- anova_test(data = df_epilepsy, dv = 'offset_z', wid = 'sub', within = 'MM', between = 'epilepsy_aed')
get_anova_table(res_aov_offset_epilepsy)

df_res_aov_offset_epilepsy <- data.frame(res_aov_offset_epilepsy)
write.csv2(df_res_aov_offset_epilepsy, 'M:\\MULTINET\\GOALS2\\06_projecten\\2023_activity_EF\\03_analysis\\01_results\\20240221_anova_offset_enhancing_tumor_non_transformed_epilepsy.csv', row.names = TRUE)


