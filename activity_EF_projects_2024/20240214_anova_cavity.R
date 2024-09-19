library(tidyverse) 
library(ggpubr)
library(rstatix)
library(MASS)

# https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/
# mixed anova: https://www.datanovia.com/en/lessons/mixed-anova-in-r/


##### ----- Data preparation ----- #####
file_path_cavity_progression_mol <- 'M:\\MULTINET\\GOALS2\\06_projecten\\2023_activity_EF\\04_analysis\\03_dataframes\\20240221_dataframe_cavity_baseline_FU_long_format_progression_epilepsy_mol.csv'
file_path_cavity_epilepsy <- 'M:\\MULTINET\\GOALS2\\06_projecten\\2023_activity_EF\\04_analysis\\03_dataframes\\20240221_dataframe_cavity_baseline_FU_long_format_epilepsy_filtered.csv'

df_cavity <- read_csv(file_path_cavity_progression_mol)
df_cavity_epilepsy <- read_csv(file_path_cavity_epilepsy)

#turn necessary categorical columns into factors to use them as covariates in the analysis (returns same results as anova already interpreted the categorical vars as factors)

#df_cavity$MM_factor <- factor(df_cavity$MM)
#df_cavity$progression_factor <- factor(df_cavity$progression)
#df_cavity$IDH_1p19q_factor <- factor(df_cavity$IDH_1p19q)
#df_cavity$epilepsy_aed_factor <- factor(df_cavity$epilepsy_aed)



#### ----- Summary statistics ----- ####

# --- BB_welch_z --- #
#MM alone 
df_cavity %>%
  group_by(MM) %>% 
  get_summary_stats(BB_welch_z, type = "mean_sd")


#MM and progression 
df_cavity %>%
  group_by(progression, MM) %>% 
  get_summary_stats(BB_welch_z, type = 'mean_sd')

##MM and epilepsy
df_cavity %>%
  group_by(epilepsy_aed, MM) %>% 
  get_summary_stats(BB_welch_z, type = 'mean_sd')

##MM and molecular status
df_cavity %>%
  group_by(IDH_1p19q, MM) %>%
  get_summary_stats(BB_welch_z, type = 'mean_sd')


# --- offset_z --- #
#MM alone
df_cavity %>%
  group_by(MM) %>% 
  get_summary_stats(offset_z, type = "mean_sd")

#MM and progression 
df_cavity %>%
  group_by(progression, MM) %>% 
  get_summary_stats(offset_z, type = 'mean_sd')

##MM and epilepsy
df_cavity %>%
  group_by(epilepsy_aed, MM) %>% 
  get_summary_stats(offset_z, type = 'mean_sd')

##MM and molecular status
df_cavity %>%
  group_by(IDH_1p19q, MM) %>%
  get_summary_stats(offset_z, type = 'mean_sd')



#### ----- Two-way repeated measures ANOVA (without taking into considertation transformations) ----- ####
# --- BB_welch_z --- #
#MM
res_aov_BB <- anova_test(data = df_cavity, dv = 'BB_welch_z', wid = 'sub', within = 'MM')
get_anova_table(res_aov_BB)

df_res_aov_BB <- data.frame(res_aov_BB)
write.csv2(df_res_aov_BB, 'M:\\MULTINET\\GOALS2\\06_projecten\\2023_activity_EF\\03_analysis\\01_results\\20240221_anova_BB_cavity_non_transformed.csv', row.names = TRUE)


#Progression
#simple implementation
res_aov_BB_progression <- anova_test(data = df_cavity, dv = 'BB_welch_z', wid = 'sub', within = 'MM', between = 'progression')
get_anova_table(res_aov_BB_progression)

df_res_aov_BB_progression <- data.frame(res_aov_BB_progression)
write.csv2(df_res_aov_BB_progression, 'M:\\MULTINET\\GOALS2\\06_projecten\\2023_activity_EF\\03_analysis\\01_results\\20240221_anova_BB_cavity_non_transformed_progression.csv', row.names = TRUE)

#Epilepsy
#simple implementation
res_aov_BB_epilepsy <- anova_test(data = df_cavity_epilepsy, dv = 'BB_welch_z', wid = 'sub', within = 'MM', between = 'epilepsy_aed')
get_anova_table(res_aov_BB_epilepsy)

df_res_aov_BB_epilepsy <- data.frame(res_aov_BB_epilepsy)
write.csv2(df_res_aov_BB_epilepsy, 'M:\\MULTINET\\GOALS2\\06_projecten\\2023_activity_EF\\03_analysis\\01_results\\20240221_anova_BB_cavity_non_transformed_epilepsy.csv', row.names = TRUE)

#Molecular status
#simple implementation

res_aov_BB_mol <- anova_test(data = df_cavity, dv = 'BB_welch_z', wid = 'sub', within = 'MM', between = 'IDH_1p19q')
get_anova_table(res_aov_BB_mol)

df_res_aov_BB_mol <- data.frame(res_aov_BB_mol)
write.csv2(df_res_aov_BB_mol, 'M:\\MULTINET\\GOALS2\\06_projecten\\2023_activity_EF\\03_analysis\\01_results\\20240221_anova_BB_cavity_non_transformed_mol.csv', row.names = TRUE)


# --- Offset_z --- #
#MM
res_aov_offset <- anova_test(data = df_cavity, dv = 'offset_z', wid = 'sub', within = 'MM')
get_anova_table(res_aov_offset)

df_res_aov_offset <- data.frame(res_aov_offset)
write.csv2(df_res_aov_offset, 'M:\\MULTINET\\GOALS2\\06_projecten\\2023_activity_EF\\03_analysis\\01_results\\20240221_anova_offset_cavity_non_transformed.csv', row.names = TRUE)

#Progression
#simple implementation
res_aov_offset_progression <- anova_test(data = df_cavity, dv = 'offset_z', wid = 'sub', within = 'MM', between = 'progression')
get_anova_table(res_aov_offset_progression)

df_res_aov_offset_progression <- data.frame(res_aov_offset_progression)
write.csv2(df_res_aov_offset_progression, 'M:\\MULTINET\\GOALS2\\06_projecten\\2023_activity_EF\\03_analysis\\01_results\\20240221_anova_offset_cavity_non_transformed_progression.csv', row.names = TRUE)

#Epilepsy
#simple implementation
res_aov_offset_epilepsy <- anova_test(data = df_cavity_epilepsy, dv = 'offset_z', wid = 'sub', within = 'MM', between = 'epilepsy_aed')
get_anova_table(res_aov_offset_epilepsy)

df_res_aov_offset_epilepsy <- data.frame(res_aov_offset_epilepsy)
write.csv2(df_res_aov_offset_epilepsy, 'M:\\MULTINET\\GOALS2\\06_projecten\\2023_activity_EF\\03_analysis\\01_results\\20240221_anova_offset_cavity_non_transformed_epilepsy.csv', row.names = TRUE)

#Molecular status
#simple implementation

res_aov_offset_mol <- anova_test(data = df_cavity, dv = 'offset_z', wid = 'sub', within = 'MM', between = 'IDH_1p19q')
get_anova_table(res_aov_offset_mol)

df_res_aov_offset_mol <- data.frame(res_aov_offset_mol)
write.csv2(df_res_aov_offset_mol, 'M:\\MULTINET\\GOALS2\\06_projecten\\2023_activity_EF\\03_analysis\\01_results\\20240221_anova_offset_cavity_non_transformed_mol.csv', row.names = TRUE)
