library(tidyverse) 
library(ggpubr)
library(rstatix)
library(MASS)

# https://www.datanovia.com/en/lessons/repeated-measures-anova-in-r/
# mixed anova: https://www.datanovia.com/en/lessons/mixed-anova-in-r/


##### ----- Data preparation ----- #####
file_path_progression <- 'M:\\MULTINET\\GOALS2\\06_projecten\\2023_activity_EF\\04_analysis\\03_dataframes\\20240221_dataframe_peritumoral_baseline_FU_long_format_progression.csv'
file_path_epilepsy <- 'M:\\MULTINET\\GOALS2\\06_projecten\\2023_activity_EF\\04_analysis\\03_dataframes\\20240221_dataframe_peritumoral_baseline_FU_long_format_epilepsy_filtered.csv'
file_path_mol <- 'M:\\MULTINET\\GOALS2\\06_projecten\\2023_activity_EF\\04_analysis\\03_dataframes\\20240221_dataframe_peritumoral_baseline_FU_long_format_mol.csv'


df_progression <- read_csv(file_path_progression)
df_epilepsy <- read_csv(file_path_epilepsy)
df_mol <- read_csv(file_path_mol)

#### ----- Summary statistics ----- ####

# --- BB_welch_z --- #
#MM alone 
df_mol %>%
  group_by(MM) %>% 
  get_summary_stats(BB_welch_z, type = "mean_sd")

#MM and progression 
df_progression %>%
  group_by(progression, MM) %>% 
  get_summary_stats(BB_welch_z, type = 'mean_sd')

##MM and epilepsy
df_epilepsy %>%
  group_by(epilepsy_aed, MM) %>% 
  get_summary_stats(BB_welch_z, type = 'mean_sd')


# --- offset_z --- #
#MM alone
df_mol %>%
  group_by(MM) %>% 
  get_summary_stats(offset_z, type = "mean_sd")

#MM and progression 
df_progression %>%
  group_by(progression, MM) %>% 
  get_summary_stats(offset_z, type = 'mean_sd')

##MM and epilepsy
df_epilepsy %>%
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

##### Assumption checks #####
### Outliers (NA here)
### Normality assumption 

# --- BB_welch_z --- #
#Progression
df_progression %>%
  group_by(progression, MM) %>% 
  shapiro_test(BB_welch_z)
#Not normally distributed: Progression, T2_FU p = 0.0004


#Epilepsy
df_epilepsy %>%
  group_by(epilepsy_aed, MM) %>% 
  shapiro_test(BB_welch_z)
#Not normally distributed: Epilepsy_aed, T1_Baseline (0.0102) and T2_FU (0.007)




# --- offset_z --- #
#Progression
df_progression %>%
  group_by(progression, MM) %>% 
  shapiro_test(offset_z)
#Not normally distributed: 
#Non-progression, T2_FU p = 0.0009
#Progression, T2_FU p = 0.0160


#Epilepsy
df_epilepsy %>%
  group_by(epilepsy_aed, MM) %>% 
  shapiro_test(offset_z)
#Not normally distributed: 
#epilepsy_aed, T2_FU p = 0.0001




### Assumption of sphericity 
#(Greenhouse-Geisser correction automatically applied to factors violating 
#the assumption when using get_anova_table)

### Assumption of homogeneity of variances 
# --- BB_welch_z --- #
#Progression
df_progression %>% 
  group_by(MM) %>%
  levene_test(BB_welch_z ~ progression)

#Epilepsy
df_epilepsy %>% 
  group_by(MM) %>%
  levene_test(BB_welch_z ~ epilepsy_aed)


# --- offset_z --- # 
#Progression
df_progression %>% 
  group_by(MM) %>% 
  levene_test(offset_z ~ progression)

#Epilepsy             
df_epilepsy %>% 
  group_by(MM) %>% 
  levene_test(offset_z ~ epilepsy_aed)

### Assumption of homogeneity of covariances (covariances should be equal
### across the cells formed by the between-subject factors)
#Progression

# --- BB_welch_z --- #
box_m(df_progression[, "BB_welch_z", drop = FALSE], df_progression$progression)

# --- offset_z --- #
box_m(df_progression[, "offset_z", drop = FALSE], df_progression$progression)

#Epilepsy
# --- BB_welch_z --- #
box_m(df_epilepsy[, "BB_welch_z", drop = FALSE], df_epilepsy$epilepsy_aed)

# --- offset_z --- #
box_m(df_epilepsy[, "offset_z", drop = FALSE], df_epilepsy$epilepsy_aed)


#### ----- Two-way repeated measures ANOVA (without taking into considertation transformations) ----- ####
# --- BB_welch_z --- #
#MM
res_aov_BB <- anova_test(data = df_progression, dv = 'BB_welch_z', wid = 'sub', within = 'MM')
get_anova_table(res_aov_BB)

df_res_aov_BB <- data.frame(res_aov_BB)
write.csv2(df_res_aov_BB, 'M:\\MULTINET\\GOALS2\\06_projecten\\2023_activity_EF\\04_analysis\\01_results\\20240327_anova_BB_peritumor_non_transformed_simple_implementation.csv', row.names = TRUE)



#Progression
#simple implementation
res_aov_BB_progression <- anova_test(data = df_progression, dv = 'BB_welch_z', wid = 'sub', within = 'MM', between = 'progression')
get_anova_table(res_aov_BB_progression)

df_res_aov_BB_progression <- data.frame(res_aov_BB_progression)
write.csv2(df_res_aov_BB_progression, 'M:\\MULTINET\\GOALS2\\06_projecten\\2023_activity_EF\\03_analysis\\01_results\\20240221_anova_BB_peritumor_non_transformed_progression.csv', row.names = TRUE)


#Epilepsy
#simple implementation
res_aov_BB_epilepsy <- anova_test(data = df_epilepsy, dv = 'BB_welch_z', wid = 'sub', within = 'MM', between = 'epilepsy_aed')
get_anova_table(res_aov_BB_epilepsy)

df_res_aov_BB_epilepsy <- data.frame(res_aov_BB_epilepsy)
write.csv2(df_res_aov_BB_epilepsy, 'M:\\MULTINET\\GOALS2\\06_projecten\\2023_activity_EF\\03_analysis\\01_results\\20240221_anova_BB_peritumor_non_transformed_epilepsy.csv', row.names = TRUE)



#Molecular marker 
#simple implementation
res_aov_BB_mol <- anova_test(data = df_mol, dv = 'BB_welch_z', wid = 'sub', within = 'MM', between = 'IDH_1p19q')
get_anova_table(res_aov_BB_mol)

df_res_aov_BB_mol <- data.frame(res_aov_BB_mol)
write.csv2(df_res_aov_BB_mol, 'M:\\MULTINET\\GOALS2\\06_projecten\\2023_activity_EF\\03_analysis\\01_results\\20240221_anova_BB_peritumor_non_transformed_mol.csv', row.names = TRUE)



# --- offset_z --- #
#MM
res_aov_BB <- anova_test(data = df_progression, dv = 'offset_z', wid = 'sub', within = 'MM')
get_anova_table(res_aov_BB)

df_res_aov_BB <- data.frame(res_aov_BB)
write.csv2(df_res_aov_BB, 'M:\\MULTINET\\GOALS2\\06_projecten\\2023_activity_EF\\03_analysis\\01_results\\20240221_anova_offset_peritumor_non_transformed.csv', row.names = TRUE)

#Progression
#simple implementation
res_aov_offset_progression <- anova_test(data = df_progression, dv = 'offset_z', wid = 'sub', within = 'MM', between = 'progression')
get_anova_table(res_aov_offset_progression)

df_res_aov_offset_progression <- data.frame(res_aov_offset_progression)
write.csv2(df_res_aov_offset_progression, 'M:\\MULTINET\\GOALS2\\06_projecten\\2023_activity_EF\\03_analysis\\01_results\\20240221_anova_offset_peritumor_non_transformed_progression.csv', row.names = FALSE)


#Epilepsy
#simple implementation
res_aov_offset_epilepsy <- anova_test(data = df_epilepsy, dv = 'offset_z', wid = 'sub', within = 'MM', between = 'epilepsy_aed')
get_anova_table(res_aov_offset_epilepsy)

df_res_aov_offset_epilepsy <- data.frame(res_aov_offset_epilepsy)
write.csv2(df_res_aov_offset_epilepsy, 'M:\\MULTINET\\GOALS2\\06_projecten\\2023_activity_EF\\03_analysis\\01_results\\20240221_anova_offset_peritumor_non_transformed_epilepsy.csv', row.names = FALSE)


#Molecular marker 
#simple implementation
res_aov_offset_mol <- anova_test(data = df_mol, dv = 'offset_z', wid = 'sub', within = 'MM', between = 'IDH_1p19q')
get_anova_table(res_aov_offset_mol)

df_res_aov_offset_mol <- data.frame(res_aov_offset_mol)
write.csv2(df_res_aov_offset_mol, 'M:\\MULTINET\\GOALS2\\06_projecten\\2023_activity_EF\\03_analysis\\01_results\\20240221_anova_offset_peritumor_non_transformed_mol.csv', row.names = TRUE)



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
T1 <- df %>%
  filter(MM == 'T1_Baseline')%>%
  dplyr::select(BB_welch_z) %>% 
  pull() %>%
  as.array()

T2 <- df %>%
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
T1_no_prog <- df %>%
  filter((MM == 'T1_Baseline') & (progression == 'no_progression'))%>%
  dplyr::select(BB_welch_z) %>% 
  pull() %>%
  as.array()

T1_prog <-  df %>%
  filter((MM == 'T1_Baseline') & (progression == 'progression'))%>%
  dplyr::select(BB_welch_z) %>% 
  pull() %>%
  as.array()

T2_no_prog <- df %>%
  filter((MM == 'T2_FU') & (progression == 'no_progression'))%>%
  dplyr::select(BB_welch_z) %>% 
  pull() %>%
  as.array()

T2_prog <- df %>%
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