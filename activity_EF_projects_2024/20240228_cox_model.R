### COX proportional hazards analysis ### 
#following this tutorial: http://www.sthda.com/english/wiki/cox-proportional-hazards-model

install.packages(c("survival", "survminer"))

library("survival")
library("survminer")
library(tidyverse)

### Load Data
path <- 'M:\\MULTINET\\GOALS2\\06_projecten\\2023_activity_EF\\04_analysis\\03_dataframes\\20240228_dataframe_cavity_delta_perc_change_no_prog_covs.csv'

df_cavity <- read_csv(path)

### Factorize categorical and ordinal covariates
df_cavity$sex_factor <- factor(df_cavity$sex)
df_cavity$IDH_1p19q_factor <- factor(df_cavity$IDH_1p19q)
df_cavity$kps_total_factor <- factor(df_cavity$kps_total, ordered = TRUE)

#make event variable binary 
df_cavity$progr_bin <- ifelse(df_cavity$progr == 'yes', 1, 0)

### Cox model ###
#Absolute change
#BB_welch_z
cox_model_BB <- coxph(Surv(PFS_weeks_MEG_FU, progr_bin) ~ BB_welch_z_delta + age_at_diagnosis + sex_factor+ IDH_1p19q_factor, data = df_cavity)
summary(cox_model_BB)

#Offset_z
cox_model_offset <- coxph(Surv(PFS_weeks_MEG_FU, progr_bin) ~ offset_z_delta + age_at_diagnosis + sex_factor + IDH_1p19q_factor, data = df_cavity)
summary(cox_model_offset)

#Percentage change
#BB_welch_z
cox_model_BB_perc_change <- coxph(Surv(PFS_weeks_MEG_FU, progr_bin) ~ BB_welch_z_perc_change + age_at_diagnosis + sex_factor + IDH_1p19q_factor, data = df_cavity)
summary(cox_model_BB_perc_change)

#Offset_z
cox_model_offset_perc_change <- coxph(Surv(PFS_weeks_MEG_FU, progr_bin) ~ offset_z_perc_change + age_at_diagnosis + sex_factor + IDH_1p19q_factor, data = df_cavity)
summary(cox_model_offset_perc_change)
