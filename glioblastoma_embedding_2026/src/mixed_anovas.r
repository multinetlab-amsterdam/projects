#### Load libraries
library(tidyverse) 
library(ggpubr)
library(rstatix)
library(MASS)
library(dplyr)



file_path <- 'M:\\path\\to\\data\\tumor_conn_analysis.csv'

df <- read_csv(file_path)

df$split_tumor_activity_factor <- factor(df$split_tumor_activity)
df$tumor_conn_binary_factor <- factor(df$tumor_conn_binary)
df$median_split_PATNET_factor <- factor(df$median_split_PATNET)

### TUMOR ACTIVITY as between-subject factor
res_aov_type_3<- anova_test(data = df, dv = 'mean_BBP_z', wid = 'sub', within = 'tumor_conn_binary_factor', between = 'split_tumor_activity_factor', type=3)
get_anova_table(res_aov_type_3)


# --- Assumption tests
# -- Normality ---
# Fit the equivalent lm model to extract residuals
library(nlme)
model_lme <- lme(mean_BBP_z ~ tumor_conn_binary_factor * split_tumor_activity_factor,
                 random = ~1 | sub, data = df)

ggplot(data.frame(residuals = residuals(model_lme)), aes(sample = residuals)) +
  stat_qq() + stat_qq_line(color = "red") +
  labs(title = "Q-Q Plot of Residuals (mixed model)")

library(car)
# --- homogeneity of variance (between subjects) ---
# Levene's test — per level of the between-subjects factor
df %>%
  group_by(tumor_conn_binary_factor) %>%
  levene_test(mean_BBP_z ~ split_tumor_activity_factor)

# Bartlett's test (assumes normality)
bartlett.test(mean_BBP_z ~ split_tumor_activity_factor, data = df)





### PATNET as between subject factor
res_aov_type_3<- anova_test(data = df, dv = 'mean_BBP_z', wid = 'sub', within = 'tumor_conn_binary_factor', between = 'median_split_PATNET_factor', type=3)
get_anova_table(res_aov_type_3)

# --- Assumption tests
# -- Normality ---
# Fit the equivalent lm model to extract residuals
library(nlme)
model_lme <- lme(mean_BBP_z ~ tumor_conn_binary_factor * median_split_PATNET_factor,
                 random = ~1 | sub, data = df)

ggplot(data.frame(residuals = residuals(model_lme)), aes(sample = residuals)) +
  stat_qq() + stat_qq_line(color = "red") +
  labs(title = "Q-Q Plot of Residuals (mixed model)")


# --- homogeneity of variance (between subjects) ---
# Levene's test — per level of the within-subjects factor
df %>%
  group_by(tumor_conn_binary_factor) %>%
  levene_test(mean_BBP_z ~ median_split_PATNET_factor)

# Bartlett's test (assumes normality)
bartlett.test(mean_BBP_z ~ median_split_PATNET_factor, data = df)

