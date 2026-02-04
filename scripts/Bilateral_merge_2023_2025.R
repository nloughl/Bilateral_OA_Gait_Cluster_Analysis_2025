# Combine the 2023 and 2025 pre-processed data sets 
rm(list = setdiff(ls(), c("script_files", "scripts_dir", "log_file", "script_path")))

library(tidyr)
library(dplyr)
library(here)
library(ggplot2)
library(stringr)

load("data/df_one_leg_2025.RData") #2025 data 
load("data/df_one_leg_2023.RData") #2023 data 
#load("data/df_one_leg_2023_w_controls.RData") #2023 data with healthy/asym controls
#df_one_leg_2023 <- df_one_leg_2023_with_controls # rename dataset with controls to go through script 

# Check for common subject ids in the 2025 and 2023 dfs (should be none)
common_subjects <- intersect(df_one_leg_2023$subject, df_one_leg_2025$subject)

# Ensure same column names 
colnames <- c("subject", "action", "signal_names", "signal_side", "signal_components", "item", "value", "speed", "age", "sex", "oa_location", "severity", "severity_contra", "kl_score", "kl_contra")

df_one_leg_2023 <- df_one_leg_2023 %>%
  dplyr::select(all_of(colnames))  # Keep only specified columns

df_one_leg_2025 <- df_one_leg_2025 %>%
  dplyr::select(all_of(colnames))  # Keep only specified columns

df_one_leg_2025 <- df_one_leg_2025 %>%
  mutate(age = as.numeric(age)) %>% # Convert to numeric
  mutate(kl_score = as.numeric(kl_score)) %>%
  mutate(kl_contra = as.numeric(kl_contra))

# Merge 2025 RedCap and V3d data 
df_combined <- bind_rows(df_one_leg_2023, df_one_leg_2025)  

print(length(unique(df_combined$subject)))

# Save data frame 
save(df_combined, file="data/df_2023_2025_combined.RData")
