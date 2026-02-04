# Merge 2025 RC and v3d Data (drop subjects already in the 2023 dataset)
rm(list = setdiff(ls(), c("script_files", "scripts_dir", "log_file", "script_path")))

library(stringr)
library(fst)
library(tidyr)
library(dplyr)
library(here)

# Load data 
load("data/df_v3d_2025_preprocessed.RData") # 2025 v3d data 
load("data/df_RC_2025_preprocessed.RData") # 2025 RC data 
load("data/df_v3d_2023_preprocessed.RData") # 2023 v3d data 

# NOTE: there are some subjects in the 2025 data set that are also in the 2023 data set. We will identify the overlapping subjects and remove them from this dataset so there are no duplicates when we merge the 2025 and 2023 data sets.

# Check for common subject ids in the v3d files of the 2025 and 2023 dfs
common_subjects <- intersect(df_MWF_2025$subject, df_MWF_2023$subject)

# Filter out common ids 
df_MWF_2025_unique <- df_MWF_2025 %>%
  filter(!(subject %in% common_subjects)) %>%
  filter(str_split(subject, "-", simplify = TRUE)[,3] == "1") #selects only ids with "1" in the third section of the id code (i.e. visit number)

# Merge 2025 RedCap and V3d data 
df_2025 <- merge(df_MWF_2025_unique, df_RC_2025, by=c('subject','signal_side'), all.x = TRUE)  

# This data set is already filtered to only include 1 leg (the worse knee) per subject. 
 
# Check for  subjects with missing OA severity measures
missing_severity <- df_2025 %>%
  filter(is.na(severity)) %>%        # Focus on rows where kl_obs2 is NA
  group_by(subject) %>%             # Ensure one row per subject
  summarise(
    severity = first(severity),  
    signal_side = first(signal_side), 
    n = n()                         # Number of rows per subject (for debugging)
  ) %>%
  distinct(subject, .keep_all = TRUE)  # Keep unique subjects

print(missing_severity) # no missing severity measures  

# Check for  subjects with missing kl_grades
missing_kl <- df_2025 %>%
  filter(is.na(kl_score)) %>%        # Focus on rows where kl_obs2 is NA
  group_by(subject) %>%             # Ensure one row per subject
  summarise(
    kl_score = first(kl_score),  
    signal_side = first(signal_side), 
    n = n()                         # Number of rows per subject (for debugging)
  ) %>%
  distinct(subject, .keep_all = TRUE)  # Keep unique subjects

print(missing_kl) # no missing severity measures 

df_2025 <- df_2025 %>%
  filter(!(subject %in% missing_severity$subject)) # exclude subjects with missing data 
  #filter(!(subject %in% missing_kl$subject))

print(length(unique(df_2025$subject)))

# Save output 
save(df_2025, file="data/df_merged_RC_v3d_2025.RData")

