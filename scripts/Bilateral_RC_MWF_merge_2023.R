# 2023 data set RedCap preprocessing and merging with df_MWF_2023

rm(list = setdiff(ls(), c("script_files", "scripts_dir", "log_file", "script_path")))

library(fst)
library(tidyr)
library(dplyr)
library(here)
library(ggplot2)
library(stringr)

load("data/df_RC_offline_2025-03-13.RData")
load("data/df_v3d_2023_preprocessed.RData")

#Filter RedCap data for columns we want 
df_RC_2023 <- df_RC %>%
  select(signal_side, subject, clothing_bottom, affected_knee, 
         age, sex, kl_obs2, oa_location_obs2, severity) %>% #Using 'severity' as severity measure 
  filter(clothing_bottom != "Skirt/dress below knees") %>% #Filter out entries with long skirts
  rename(oa_location = oa_location_obs2) %>% 
  rename(kl_score = kl_obs2)

# Merge 2023 RedCap and Mean wave form data 
df_2023 <- df_MWF_2023 %>%
  left_join(df_RC_2023, by=c('subject','signal_side'))

# Filter out any patients with only 1 leg entry (need left and right)
table(df_2023$signal_side) #Uneven left and right entries 

df_2023 <- df_2023 %>%
  group_by(subject) %>% # Group by subject ID
  filter(all(c("LEFT", "RIGHT") %in% signal_side)) %>%  # Keep only subjects with both sides
  ungroup()

table(df_2023$signal_side) # check they are now equal!! 

# Check for  subjects with missing OA severity or KL measures
missing_severity <- df_2023 %>%
  filter(is.na(severity)) %>%        
  group_by(subject) %>%             # Ensure one row per subject
  summarise(
    severity = first(severity),  
    kl_score = first(kl_score),# Will be NA
    affected_knee = first(affected_knee), 
    signal_side = first(signal_side), 
    n = n()                         # Number of rows per subject (for debugging)
  ) %>%
  distinct(subject, .keep_all = TRUE)  # Keep unique subjects

print(missing_severity) # all missing severity, kl score, and affected knee data. Likely asymptomatic controls. See docs. Will be exlcuded. 

missing_kl <- df_2023 %>%
  filter(is.na(kl_score)) %>%        
  group_by(subject) %>%             # Ensure one row per subject
  summarise(
    severity = first(severity),  
    kl_score = first(kl_score),
    affected_knee = first(affected_knee), 
    signal_side = first(signal_side), 
    n = n()                         # Number of rows per subject (for debugging)
  ) %>%
  distinct(subject, .keep_all = TRUE) 

print(missing_kl)

df_2023 <- df_2023 %>%
  filter(!(subject %in% missing_severity$subject)) # exclude subjects with missing data
  #filter(!(subject %in% missing_kl$subject))

# Assign severity of contralateral knee
df_contra_2023 <- df_2023 %>%
  select(subject, signal_side, severity) %>%
  distinct(subject, signal_side, severity) %>%  # Ensure one row per subject-side
  mutate(signal_side = ifelse(signal_side == "LEFT", "RIGHT", "LEFT")) %>%
  rename(severity_contra = severity)

df_contra_kl_2023 <- df_2023 %>%
  select(subject, signal_side, kl_score) %>%
  distinct(subject, signal_side, kl_score) %>%  # Ensure one row per subject-side
  mutate(signal_side = ifelse(signal_side == "LEFT", "RIGHT", "LEFT")) %>%
  rename(kl_contra = kl_score)

# Merge on subject and signal_side only
df_2023 <- df_2023 %>%
  left_join(df_contra_2023, by = c("subject", "signal_side")) %>%
  left_join(df_contra_kl_2023, by = c("subject", "signal_side"))

print(head(df_2023 %>% 
             select(subject, signal_side, severity, severity_contra, kl_score, kl_contra)))

# Save final data frame 
save(df_2023, file = "data/df_RC_MWF_merged_2023.RData")
