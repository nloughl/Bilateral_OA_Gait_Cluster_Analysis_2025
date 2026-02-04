rm(list = setdiff(ls(), c("script_files", "scripts_dir", "log_file", "script_path")))

library(fst)
library(tidyr)
library(dplyr)
library(here)

# Purpose: This script takes the patient OA severity information, currently rated using numeric KL scores (1-4) and converts them to descriptive strings (e.g. Severe, mild, moderate) to match the severity measure used in other data sets for this project. 

# Load data to be recoded 
df_RC_2025 <- read.fst("data/2025-03-27-140440_df_hmrl_severity.fst")

# View data and identify columns to be renamed/recoded 
#View(df_RC_2025)
colnames(df_RC_2025)

# Recode right and left knee severity from numeric KL to string description 
df_RC_2025 <- df_RC_2025 |>
  mutate(right_sev_str = recode(df_RC_2025$right_knee_severity,
         "1" = "Asymptomatic",
         "2" = "Mild",
         "3" = "Moderate",
         "4" = "Severe"
  )) |>
  
  mutate(left_sev_str = recode(df_RC_2025$left_knee_severity,
                                "1" = "Asymptomatic",
                                "2" = "Mild",
                                "3" = "Moderate",
                                "4" = "Severe"
  ))

# Ensure correct recoding 
#View(df_RC_2025)
table(df_RC_2025$left_sev_str)
table(df_RC_2025$right_sev_str)

# Check for replacements 
df_RC_2025 %>%
  filter(previous_replacements != "none") %>%
  distinct(subject, .keep_all = TRUE) %>%
  dplyr::select(subject, previous_replacements, right_sev_str, left_sev_str)

df_RC_2025 <- df_RC_2025 %>%
  mutate(right_sev_str = case_when(
    subject %in% c("E21-337-1-100", "E21-472-1-100") ~ "Replaced",  # Replace for specific subjects
    TRUE ~ right_sev_str  # Keep existing values for others
  )) %>%
  mutate(left_sev_str = case_when(
    subject == "E21-337-1-100" ~ "Replaced",
    TRUE ~ left_sev_str
  ))

# Contralateral Recoding (change right and left to ipsi and contra - lateral)
df_RC_2025 <- df_RC_2025 %>%
  mutate(
    severity = case_when(
      signal_side == "RIGHT" ~ right_sev_str,
      signal_side == "LEFT" ~ left_sev_str,
      TRUE ~ NA_character_
    ),
    severity_contra = case_when(
      signal_side == "RIGHT" ~ left_sev_str,
      signal_side == "LEFT" ~ right_sev_str,
      TRUE ~ NA_character_
    )
  )

df_RC_2025 <- df_RC_2025 %>%
  mutate(
    kl_score = case_when(
      signal_side == "RIGHT" ~ kl_right_o2,
      signal_side == "LEFT" ~ kl_left_o2,
      TRUE ~ NA_character_
    ),
    kl_contra = case_when(
      signal_side == "RIGHT" ~ kl_left_o2,
      signal_side == "LEFT" ~ kl_right_o2,
      TRUE ~ NA_character_
    )
  )

# Check the result
df_RC_2025 %>%
  filter(right_sev_str != left_sev_str) %>%
  dplyr::select(subject, signal_side, right_sev_str, left_sev_str, severity, severity_contra) %>%
  head(5)
df_RC_2025 %>%
  filter(kl_right_o2 != kl_left_o2) %>%
  dplyr::select(subject, signal_side, kl_right_o2, kl_left_o2, kl_score, kl_contra) %>%
  head(5)

# OA location Recoding
df_RC_2025 <- df_RC_2025 %>%
  mutate(
    oa_location = case_when(
      signal_side == "RIGHT" ~ right_knee_oa_location,
      signal_side == "LEFT" ~ left_knee_oa_location,
      TRUE ~ NA_character_
    )
  ) |>
  
  mutate(oa_location = recode(oa_location,
                               "patellofemoral" = "isolated patellofemoral",
                               "medial, lateral, patellofemoral" = "tricompartmental"
  ))

# Check result 
df_RC_2025 %>%
  dplyr::select(subject, signal_side, right_knee_oa_location, left_knee_oa_location, oa_location) %>%
  head(5)

# Sex Recoding
df_RC_2025 <- df_RC_2025 |>
  mutate(sex = recode(sex,
                                "male" = "Male",
                                "female" = "Female",
                      "prefer_not_to_disclose" = "Prefer not to disclose")
  )
print(unique(df_RC_2025$sex))

# Filter for Columns we want 
df_RC_2025 <- df_RC_2025 %>%
  dplyr::select(subject,signal_side, age, sex, oa_location, severity, severity_contra, kl_score, kl_contra) 

head(df_RC_2025)

# Save Output 
save(df_RC_2025, file="data/df_RC_2025_preprocessed.RData")

