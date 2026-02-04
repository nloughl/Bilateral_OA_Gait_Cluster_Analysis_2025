# Filter for at least one 'severe' leg 2025 
rm(list = setdiff(ls(), c("script_files", "scripts_dir", "log_file", "script_path")))

library(fst)
library(tidyr)
library(dplyr)
library(here)
library(ggplot2)
library(stringr)

load("data/df_merged_RC_v3d_2025.RData")

# Filter subjects with at least one 'severe' knee 
df_one_leg_2025 <- df_2025 %>%
  filter(!is.na(severity)) %>%  # Remove cases with NA
  #filter(item %in% (1:60)) %>%  # Stance phase only 
  group_by(subject) %>%
  filter(
    any(severity == "Severe")  # At least one leg is severe
  ) %>%
  ungroup()

print(length(unique(df_one_leg_2025$subject)))


# Save final data frame 
save(df_one_leg_2025, file="data/df_one_leg_2025.RData")
