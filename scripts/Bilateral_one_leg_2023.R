# Filter for at least one 'severe' leg (or also include healthy controls)
rm(list = setdiff(ls(), c("script_files", "scripts_dir", "log_file", "script_path")))

library(fst)
library(tidyr)
library(dplyr)
library(here)
library(ggplot2)
library(stringr)

load("data/df_RC_MWF_merged_2023.RData")

# Filter subjects with at least one 'severe' knee 
df_filtered_2023 <- df_2023 %>%
  filter(!is.na(severity)) %>%  # Remove cases with NA
  #filter(item %in% (1:60)) %>%  # Stance phase only 
  group_by(subject) %>%
  filter(
    any(severity == "Severe")  # At least one leg is severe
  ) %>%
  ungroup()

# Select one leg per subject based on severity, then kl_score
set.seed(123)  # For reproducibility
df_one_leg_2023 <- df_filtered_2023 %>%
  group_by(subject) %>%
  mutate(
    severe_sides = sum(unique(signal_side[severity == "Severe"]) %in% c("LEFT", "RIGHT"), na.rm = TRUE),
    max_kl = max(kl_score[severity == "Severe"], na.rm = TRUE, default = -Inf),  # Default if all NA
    n_max_kl = sum(kl_score == max_kl & severity == "Severe", na.rm = TRUE),
    all_na_kl = all(is.na(kl_score[severity == "Severe"])),
    # Ensure chosen_side is always valid
    chosen_side = case_when(
      severe_sides == 1 ~ first(signal_side[severity == "Severe"], default = first(signal_side)),  # One Severe leg
      severe_sides == 2 & (all_na_kl | n_max_kl == sum(severity == "Severe")) ~ sample(unique(signal_side[severity == "Severe"]), 1),  # Both Severe, no KL or equal
      severe_sides == 2 ~ first(signal_side[kl_score == max_kl & severity == "Severe"], default = first(signal_side[severity == "Severe"])),  # Both Severe, max KL
      severe_sides == 0 ~ first(signal_side[severity == "Severe"], default = first(signal_side)),  # Fallback if severe_sides misfires
      TRUE ~ first(signal_side[severity == "Severe"], default = first(signal_side))  # Catch-all
    )
  ) %>%
  filter(!is.na(chosen_side) & signal_side == chosen_side) %>%  # Keep only chosen leg, handle NA
  select(-c(severe_sides, max_kl, n_max_kl, all_na_kl, chosen_side)) %>%
  ungroup()

#-------------------------------------------------------------------------------
# ### Add healthy controls (comment out if not wanted) 
# asymp_controls <- df_2023 %>%
#   #filter(item %in% 1:60) %>%  #stance phase filter
#   group_by(subject) %>%
#   summarise(all_asymp = all(severity == "Asymptomatic")) %>%  # Check if all rows (both legs) are Asymptomatic
#   filter(all_asymp) %>%
#   pull(subject)  # Get subject IDs
# 
# ### Filter data & Select 1 leg at random
# set.seed(123)
# df_healthy <- df_2023 %>%
#   filter(subject %in% asymp_controls) %>%
#   distinct(subject, signal_side) %>%
#   group_by(subject) %>%
#   slice_sample(n = 1) %>%
#   ungroup()
# 
# ### Extract Data 
# df_healthy_one_leg <- df_2023 %>%
#   semi_join(df_healthy, by = c("subject", "signal_side"))
# 
# ### Merge into  severe OA data frame
# df_one_leg_2023_with_controls <- bind_rows(df_one_leg_2023, 
#                              df_healthy_one_leg %>% dplyr::select(any_of(names(df_one_leg_2023)))
#                              )
# ###Save df
# save(df_one_leg_2023_with_controls, file="data/df_one_leg_2023_w_controls.RData")
#-------------------------------------------------------------------------------

# Verify/Ensure no samples lost 
print(paste("Unique subjects in df_filtered:", length(unique(df_filtered_2023$subject))))
print(paste("Unique subjects in df_one_leg:", length(unique(df_one_leg_2023$subject))))

# Save final data frame 
save(df_one_leg_2023, file="data/df_one_leg_2023.RData")
