rm(list = setdiff(ls(), c("script_files", "scripts_dir", "log_file", "script_path")))

library(fst)
library(tidyr)
library(dplyr)
library(here)
library(ggplot2)
library(stringr)

# Purpose: This script reads in the v3d data, filters it for relevant subjects, and calculates the mean wave forms and walking speed for each subject. The script returns a processed data frame. 

# Load data 
df_v3d_2025 <- read.fst("data/2025-03-27-140338_df_hmrl_v3d.fst")
df_v3d_2025 <- as.data.frame(df_v3d_2025)

# Ensure only 1 trial per encounter id and only the subjects first visit is present
## Filter for only first visit
df_v3d_2025 <- df_v3d_2025 %>%
  filter(action =='WALK') %>% #Filter for Walk only 
  filter(str_split(subject, "-", simplify = TRUE)[,3] == "1") #selects only ids with "1" in the third section of the id code (i.e. visit number)

## Check for subject IDs with >1 trials (duplicate encounter id entries)
multiple_trials <- df_v3d_2025 %>%
  filter(trial %in% c("00001", "00002")) %>%
  distinct(subject, trial) %>%
  group_by(subject) %>%
  summarise(count = n()) %>%
  filter(count > 1) %>%
  pull(subject)

print(multiple_trials) # Look up these encounter ids in Redcap to see which trial to use.  


# Walking Speeds
df_Metric_2025 <- df_v3d_2025 %>%
  filter(signal_types =='METRIC') %>%
  filter(signal_side =='OVERALL') %>%
  filter(signal_names == "SPEED_MEAN")

df_Metric_2025 %>%
  group_by(action) %>%
  summarise(Speed = mean(value),
            Speed_sd = sd(value))

df_Metric_2025 <- df_Metric_2025 %>%  
  rename("speed" = "value") %>%
  select(subject, speed)

# Caclulate mean waveforms (MWF)
df_MWF_2025 <- df_v3d_2025 %>%
  filter(signal_types =='LINK_MODEL_BASED') %>%
  group_by(subject, action, signal_names, signal_side, signal_components, item) %>%
  summarise(value  = mean(value)) 

  # Merge walking speed data with meanwaveforms
df_MWF_2025 <- df_MWF_2025 %>%
  left_join(df_Metric_2025, by = c("subject"))

# Plot Waveforms 
## Rename components
component_names <- c(
  `X`="Sagittal  Plane",
  `Y`="Frontal Plane",
  `Z`="Transverse Plane")

joint_angles <- c("HIP_ANGLE","KNEE_ANGLE","ANKLE_ANGLE")

## Ensemble Group Mean waveforms  ALL SUBJECTS
df_MWF_2025 %>%
  filter(signal_names %in% joint_angles) %>%
  ggplot(aes(x = item, y = value, group = action)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.25) +
  stat_summary(fun = mean, geom = "line",  linewidth=1)+
  stat_summary(
    fun.data = "mean_sdl",
    fun.args = list(mult = 1),
    #mapping = aes(color = action, fill = action),
    geom = "ribbon",
    alpha = 0.25)+
  facet_wrap(signal_names ~ signal_components, scales = "free",labeller = labeller(signal_components = as_labeller(component_names)), ncol = 3) +
  theme_minimal() +
  theme(axis.line = element_line(linewidth = 1, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        legend.position = 'right') +
  scale_x_continuous(expand = c(0, 0))+
  scale_color_manual(values=c("Walk" = "#b71313"))+
  xlab('Percent Gait Cycle (%)') +
  ylab('Joint Angle (Degrees)')


# Save Final data frame 
save(df_MWF_2025, file = "data/df_v3d_2025_preprocessed.RData")

