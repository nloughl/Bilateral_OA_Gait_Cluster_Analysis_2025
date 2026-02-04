rm(list = ls())

#load libraries 
library(here)
library(tidyverse)
library(eatTools) #recodeLookup is helpful here.
library(stringr)
library(fst)
library(Hmisc)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork) # needed for plots
library(factoextra) # needed for fviz_eig (nice scree plots)


# Load the R workspace files
load("data/df_2023-08-23-153421_ORTHO_offline.RData")
load("data/df_RC_offline_2025-03-13.RData")


# Check the loaded objects
ls()

################################################################################
#Plot raw data first!

#Plot completeness of the dataset (check for missing values)
library(naniar) #CRAN library for feature completeness plotting
gg_miss_var(df_RC, facet=sex, show_pct = TRUE) + labs(title = "Missing data from j=952 observations for males (1) and females (2)")
gg_miss_var(df_RC, facet=severity, show_pct = TRUE) + labs(title = "Missing data from j=952 observations by OA severity)")

#Plot sex distribution 
ggplot(df_RC, aes(x = sex, fill = sex)) +
  geom_bar() +
  scale_x_discrete(labels = c("Male", "Female")) +  # Label axis properly
  labs(
    title = "Distribution of Sex in Dataset",
    x = "Sex",
    y = "Count"
  ) +
  theme_minimal()

#Plot severity/KL grade distribution 
ggplot(df_RC, aes(x = factor(kl_obs2), fill = factor(kl_obs2))) +
  geom_bar() +
  scale_fill_brewer(palette = "Set2") +  # Nice color scheme
  labs(
    title = "Distribution of KL Grades in Dataset",
    x = "KL Grade",
    y = "Count",
    fill = "KL Grade"
  ) +
  theme_minimal()

#Severity by Sex 
ggplot(df_RC, aes(x = factor(kl_obs2), fill = factor(kl_obs2))) +
  geom_bar() +
  facet_wrap(~ sex) +  # Split by sex
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "KL Grade Distribution by Sex",
    x = "KL Grade",
    y = "Count"
  ) +
  theme_minimal()

#Calculate descriptive statistics on demographics 
summary(df_RC)

################################################################################

## Calculate Means for Waveforms
dfMeanWaveforms <- df_v3d_in %>%
  filter(signal_types =='LINK_MODEL_BASED') %>%
  group_by(subject, action, signal_names, signal_side, signal_components, item) %>%
  summarise(value  = mean(value)) %>%
  filter(action =='Walk') #Filter for Walk only 

dfMeanWaveforms <-as.data.frame(dfMeanWaveforms)

# unfactor item (i.e. percent gait cycle in this case for waveforms)
#dfMeanWaveforms$item <- as.numeric(levels(dfMeanWaveforms$item))[dfMeanWaveforms$item]


## Rename components
component_names <- c(
  `X`="Sagittal  Plane",
  `Y`="Frontal Plane",
  `Z`="Transverse Plane")

joint_angles <- c("HIP_ANGLE","KNEE_ANGLE","ANKLE_ANGLE")

## Ensemble Group Mean waveforms  ALL SUBJECTS
dfMeanWaveforms %>%
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


# Walking Speeds----
  ## Calculate Means for Waveforms
  dfMetric <- df_v3d_in %>%
    filter(signal_types =='METRIC')
  
  # length(unique(side_mergeLookup$subject))
  
  
  #Filtered on Walk only
  dfMetric <- dfMetric %>%
    filter(signal_side =='OVERALL') %>%
    filter(signal_names == "SPEED_MEAN") %>%
    filter(action =='Walk')
  
  dfMetric %>%
    group_by(action) %>%
    summarise(Speed = mean(value),
              Speed_sd = sd(value))
  
  dfMetric<- dfMetric %>%  
    rename("speed" = "value")

################################################################################
#Filter RedCap data for columns we want 
  
  #Filter for relevant columns 
  df_RC_filtered <- df_RC %>%
    select(signal_side, subject, clothing_bottom, affected_knee, 
           age, sex, kl_obs2, oa_location_obs2, severity, None) %>% #Using 'severity' as severity measure 
    filter(clothing_bottom != "Skirt/dress below knees") #%>% #Filter out entries with long skirts

#################################################################################

## Merge Mean Waveforms with demographic data
  
dfMeanWaveforms_RC <- merge(dfMeanWaveforms, df_RC_filtered, by=c('subject','signal_side'), all.x = TRUE)  
  
  #Filter for only first visits   
  library(stringr)
  
  length(unique(dfMeanWaveforms_RC$subject))
  
  dfMeanWaveforms_RC <- dfMeanWaveforms_RC %>%
    filter(str_split(subject, "-", simplify = TRUE)[,3] == "1") #selects only ids with "1" in the third section of the id code (i.e. visit number)

  
  length(unique(dfMeanWaveforms_RC$subject))
  
  ##########################################
## COmpare dataset subjects unique 
  unique_subjects_2023 <- df_RC %>%
    filter(str_split(subject, "-", simplify = TRUE)[,3] == "1") #selects only ids with "1" in the third section of the id code (i.e. visit number)
  unique_subjects_2023 <- unique(unique_subjects_2023$subject)
  unique_subjects_2025 <- data_RC %>%
    filter(str_split(subject, "-", simplify = TRUE)[,3] == "1") #selects only ids with "1" in the third section of the id code (i.e. visit number)
  unique_subjects_2025 <- unique(unique_subjects_2025$subject)
  
  common_subjects <- intersect(unique_subjects_2023, unique_subjects_2025)
 
  
  
  # Assuming prior_subjects and recent_subjects are the vectors of subject IDs
  
  missing_subjects <- setdiff(unique_subjects_2025, unique_subjects_2023)
  print(missing_subjects)
  
  # Check if all prior subjects exist in the recent dataset
  if (length(missing_subjects) == 0) {
    print("All prior subjects are present in the recent dataset.")
  } else {
    print(paste("Missing", length(missing_subjects), "subjects in the recent dataset."))
  }
  
  #####################################
  
  
  
  
# Filter out any patients with only 1 leg entry (need left and right)
table(dfMeanWaveforms_RC$signal_side) #Uneven left and right entries 

library(dplyr)

dfMeanWaveforms_RC <- dfMeanWaveforms_RC %>%
  group_by(subject) %>%                      # Group by subject ID
  filter(all(c("LEFT", "RIGHT") %in% signal_side)) %>%  # Keep only subjects with both sides
  ungroup()

table(dfMeanWaveforms_RC$signal_side) #equal!!

##Save dataframe 
save(dfMeanWaveforms_RC, file="data/df_2023_ortho_OFFLINE_RC_merged.RData")


################################################################################
# Plot filtered sample by severity
library(dplyr)
library(ggplot2)

# Summarize the number of unique subjects per severity group
df_summary_sev <- dfMeanWaveforms_RC %>%
  group_by(severity) %>%
  summarise(unique_subjects = n_distinct(subject))

# Plot the summarized Severity data
ggplot(df_summary_sev, aes(x = factor(severity), y = unique_subjects, fill = factor(severity))) +
  geom_col() +  # Use geom_col since we’ve pre-calculated the counts
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Number of Unique Subjects by severity",
    x = "Severity",
    y = "Number of Unique Subjects"
  ) +
  theme_minimal()

# Summarize the number of unique subjects per severity group in 2025 dataset 
df_summary_sev_2025 <- df_severity_recoded %>%
  group_by(severity) %>%
  summarise(unique_subjects = n_distinct(subject))

# Plot the summarized Severity data
ggplot(df_summary_sev, aes(x = factor(severity), y = unique_subjects, fill = factor(severity))) +
  geom_col() +  # Use geom_col since we’ve pre-calculated the counts
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Number of Unique Subjects by severity",
    x = "Severity",
    y = "Number of Unique Subjects"
  ) +
  theme_minimal()


# Plot mean waveforms by Severity
plot_WFbySev <- dfMeanWaveforms_RC %>%
  filter(signal_names %in% joint_angles) %>%  
  filter(item %in% (1:60)) %>%
  ggplot(aes(x = item, y = value, color = factor(severity), fill = factor(severity))) +  
  geom_hline(yintercept = 0, color = "black", linewidth = 0.25) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  # stat_summary(
  #   fun.data = "mean_sdl",
  #   fun.args = list(mult = 1),
  #   geom = "ribbon",
  #   alpha = 0.25) +  # Fill will match color
  facet_wrap(signal_names ~ signal_components, 
             scales = "free", 
             labeller = labeller(signal_components = as_labeller(component_names)), 
             ncol = 3) +  # Assuming component_names is defined
  theme_minimal() +
  theme(
    axis.line = element_line(linewidth = 1, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(),
    legend.position = 'right'
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Mean Waveforms by Severity",
       x = "Gait Cycle (%)",
       y = "Mean Waveform Value",
       color = "Severity",
       fill = "Severity") +  # Update labels for kl_obs2
  scale_color_viridis(discrete = TRUE, option = "D") +  # Viridis for lines
  scale_fill_viridis(discrete = TRUE, option = "D")

plot_WFbySev

####################################################################################
## Check missingness in Severity

library(dplyr)

# Identify subjects with missing severity
missing_severity <- dfMeanWaveforms_RC %>%
  filter(is.na(severity)) %>%        # Focus on rows where kl_obs2 is NA
  group_by(subject) %>%             # Ensure one row per subject
  summarise(
    severity = first(severity),  
    kl_obs2 = first(kl_obs2),# Will be NA
    affected_knee = first(affected_knee), #Include this column 
    signal_side = first(signal_side), #Include this column # Include severity column
    n = n()                         # Number of rows per subject (for debugging)
  ) %>%
  distinct(subject, .keep_all = TRUE)  # Keep unique subjects

# Count the number of subjects with missing kl_obs2
n_missing_severity <- nrow(missing_severity)
print(paste("Number of unique subjects with missing severity:", n_missing_severity))

# View the results
print(missing_severity)

#########################################################################################
## Assign contralateral severity
library(dplyr)

#Assign severity of contralteral knee
df_SevContra <- dfMeanWaveforms_RC %>%
  select(subject, signal_side, severity) %>%
  distinct(subject, signal_side, severity) %>%  # Ensure one row per subject-side
  mutate(signal_side = ifelse(signal_side == "LEFT", "RIGHT", "LEFT")) %>%
  rename(severity_contra = severity)

# Merge on subject and signal_side only
dfMeanWaveforms_RC_sev <- dfMeanWaveforms_RC %>%
  left_join(df_SevContra, by = c("subject", "signal_side"))


# Verify dimensions
dim(dfMeanWaveforms_RC)  # Original size
dim(dfMeanWaveforms_RC_sev)

# Check the result
print(head(dfMeanWaveforms_RC_sev %>% 
             select(subject, signal_side, severity, severity_contra)))

# Check replaced knees 
# print(dfMeanWaveforms_RC_sev %>% 
#         filter(severity == "Replaced") %>% 
#         select(subject, signal_side, kl_obs2, affected_knee, severity, severity_contra))



# Count severe cases 
n_severe <- dfMeanWaveforms_RC_sev %>%
  filter(severity == "Severe") %>%
  distinct(subject) %>%
  nrow()
print(paste("Number of unique subjects with at least one severe knee:", n_severe))

# #Count moderate cases
# n_mod <- dfMeanWaveforms_RC_sev %>%
#   filter(severity == "Moderate") %>%
#   distinct(subject) %>%
#   nrow()
# print(paste("Number of unique subjects with at least one moderate knee:", n_mod))


### Filter data --------------------------------------------------------------------
# Filter subjects with at least one severe knee 
df_filtered <- dfMeanWaveforms_RC_sev %>%
  filter(!is.na(severity)) %>%  # Remove cases with NA
  filter(item %in% (1:60)) %>%  # Stance phase only 
  group_by(subject) %>%
  filter(
    any(severity == "Severe")  # At least one leg is severe
  ) %>%
  ungroup()

########## COLNAMES PREAPRE TO MERGE WITH 2025 DATA #######
colnames<- c("subject", "item", "sex", "severity_contra", "signal_side", "value", "severity", "age", )
#############################################################

# Summarize subjects with 2 Severe legs
severe_legs_analysis <- df_filtered %>%
  # Get unique severity and kl_obs2 per subject and signal_side
  distinct(subject, signal_side, severity, kl_obs2) %>%
  # Filter to only Severe legs
  filter(severity == "Severe") %>%
  group_by(subject) %>%
  summarise(
    n_severe_legs = n_distinct(signal_side),  # Count unique Severe legs
    n_unique_kl = n_distinct(kl_obs2, na.rm = TRUE),  # Count unique non-NA kl_obs2 values
    all_na_kl = all(is.na(kl_obs2)),  # Are all kl_obs2 NA?
    kl_values = paste(unique(kl_obs2), collapse = ", ")  # List kl_obs2 values
  ) %>%
  filter(n_severe_legs == 2)  # Only subjects with 2 Severe legs

# Cases with different kl_obs2
severe_kl_diff <- severe_legs_analysis %>%
  filter(n_unique_kl > 1)
print("Subjects with 2 Severe legs and different kl_obs2:")
print(severe_kl_diff)
print(paste("Number of subjects with 2 Severe legs and different kl_obs2:", nrow(severe_kl_diff)))

# Cases with all kl_obs2 missing
severe_kl_missing <- severe_legs_analysis %>%
  filter(all_na_kl)
print("Subjects with 2 Severe legs and all kl_obs2 Ascertainmentally missing kl_obs2:")
print(severe_kl_missing)
print(paste("Number of subjects with 2 Severe legs and all kl_obs2 missing:", nrow(severe_kl_missing)))

# Total subjects with 2 Severe legs for context
print(paste("Total subjects with 2 Severe legs:", nrow(severe_legs_analysis)))


# Select one leg per subject based on severity, then kl_obs2
set.seed(123)  # For reproducibility
df_one_leg <- df_filtered %>%
  group_by(subject) %>%
  mutate(
    severe_sides = sum(unique(signal_side[severity == "Severe"]) %in% c("LEFT", "RIGHT"), na.rm = TRUE),
    max_kl = max(kl_obs2[severity == "Severe"], na.rm = TRUE, default = -Inf),  # Default if all NA
    n_max_kl = sum(kl_obs2 == max_kl & severity == "Severe", na.rm = TRUE),
    all_na_kl = all(is.na(kl_obs2[severity == "Severe"])),
    # Ensure chosen_side is always valid
    chosen_side = case_when(
      severe_sides == 1 ~ first(signal_side[severity == "Severe"], default = first(signal_side)),  # One Severe leg
      severe_sides == 2 & (all_na_kl | n_max_kl == sum(severity == "Severe")) ~ sample(unique(signal_side[severity == "Severe"]), 1),  # Both Severe, no KL or equal
      severe_sides == 2 ~ first(signal_side[kl_obs2 == max_kl & severity == "Severe"], default = first(signal_side[severity == "Severe"])),  # Both Severe, max KL
      severe_sides == 0 ~ first(signal_side[severity == "Severe"], default = first(signal_side)),  # Fallback if severe_sides misfires
      TRUE ~ first(signal_side[severity == "Severe"], default = first(signal_side))  # Catch-all
    )
  ) %>%
  filter(!is.na(chosen_side) & signal_side == chosen_side) %>%  # Keep only chosen leg, handle NA
  select(-c(severe_sides, max_kl, n_max_kl, all_na_kl, chosen_side)) %>%
  ungroup()

# Verify/Ensure no samples lost 
print(paste("Unique subjects in df_filtered:", length(unique(df_filtered$subject))))
print(paste("Unique subjects in df_one_leg:", length(unique(df_one_leg$subject))))

## Summarize final sample 
# Summary table 
library(dplyr)
library(tableone)
library(pander)

sample_summary <- df_one_leg %>%
group_by(subject) 

variable_list <- c("age",  "sex", "severity_contra", "kl_obs2", "oa_location_obs2", "signal_side")
factor_variables <- c("sex", "severity_contra", "kl_obs2", "oa_location_obs2", "signal_side")
continuous_vars <- c("age")
 
dat_sample_summary <- df_one_leg %>%
      group_by(subject) %>%
      summarise(age = mean(age, na.rm = TRUE))
    
print(dat_sample_summary)

tab0 <- CreateTableOne(vars = variable_list,, data = dat_sample_summary, factorVars = factor_variables, test = FALSE)
print(tab0,  formatOptions = list(big.mark = ","))
 
tab0_matrix <- print(tab0,
  includeNA = TRUE,
  showAllLevels = TRUE
)
pandoc_table <- pandoc.table(tab0_matrix,
  split.table = Inf,
  style = "rmarkdown",
  caption = "Overall"
)
print(sample_summary)


# Summarize the number of unique subjects per contra severity group
df_one_leg_plot <- df_one_leg %>%
  group_by(severity_contra) %>%
  summarise(unique_subjects = n_distinct(subject))


# Plot the summarized Severity data
ggplot(df_one_leg_plot, aes(x = factor(severity_contra), y = unique_subjects, fill = factor(severity_contra))) +
  geom_col() +  # Use geom_col since we’ve pre-calculated the counts
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Number of Unique Subjects by Contralateral severity",
    x = "Severity",
    y = "Number of Unique Subjects"
  ) +
  theme_minimal()



### Plot final sample data 
library(ggplot2)
library(dplyr)
library(viridis)

# Plot mean waveforms colored by severity_contra (stance only)
plot_WFbySev_contra_stance <- df_one_leg %>%
  filter(signal_names %in% joint_angles) %>%
  filter(item %in% (1:60)) %>%
  ggplot(aes(x = item, y = value, color = severity_contra, fill = severity_contra)) +  # Map both color and fill
  geom_hline(yintercept = 0, color = "black", linewidth = 0.25) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(
    fun.data = "mean_sdl",
    fun.args = list(mult = 1),
    geom = "ribbon",
    alpha = 0.25) +  # Fill will match color
  facet_wrap(signal_names ~ signal_components, 
             scales = "free", 
             labeller = labeller(signal_components = as_labeller(component_names)), 
             ncol = 3) +
  theme_minimal() +
  theme(
    axis.line = element_line(linewidth = 1, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(),
    legend.position = 'right'
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Mean Stance Phase Waveforms Colored by Contralateral Severity",
       x = "Gait Cycle (%)",
       y = "Mean Waveform Value",
       color = "Contralateral Severity",
       fill = "Contralateral Severity") +  # Ensure fill is also labeled
  scale_color_viridis(discrete = TRUE, option = "D") +  # Viridis for lines
  scale_fill_viridis(discrete = TRUE, option = "D")

plot_WFbySev_contra_stance

# Plot mean waveforms colored by kl_contralateral (stance only) -- Power point Size 
plot1 <- df_one_leg %>%
  filter(signal_names %in% joint_angles) %>%
  ggplot(aes(x = item, y = value, color = severity_contra, fill = severity_contra)) +  # Map both color and fill
  geom_hline(yintercept = 0, color = "black", linewidth = 0.25) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  # stat_summary(
  #   fun.data = "mean_sdl",
  #   fun.args = list(mult = 1),
  #   geom = "ribbon",
  #   alpha = 0.25) +  # Fill will match color
  facet_wrap(signal_names ~ signal_components, 
             scales = "free", 
             labeller = labeller(signal_components = as_labeller(component_names)), 
             ncol = 3) +
  theme_minimal() +
  theme(
    axis.line = element_line(linewidth = 1, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), 
    panel.background = element_blank(),
    text = element_text(size = 20),
    legend.position = 'right'
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Mean Waveforms by Contralateral Severity",
       x = "Gait Cycle (%)",
       y = "Mean Waveform Value",
       color = "Contralateral Severity",
       fill = "Contralateral Severity") +  # Ensure fill is also labeled
  scale_color_viridis(discrete = TRUE, option = "D") +  # Viridis for lines
  scale_fill_viridis(discrete = TRUE, option = "D")
 
plot1

# Save the plot with ggsave
ggsave("plot_WFbySev_contra_stance_PP.png", plot = plot1, width = 12, height = 8, dpi = 300)


## Summarize Sample Characteristics
library(dplyr)
library(tableone)
library(pander)

variable_list <- c("sex", "kl_obs2", "severity_contra", "age", "oa_location_obs2")
factor_variables <- c("sex", "kl_obs2", "severity_contra", "oa_location_obs2")

# Sumamrize Characteristics 
df_sample_summary <- df_one_leg %>%
  group_by(subject) %>%
  summarise(
    age = first(age),       # Keep first occurrence (all are the same)
    sex = first(sex),       # Keep first occurrence
    kl_obs2 = first(kl_obs2),  # Keep first occurrence
    severity_contra = first(severity_contra),
    oa_location_obs2 = first(oa_location_obs2)
  ) %>%
  ungroup()

# Create Table 
tbl_sample_summary <- CreateTableOne(vars = variable_list,, data = df_sample_summary, factorVars = factor_variables, test = FALSE)
print(tbl_sample_summary,  formatOptions = list(big.mark = ","))

tbl_sample_summary_matrix <- print(tbl_sample_summary,
                     includeNA = TRUE,
                     showAllLevels = TRUE
)
pandoc_table <- pandoc.table(tbl_sample_summary_matrix,
                             split.table = Inf,
                             style = "rmarkdown",
                             caption = "Overall"
)
length(unique(df_one_leg$subject))

# Histogram of sex distribution within severity grades 

  ggplot(data=df_sample_summary, aes(x=severity_contra)) + theme_bw() + 
    facet_grid(. ~ sex) +
    geom_histogram(stat="count") + 
    labs(title = "Number of Observations for males and females by Contralateral Severity Rating", x="Contralateral Severity", y="Count")
