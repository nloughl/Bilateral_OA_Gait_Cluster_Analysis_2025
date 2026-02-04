# MDS for clustering 

## Load Libraries 
library(dplyr)
library(ggplot2)
library(cluster)  # For dist() function with Manhattan distance
library(stats)

# Load in the data 
load("data/df_2023_2025_combined.RData") # for demographics 
load("data/PCA_frontal_plane_results.RData")
load("data/PCA_sagittal_plane_results.RData")

# load("data/PCA_frontal_plane_results_no_isopatell.RData") # If you want to exclude iso patellofemoral OA
# load("data/PCA_sagittal_plane_results_no_isopatell.RData") # If you want to exclude iso patellofemoral OA

## Prepare the data (Need to combine all numeric data to be clustered into 1 table)
# Rename frontal plane PCA columns
frontal_pca_results <- KneeFrontalPlanePCAll %>%
  rename(
    PC1_front = PC1,
    PC2_front = PC2,
    PC3_front = PC3
  )

# Rename sagittal plane PCA columns
sagittal_pca_results <- KneeSagPlanePCAll %>%
  rename(
    PC1_sag = PC1,
    PC2_sag = PC2,
    PC3_sag = PC3
  )

# Verify
colnames(frontal_pca_results)  # Should include PC1_front, PC2_front, PC3_front
colnames(sagittal_pca_results)  # Should include PC1_sag, PC2_sag, PC3_sag

# Get Demographic data (age and sex)
demographic_data <- df_combined %>%
  dplyr::select(subject, sex, age, speed) %>%
  group_by(subject) %>%
  summarise(
    age = first(age),       # Keep first occurrence (all are the same)
    sex = first(sex),       
    speed = first(speed)
  ) %>%
  ungroup()

# Select other relevant columns for MDS 
frontal_data <- frontal_pca_results %>%
  dplyr::select(subject, signal_side, PC1_front, PC2_front, PC3_front)
sagittal_data <- sagittal_pca_results %>%
  dplyr::select(subject, signal_side, PC1_sag, PC2_sag, PC3_sag)

# Merge dfs to inlcude frontal and sag plane PCs and walking speed and demographics 
mds_combined_data <- frontal_data %>%
  left_join(sagittal_data, by = c("subject", "signal_side")) %>%
  inner_join(demographic_data, by = c("subject"))

# Verify dimensions 
dim(mds_combined_data)  # Should be ~179 (or 145 if lateral cases removed)
colnames(mds_combined_data)  # Should include PC1_front, PC2_front, PC1_sag, PC2_sag, speed, age, sex etc.

# Identify PCS with >=80% culmulative proportion to include for analysis 
lines <- readLines("outputs/pca_summary_frontal.txt")
cumulative_line <- grep("Cumulative Proportion", lines, value = TRUE)
cum_props <- as.numeric(stringr::str_extract_all(cumulative_line, "\\d+\\.\\d+", simplify = TRUE))
num_pcs_to_keep <- which(cum_props >= 0.8)[1]
front_pcs <- paste0("PC", 1:num_pcs_to_keep, "_front")
print(front_pcs)

lines_sag <- readLines("outputs/pca_summary_sagittal.txt")
cumulative_line_sag <- grep("Cumulative Proportion", lines_sag, value = TRUE)
cum_props_sag <- as.numeric(stringr::str_extract_all(cumulative_line_sag, "\\d+\\.\\d+", simplify = TRUE))
num_pcs_to_keep_sag <- which(cum_props_sag >= 0.8)[1]
sag_pcs <- paste0("PC", 1:num_pcs_to_keep_sag, "_sag")
print(sag_pcs)

# Select numeric variables for MDS
mds_data <- mds_combined_data %>%
  dplyr::select(all_of(c(front_pcs, sag_pcs, "speed", "age"))) %>%
  scale() %>%  # Standardize
  as.data.frame() 

# Save mds input data with subject info 
mds_input_data <- mds_combined_data %>%
  dplyr::select(all_of(c(front_pcs, sag_pcs, "speed", "age"))) %>%
  mutate(subject = mds_combined_data$subject[complete.cases(mds_data)],
         signal_side = mds_combined_data$signal_side[complete.cases(mds_data)])
save(mds_input_data, file="data/mds_input_data.Rdata")

# # Handle any NAs (if any)
# mds_data <- na.omit(mds_data)


# Compute Manhattan distance matrix
dist_matrix <- dist(mds_data, method = "manhattan")

# Perform classical MDS
mds_result <- cmdscale(dist_matrix, k = 2, eig = TRUE)  # k = 2 for 2D visualization

df_mds_result <- as.data.frame(mds_result$points)
save(df_mds_result, file="data/df_mds_result.RData")

# Create data frame with MDS coordinates
mds_coords <- as.data.frame(mds_result$points) %>%
  rename(Dim1 = V1, Dim2 = V2) %>%
  mutate(subject = mds_combined_data$subject[complete.cases(mds_data)],
         signal_side = mds_combined_data$signal_side[complete.cases(mds_data)])

# Add other subject id and sex for context
mds_coords <- mds_coords %>%
  left_join(dplyr::select(demographic_data, subject, sex), by = "subject")

save(mds_coords, file="data/df_mds_coords.RData")

# Scatter plot of MDS results
ggplot(mds_coords, aes(x = Dim1, y = Dim2)) +
  geom_point(size = 2, color = "blue") +
  theme_minimal() +
  labs(title = "MDS of Gait Patterns (Frontal/Sagittal PCs, Walking Speed) with Manhattan Distance",
       x = "MDS Dimension 1", y = "MDS Dimension 2") +
  theme(legend.position = "right")


# Plot MDS results by severity_contra 
sev_data <- df_combined %>%
  group_by(subject, signal_side) %>%
  summarise(
    severity_contra = first(severity_contra),  # Take first if consistent
    kl_contra = first(kl_contra),
    .groups = "drop"
  )

dim(sev_data)

mds_sev <- mds_coords %>%
  inner_join(sev_data, by = c("subject", "signal_side"), relationship = "one-to-one")

# Verify
dim(mds_sev)  # Should match dim(mds_coords), 179 or 145 subjects 
colnames(mds_sev)

ggplot(mds_sev, aes(x = Dim1, y = Dim2, color = factor(severity_contra))) +
  geom_point(size = 2) +
  labs(title = "MDS with Manhattan Distance", x = "Dimension 1", y = "Dimension 2") +
  theme_minimal()

# Colour by kl_contra
ggplot(mds_sev, aes(x = Dim1, y = Dim2, color = kl_contra)) +
  geom_point(size = 2) +
  labs(title = "MDS with Manhattan Distance", x = "Dimension 1", y = "Dimension 2") +
  theme_minimal()

