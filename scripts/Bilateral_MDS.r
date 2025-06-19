# MDS for clustering 

## Load Libraries 
library(dplyr)
library(ggplot2)
library(cluster)  # For dist() function with Manhattan distance
library(stats)

# Load in the data 
load("data/PCA_frontal_plane_results.RData")
load("data/PCA_sagittal_plane_results.RData")

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
  dplyr::select(subject, sex, age, speed, oa_location) %>%
  group_by(subject) %>%
  summarise(
    age = first(age),       # Keep first occurrence (all are the same)
    sex = first(sex),       # Keep first occurrence
    speed = first(speed),       # Keep first occurrence
  ) %>%
  ungroup()

# Select other relevant columns for MDS 
frontal_data <- frontal_pca_results %>%
  dplyr::select(subject, signal_side, PC1_front, PC2_front, PC3_front)
sagittal_data <- sagittal_pca_results %>%
  dplyr::select(subject, signal_side, PC1_sag, PC2_sag, PC3_sag)

# Merge dfs to inlcude frontal and sag plane PCs and walking speed and demographics 
combined_data <- frontal_data %>%
  left_join(sagittal_data, by = c("subject", "signal_side")) %>%
  inner_join(demographic_data, by = c("subject"))

# Verify dimensions 
dim(combined_data)  # Should be ~179
colnames(combined_data)  # Should include PC1_front, PC2_front, PC1_sag, PC2_sag, speed, age, sex etc.

# Select numeric variables for MDS
mds_data <- combined_data %>%
  dplyr::select(PC1_front, PC2_front, PC1_sag, PC2_sag, speed, age) %>%
  scale() %>%  # Standardize
  as.data.frame()

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
  mutate(subject = combined_data$subject[complete.cases(mds_data)],
         signal_side = combined_data$signal_side[complete.cases(mds_data)])

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
    .groups = "drop"
  )

dim(sev_data)

mds_sev <- mds_coords %>%
  inner_join(sev_data, by = c("subject", "signal_side"), relationship = "one-to-one")

# Verify
dim(mds_sev)  # Should match dim(mds_coords), 179 subjects 
colnames(mds_sev)

ggplot(mds_sev, aes(x = Dim1, y = Dim2, color = factor(severity_contra))) +
  geom_point(size = 2) +
  labs(title = "MDS with Manhattan Distance", x = "Dimension 1", y = "Dimension 2") +
  theme_minimal()

