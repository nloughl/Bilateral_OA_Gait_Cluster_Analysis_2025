## Hierarchical Clustering (agglom), using Ward's minimum variance criteria
### Evaluate clustering tendency 
# If the value of Hopkins statistic is near 0.5, the data are uniformly distributed. If the value of the Hopkins statistic is far from 0.5, the data are clusterable


# Load data 
load("data/df_mds_coords.RData")
load("data/df_mds_result.RData")

#MDS results
res_mds <- get_clust_tendency(df_mds_result, n = nrow(mds_data)-1, graph = FALSE)
print(res_mds) 

### clustering (Using MDS results)
library(clusterSim)
library(cluster)     # For clustering and silhouette
library(factoextra)  # For visualization
library(FactoMineR)  # For HCPC clustering
library(dplyr)       # For data manipulation
library(ggplot2)     # For visualization
library(mclust)      # For ARI


# Compute Manhattan distance on MDS results 
dist_matrix_mds <- dist(df_mds_result, method = "manhattan")

# Perform hierarchical clustering using Ward's method
hc_ward <- hclust(dist_matrix_mds, method = "ward.D2")

#Determine number of clsuters 
  # Computes average silhouette of observations for different values of k
  fviz_nbclust(
    df_mds_result,                   
    FUNcluster = hcut,            # Hierarchical clustering 
    method = "silhouette",        
    k.max = 10,                   # Maximum number of clusters to test
    hc_method = "ward.D2",        # Linkage method 
  )

# Plot dendrogram
plot(hc_ward, main = "Dendrogram (Ward's Method)", sub = "", xlab = "", cex = 0.6)

#Dendrogram to inspect number of clusters 
# Cut dendrogram into k clusters (can edit number of clusters)
clusters <- cutree(hc_ward, k = 3)

# Add clusters to data
mds_coords$cluster <- factor(clusters)

# Plot clusters
ggplot(mds_coords, aes(x = Dim1, y = Dim2, color = cluster)) +
  geom_point(size = 3) +
  labs(title = "Hierarchical Clustering on MDS Coordinates") +
  theme_minimal()

# Compute silhouette scores
sil_scores <- silhouette(clusters, dist_matrix_mds)

# Plot silhouette scores
sil_width_plot <- fviz_silhouette(sil_scores)
ggsave("outputs/sil_width_plot.png", plot = sil_width_plot, width = 8, height = 6, dpi = 300)

### Outlier Detection 
conVars <- c("Dim1", "Dim2")#
IsOutlier <- function(x) {
  outlier <- (x < quantile(x, 0.25, na.rm = TRUE) - 3 * IQR(x, na.rm = TRUE)) |
    (x > quantile(x, 0.75, na.rm = TRUE) + 3 * IQR(x, na.rm = TRUE))
  as.numeric(outlier)
}

# Identify outliers for continuous variables
for (i in names(conVars)) {
  mds_coords[[paste0("Outlier_", i)]] <- IsOutlier(Data.Cluster[[i]])
}

# Create Exclude column
mds_coords$Exclude <- rowSums(mds_coords[, grepl("Outlier", names(mds_coords))], na.rm = TRUE)

# Check for outliers
has_outliers <- any(mds_coords$Exclude != 0)
print(has_outliers)


### Evaluate with adjusted Rand Index (ARI)

# Method from Young-Shand et al. 2023:
#"The cluster model of the full data set was compared with 100 subsets of the full data set (fraction of 0.8) for k = 2:10, and the mean ARI was reported. The selected clustering  solution (k clusters) with the greatest mean s and ARI was chosen."

library(dplyr)
library(cluster)  # For silhouette
library(mclust)   # For adjustedRandIndex
library(stats)    # For hclust

# Hierarchical clustering function
hclust_function <- function(data, k) {
  dist_matrix <- dist(data, method = "manhattan")
  hc <- hclust(dist_matrix, method = "ward.D2")
  clusters <- cutree(hc, k = k)
  return(clusters)
}

# Parameters
n_iterations <- 100  # Number of subsamples
fraction <- 0.8      # 80% of data
k_range <- 2:10      # Test k from 2 to 10
n_samples <- round(nrow(df_mds_result) * fraction)  # Number of rows to sample

# Store results
ari_results <- matrix(NA, nrow = n_iterations, ncol = length(k_range))
colnames(ari_results) <- paste("k =", k_range)
mean_silhouette <- numeric(length(k_range))  # Mean silhouette for full data

set.seed(123)  # For reproducibility
# Cluster full dataset for each k
full_clusters <- lapply(k_range, function(k) hclust_function(df_mds_result, k))

# Subsampling and ARI calculation
for (i in 1:n_iterations) {
  # Subsample 80% of the data
  sample_idx <- sample(1:nrow(df_mds_result), size = n_samples, replace = FALSE)
  sampled_data <- df_mds_result[sample_idx, ]
  
  # Cluster subsample for each k and compute ARI
  for (j in seq_along(k_range)) {
    k <- k_range[j]
    subsample_clusters <- hclust_function(sampled_data, k)
    full_clusters_sub <- full_clusters[[j]][sample_idx]  # Subset full clusters to match
    ari_results[i, j] <- adjustedRandIndex(full_clusters_sub, subsample_clusters)
  }
}

# Compute mean ARI for each k
mean_ari <- colMeans(ari_results, na.rm = TRUE)
names(mean_ari) <- paste("k =", k_range)

# Compute mean silhouette for full dataset
for (j in seq_along(k_range)) {
  k <- k_range[j]
  sil <- silhouette(full_clusters[[j]], dist(df_mds_result, method = "euclidean"))
  mean_silhouette[j] <- mean(sil[, 3], na.rm = TRUE)
}
names(mean_silhouette) <- paste("k =", k_range)

# Print results
cat("\nMean Adjusted Rand Index (Stability) for each k:\n")
print(round(mean_ari, 3))

cat("\nMean Silhouette Width for each k (Full Dataset):\n")
print(round(mean_silhouette, 3))

# Save results to csv 
subsampling_results <- data.frame(
  k = k_range,
  Mean_ARI = round(mean_ari, 3),
  Mean_Silhouette = round(mean_silhouette, 3)
)

write.csv(subsampling_results, file = "outputs/subsampling_results.csv", row.names = FALSE)


# Select optimal k (highest mean ARI and silhouette)
combined_score <- mean_ari + mean_silhouette  # Simple sum (could adjust weighting)
optimal_k_idx <- which.max(combined_score)
optimal_k <- k_range[optimal_k_idx]
cat("\nSelected k with highest combined mean ARI and silhouette:", optimal_k, "\n")
cat("Mean ARI for selected k:", round(mean_ari[optimal_k_idx], 3), "\n")
cat("Mean Silhouette for selected k:", round(mean_silhouette[optimal_k_idx], 3), "\n")

# Final clustering with optimal k
final_clusters <- hclust_function(df_mds_result, optimal_k)

# Add final clusters to data
mds_coords$final_clusters <- factor(final_clusters)
save(mds_coords, file="data/df_HC_mds_coords.RData")

## Create cluster plot based on HCPC
HC_MDS_plot <- fviz_cluster(
  list(data = mds_coords[, c("Dim1", "Dim2")], cluster = mds_coords$final_clusters),
  repel = TRUE,
  show.clust.cent = FALSE,
  show.legend.text = FALSE,
  pointsize = 1.8,
  geom = "point",
  ggtheme = theme_minimal() +
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)),
  main = "Hierarchical Clustering on MDS Coordinates"
) +
  theme_bw()
HC_MDS_plot
ggsave("outputs/HC_MDS_plot.png", plot = HC_MDS_plot, width = 8, height = 6, dpi = 300)

HC_MDS_plot_v2 <- ggplot(mds_coords, aes(x = Dim1, y = Dim2, color = final_clusters)) +
  geom_point(size = 3) +
  labs(title = "Hierarchical Clustering on MDS Coordinates") +
  theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
HC_MDS_plot_v2
ggsave("outputs/HC_MDS_plot_v2.png", plot = HC_MDS_plot_v2, width = 8, height = 6, dpi = 300)

####################################################################################3

### Plot Waveforms by cluster 
# Merge cluster assignment to waveform data 

dat_clust_plot <- df_combined %>%
  left_join(mds_coords[, c("subject", "final_clusters")], join_by("subject"))

component_names <- c(
  `X`="Sagittal  Plane",
  `Y`="Frontal Plane",
  `Z`="Transverse Plane")


# Plot
library(viridis)

plot_WFbyClust_mds_front <- dat_clust_plot %>%
  filter(signal_names == "KNEE_ANGLE") %>%
  filter(signal_components == "Y") %>%
  filter(item %in% (1:60)) %>%
  ggplot(aes(x = item, y = value, color = final_clusters, fill = final_clusters)) +  # Map both color and fill
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
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Mean Waveforms Coloured by Cluster",
       x = "Gait Cycle (%)",
       y = "Mean Waveform Value",
       color = "Cluster",
       fill = "Cluster") +  # Ensure fill is also labeled
  scale_color_viridis(discrete = TRUE, option = "D") +  # Viridis for lines
  scale_fill_viridis(discrete = TRUE, option = "D")

plot_WFbyClust_mds_front
ggsave("outputs/HC_WFbyClust_plot_front.png", plot = plot_WFbyClust_mds_front, width = 8, height = 6, dpi = 300)

plot_WFbyClust_mds_front_PP <- dat_clust_plot %>%
  filter(signal_names == "KNEE_ANGLE") %>%
  filter(signal_components == "Y") %>%
  filter(item %in% (1:60)) %>%
  ggplot(aes(x = item, y = value, color = final_clusters, fill = final_clusters)) +  # Map both color and fill
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
  labs(title = "Mean Waveforms Coloured by Cluster",
       x = "Gait Cycle (%)",
       y = "Mean Waveform Value",
       color = "Cluster",
       fill = "Cluster") +  # Ensure fill is also labeled
  scale_color_viridis(discrete = TRUE, option = "D") +  # Viridis for lines
  scale_fill_viridis(discrete = TRUE, option = "D")

plot_WFbyClust_mds_front_PP
ggsave("outputs/HC_WFbyClust_plot_front_PP.png", plot = plot_WFbyClust_mds_front_PP, width = 8, height = 6, dpi = 300)


plot_WFbyClust_mds_sag <- dat_clust_plot %>%
  filter(signal_names == "KNEE_ANGLE") %>%
  filter(signal_components == "X") %>%
  ggplot(aes(x = item, y = value, color = final_clusters, fill = final_clusters)) +  # Map both color and fill
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
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Mean Waveforms Coloured by Cluster",
       x = "Gait Cycle (%)",
       y = "Mean Waveform Value",
       color = "Cluster",
       fill = "Cluster") +  # Ensure fill is also labeled
  scale_color_viridis(discrete = TRUE, option = "D") +  # Viridis for lines
  scale_fill_viridis(discrete = TRUE, option = "D")

plot_WFbyClust_mds_sag
ggsave("outputs/HC_WFbyClust_plot_sag.png", plot = plot_WFbyClust_mds_sag, width = 8, height = 6, dpi = 300)

plot_WFbyClust_mds_sag_PP <- dat_clust_plot %>%
  filter(signal_names == "KNEE_ANGLE") %>%
  filter(signal_components == "X") %>%
  ggplot(aes(x = item, y = value, color = final_clusters, fill = final_clusters)) +  # Map both color and fill
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
  labs(title = "Mean Waveforms Coloured by Cluster",
       x = "Gait Cycle (%)",
       y = "Mean Waveform Value",
       color = "Cluster",
       fill = "Cluster") +  # Ensure fill is also labeled
  scale_color_viridis(discrete = TRUE, option = "D") +  # Viridis for lines
  scale_fill_viridis(discrete = TRUE, option = "D")

plot_WFbyClust_mds_sag_PP
ggsave("outputs/HC_WFbyClust_plot_sag_PP.png", plot = plot_WFbyClust_mds_sag_PP, width = 8, height = 6, dpi = 300)

############################################################################################3

# Evaluate Cluster Stability 
library(fpc)

kbest = 3 #can edit according to number of clusters chosen 
perc = 0.8
s = 987

# Run clusterboot for stability
cboot_HCPC = clusterboot(dist_matrix_mds, clustermethod = hclustCBI, method = "ward.D2", k = kbest, bootmethod = "boot",B=200, seed= s)

print(cboot_HCPC)
# Save results
# Extract components
jaccard <- cboot_HCPC$bootmean
dissolved <- cboot_HCPC$bootbrd
recovered <- cboot_HCPC$bootrecover

# Combine into a data frame
cboot_df <- data.frame(
  Cluster = seq_along(jaccard),
  Jaccard = jaccard,
  Dissolved = dissolved,
  Recovered = recovered
)

overall_mean_jaccard_HCPC <- mean(cboot_HCPC$bootmean, na.rm = TRUE)

jaccard_summary_HCPC <- rbind(
  cboot_df,
  data.frame(
    Cluster = "Overall_Mean",
    Jaccard = overall_mean_jaccard_HCPC,
    Dissolved = NA,
    Recovered = NA
  )
)
# Save to CSV
write.csv(jaccard_summary_HCPC, "outputs/jaccard_HCPC_with_overall.csv", row.names = FALSE)


## Cluster Characteristics
library(dplyr)
library(tableone)
library(pander)

variable_list <- c("age", "PC1_front", "PC1_sag", "speed", "sex", "severity_contra", "oa_location", "signal_side")
factor_variables <- c("sex", "severity_contra", "oa_location", "signal_side")
continuous_vars <- c("age", "PC1_front", "PC1_sag", "speed")
 
dat_clust_tables <- combined_data %>%
  left_join(mds_coords[, c("subject", "final_clusters")], join_by("subject"))
dat_clust_tables <- dat_clust_tables %>%
  left_join(
    df_combined %>%
      group_by(subject) %>%
      summarise(
        oa_location = first(oa_location),
        severity_contra = first(severity_contra)
      ) %>%
      dplyr::select(subject, oa_location, severity_contra), 
    by = "subject"
    
    )


tab0 <- CreateTableOne(vars = variable_list,, data = dat_clust_tables, factorVars = factor_variables, test = FALSE)
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

write.csv(tab0_matrix, file = "outputs/HCPC_sample_summary.csv")

## cluster tables
tab1 <- CreateTableOne(vars = variable_list, strata = "final_clusters" , data = dat_clust_tables, factorVars = factor_variables, test = FALSE)
print(tab1)

tab1_matrix <- print(tab1,
  includeNA = TRUE,
  showAllLevels = TRUE
)
pandoc_table <- pandoc.table(tab1_matrix,
  split.table = Inf,
  style = "rmarkdown",
  caption = "Overall"
)

#convert to .csv for report 
write.csv(tab1_matrix, file = "outputs/HCPC_cluster_summaries.csv")



## Statistical tests on Clusters
library(broom)
# Ensure cluster is a factor
dat_clust_tables <- dat_clust_tables %>% mutate(cluster = as.factor(final_clusters))

test_continuous_var <- function(data, variable, cluster_col = "final_clusters", file_prefix = NULL) {
  var_sym <- rlang::sym(variable)
  cluster_sym <- rlang::sym(cluster_col)
  
  # Normality check
  shapiro_test <- data %>%
    group_by(!!cluster_sym) %>%
    filter(n() >= 3) %>% # exclude clusters with less than 3 observations 
    summarise(p_value = shapiro.test(!!var_sym)$p.value, .groups = "drop")
  print(shapiro_test)
  
  # Check excluded clusters
  excluded <- data %>%
    group_by(!!cluster_sym) %>%
    filter(n() < 3) %>%
    summarise(n = n(), .groups = "drop")
  print(excluded)
  
  # Determine which test to run
  if (any(shapiro_test$p_value < 0.05)) {
    test_result <- kruskal.test(reformulate(cluster_col, response = variable), data = data)
    test_type <- "Kruskal-Wallis"
    test_p <- test_result$p.value
  } else {
    test_result <- aov(reformulate(cluster_col, response = variable), data = data)
    test_type <- "ANOVA"
    test_summary <- summary(test_result)
    test_p <- test_summary[[1]][["Pr(>F)"]][1]
  }
  
  cat(paste0("\n", test_type, " test for ", variable, ":\n"))
  print(test_result)
  
  # Save main test result (optional -- if signif, pairwise tests will be saved)
  if (!is.null(file_prefix)) {
    tidy_test <- tidy(test_result)
    write.csv(tidy_test, file = paste0("outputs/", file_prefix, "_", test_type, ".csv"), row.names = FALSE)
  }
  
  # Plot boxplot
  boxplot(as.formula(paste(variable, "~", cluster_col)), data = data,
          main = paste(variable, "Across Clusters"), xlab = "Cluster")
  
  # Pairwise test if significant
  if (test_p < 0.05) {
    if (inherits(test_result, "aov")) {
      tukey_result <- TukeyHSD(test_result)
      print(tukey_result)
      pairwise_df <- as.data.frame(tukey_result[[cluster_col]])
      pairwise_df$Comparison <- rownames(pairwise_df)
      pairwise_df <- pairwise_df[, c("Comparison", names(pairwise_df)[1:(ncol(pairwise_df)-1)])]
      method <- "TukeyHSD"
    } else {
      wilcox_result <- pairwise.wilcox.test(data[[variable]], data[[cluster_col]], p.adjust.method = "bonferroni")
      print(wilcox_result)
      pairwise_df <- tidy(wilcox_result)
      method <- "Wilcoxon"
    }
    
    # Save pairwise result
    if (!is.null(file_prefix)) {
      write.csv(pairwise_df, file = paste0("outputs/", file_prefix, "_", method, ".csv"), row.names = FALSE)
    }
  }
}

test_continuous_var(dat_clust_tables, "age", file_prefix = "HCPC_Age")
test_continuous_var(dat_clust_tables, "speed", file_prefix = "HCPC_Speed")
test_continuous_var(dat_clust_tables, "PC1_front", file_prefix = "HCPC_PC1Front")
test_continuous_var(dat_clust_tables, "PC1_sag", file_prefix = "HCPC_PC1Sag")


# Categorical Variables
test_categorical <- function(var) {
  # Create contingency table of unique subjects
  unique_data <- dat_clust_tables %>%
    distinct(subject, final_clusters, !!sym(var))
  
  # Extract vectors for table
  var_vector <- pull(unique_data, !!sym(var))
  cluster_vector <- pull(unique_data, "final_clusters")
  
  contingency_table <- table(var_vector, cluster_vector)
  
  # Print the table
  cat(sprintf("\nContingency table for %s (unique subjects):\n", var))
  print(contingency_table)
  # Save contingency table as CSV
  contingency_df <- as.data.frame.matrix(contingency_table)
  contingency_df$Category <- rownames(contingency_df)
  write.csv(contingency_df, file = sprintf("outputs/HCPC_%s_contingency_table.csv", var), row.names = FALSE)
  
  # Chi-squared test
  chi_test <- chisq.test(contingency_table)
  cat(sprintf("\nChi-squared test for %s:\n", var))
  print(chi_test)
  chi_df <- broom::tidy(chi_test)
  write.csv(chi_df, file = sprintf("outputs/HCPC_%s_chi_squared.csv", var), row.names = FALSE)
  
  # Check expected counts
  expected <- chi_test$expected
  if (any(expected < 5)) {
    cat(sprintf("\nWarning: Some expected counts < 5 for %s, performing Fisher's exact test.\n", var))
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    cat(sprintf("Fisher's exact test for %s:\n", var))
    print(fisher_test)
    fisher_df <- broom::tidy(fisher_test)
    write.csv(fisher_df, file = sprintf("outputs/HCPC_%s_fisher_exact.csv", var), row.names = FALSE)
  }
}

# Run tests for each categorical variable
lapply(factor_variables, test_categorical)

#########################################################################################

## Plot with healthy (asymptomatic) controls
library(ggplot2)
library(dplyr)
library(RColorBrewer)

load("data/df_RC_MWF_merged_2023.RData")

# Subjects where both legs are "Asymptomatic" (asymptomatic controls from df_MWF_2023)
asymp_controls <- df_2023 %>%
  filter(item %in% 1:60) %>%  #stance phase filter
  group_by(subject) %>%
  summarise(all_asymp = all(severity == "Asymptomatic")) %>%  # Check if all rows (both legs) are Asymptomatic
  filter(all_asymp) %>%
  pull(subject)  # Get subject IDs

# Extract data
healthy_df <- df_2023 %>%
  filter(subject %in% asymp_controls, 
         item %in% 1:60) %>%  # Stance phase
  mutate(cluster = "Healthy",   # Assign "Healthy" as the cluster label
         group = "Healthy")     # Unified group column

# Prepare clustered data with group column
clust_df <- dat_clust_plot %>%
  filter(item %in% 1:60) %>%  # Ensure same filter
  mutate(group = as.character(final_clusters))  # Convert cluster to character for consistency

#Combine into new data frame
dat_clust_healthy <- bind_rows(clust_df, healthy_df)

# check output 
print(paste("Unique subjects in clustered data:", n_distinct(clust_df$subject)))
print(paste("Unique subjects in healthy controls:", n_distinct(healthy_df$subject)))
print(table(dat_clust_healthy$group))

##Plot Powerpoint version 

plot_WFbyClust_with_controls_PP <- dat_clust_healthy %>%
  filter(signal_names == "KNEE_ANGLE") %>%
  filter(signal_components != "Z") %>%
  ggplot(aes(x = item, y = value, color = group, fill = group)) +  # Use group instead of cluster
  geom_hline(yintercept = 0, color = "black", linewidth = 0.25) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  # Uncomment if you want ribbons
  # stat_summary(
  #   fun.data = "mean_sdl",
  #   fun.args = list(mult = 1),
  #   geom = "ribbon",
  #   alpha = 0.25
  # ) +
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
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Mean Waveforms by Cluster with\nHealthy Comparison",
       x = "Gait Cycle (%)",
       y = "Mean Waveform Value",
       color = "Cluster",
       fill = "Cluster") +
  scale_color_viridis(discrete = TRUE, option = "D") +  # Viridis for lines
  scale_fill_viridis(discrete = TRUE, option = "D")

# Display plot
print(plot_WFbyClust_with_controls_PP)

##Plot report version 
plot_WFbyClust_with_controls <- dat_clust_healthy %>%
  filter(signal_names == "KNEE_ANGLE") %>%
  filter(signal_components != "Z") %>%
  filter(item %in% 1:60) %>%
  ggplot(aes(x = item, y = value, color = group, fill = group)) +  # Use group instead of cluster
  geom_hline(yintercept = 0, color = "black", linewidth = 0.25) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  # Uncomment if you want ribbons
  # stat_summary(
  #   fun.data = "mean_sdl",
  #   fun.args = list(mult = 1),
  #   geom = "ribbon",
  #   alpha = 0.25
  # ) +
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
    text = element_text(size = 12),
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Mean Waveforms by Cluster with Healthy Comparison",
       x = "Gait Cycle (%)",
       y = "Mean Waveform Value",
       color = "Cluster",
       fill = "Cluster") +
  scale_color_viridis(discrete = TRUE, option = "D") +  # Viridis for lines
  scale_fill_viridis(discrete = TRUE, option = "D")

# Display plot
print(plot_WFbyClust_with_controls)
ggsave("outputs/HC_WFbyClust_with_controls_plot.png", plot = plot_WFbyClust_with_controls, width = 8, height = 6, dpi = 300)

