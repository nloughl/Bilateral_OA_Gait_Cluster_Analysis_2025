# Mixed data types clustering

library(dplyr)
library(cluster)
library(factoextra)
library(fpc)
library(kamila)
library(mclust)

## Load data 
load("data/df_2023_2025_combined.RData") # raw data fro sample -- used for demographics/cat variables 
load("data/mds_input_data.RData") # PCA results and numeric data before MDS 
# load("data/df_mds_coords.RData") # scaled MDS result data with subject and signal side data 
# mds_input_data <- mds_coords # rename mds data if using mds results as clustering data instead of raw PCA results

# Create results directory
dir.create("outputs/mixed_results", showWarnings = FALSE)

# Create dataframe for clustering (including both numeric and categorical variables)
mixed_data <- df_combined %>%
  group_by(subject, signal_side) %>%
  summarise(
    severity_contra = first(severity_contra),  # Take first if consistent
    sex = first(sex),
    kl_contra = first(kl_contra),
    knee_oa_location = first(oa_location),
  )

dim(mixed_data)

mixed_data <- mixed_data %>%
  left_join(mds_input_data, by = c("subject", "signal_side")) 

dim(mixed_data)


# Define clustering variables
Data.Cluster <- mixed_data %>%
  dplyr::select(subject, age, PC1_front, PC1_sag, PC2_sag, PC3_sag, speed, sex, severity_contra, knee_oa_location) %>%
  na.omit()

#Explicitly ungroup by subject 
Data.Cluster <- Data.Cluster %>% ungroup()
groups(Data.Cluster)  # Should be NULL

# Convert categorical variables to factors
Data.Cluster[sapply(Data.Cluster, is.character)] <- lapply(Data.Cluster[sapply(Data.Cluster, is.character)], as.factor)
# Data.Cluster$kl_contra <- as.factor(Data.Cluster$kl_contra) # only if including kl_contra as a clustering variable

# Separate continuous and categorical variables
conVars <- Data.Cluster %>%
  dplyr::select(age, PC1_front, PC1_sag, PC2_sag, PC3_sag, speed) %>%
  scale() %>%  # Standardize using scale()
  as.data.frame()

catVars <- Data.Cluster %>%
  dplyr::select(sex, severity_contra, knee_oa_location) %>%
  mutate(across(everything(), as.factor))

# Verify
str(conVars)  # Should be numeric
str(catVars)  # Should be factors

# Compute Gower's distance
gower_dist <- daisy(cbind(conVars, catVars), metric = "gower")
gower_mat <- as.matrix(gower_dist)

## Determine optimal number of clusters 
# Silhouette method
fviz_nbclust(gower_mat, cluster::pam, method = "silhouette", k.max = 20)

# Prediction strength
PredictionThreshold <- 0.5
pred_strength <- prediction.strength(gower_mat, Gmin = 2, Gmax = 10, M = 50, 
                                    clustermethod = claraCBI, cutoff = PredictionThreshold, 
                                    classification = "centroid", distances = TRUE)
print(pred_strength)

# Sub-sampling method 

#MS Fucntions
# Define L1 (Manhattan) and matching distances
L1Dist <- function(v1, v2) sum(abs(v1 - v2))
matchingDist <- function(v1, v2) sum(as.integer(v1) != as.integer(v2))

# PAM wrapper for Modha-Spangler
pammix <- function(conData, catData, conWeight, nclust, ...) {
  conData <- as.data.frame(conData)
  catData <- as.data.frame(catData)
  distMat <- daisy(x = cbind(conData, catData), metric = "gower",
                   weights = rep(c(conWeight, 1 - conWeight),
                                 times = c(ncol(conData), ncol(catData))))
  clustRes <- pam(x = distMat, k = nclust, diss = TRUE, ...)
  
  return(list(cluster = clustRes$clustering, nc = nclust,
              conCenters = conData[clustRes$id.med, , drop = FALSE],
              catCenters = catData[clustRes$id.med, , drop = FALSE]))
}


# Parameters
n_iterations <- 100  # Number of subsamples
fraction <- 0.8      # 80% of data
k_range <- 2:10      # Test k from 2 to 10
k_range_ms <- 2:10    # different range due to nclustwarnings from: nclust >= the number of categorical variable level combinations
n_samples <- round(nrow(Data.Cluster) * fraction)  # Number of rows to sample
set.seed(123)

# Store results
ari_results_pam <- matrix(NA, nrow = n_iterations, ncol = length(k_range))
ari_results_kam <- matrix(NA, nrow = n_iterations, ncol = length(k_range))
ari_results_ms <- matrix(NA, nrow = n_iterations, ncol = length(k_range))
colnames(ari_results_pam) <- colnames(ari_results_kam) <- colnames(ari_results_ms) <- paste0("k_", k_range)
sil_results_pam <- matrix(NA, nrow = n_iterations, ncol = length(k_range))
sil_results_ms <- matrix(NA, nrow = n_iterations, ncol = length(k_range))

# Cluster full dataset for PAM, KAMILA, and Modha-Spangler
full_clusters_pam <- list()
full_clusters_kam <- list()
full_clusters_ms <- list()

# Scale continuous features for PAM
categorical_features <- c("knee_oa_location", "severity_contra", "sex")
continuous_features <- c("PC1_front", "PC1_sag", "PC2_sag", "PC3_sag", "age", "speed")
#cluster_data <- scale(Data.Cluster[, continuous_features, drop = FALSE])

for (j in seq_along(k_range)) {
  k <- k_range[j]
  
  # PAM
  pam_result <- pam(cluster_data, k = k)
  full_clusters_pam[[j]] <- pam_result$clustering
  
  # KAMILA
  kam_result <- kamila(
    conVar = Data.Cluster[, continuous_features],
    catFactor = Data.Cluster[, categorical_features],
    numClust = k,
    numInit = 10
  )
  full_clusters_kam[[j]] <- kam_result$finalMemb
}

for (j in seq_along(k_range_ms)) {
  k <- k_range_ms[j]
  
  # Modha-Spangler
  ms_result <- gmsClust(
    conData = Data.Cluster[, continuous_features],
    catData = Data.Cluster[, categorical_features],
    nclust = k, clustFun = pammix, conDist = L1Dist, catDist = matchingDist
  )
  full_clusters_ms[[j]] <- ms_result$results$cluster
}

# Subsampling and ARI calculation
for (i in 1:n_iterations) {
  # Subsample 80% of the data
  sample_idx <- sample(1:nrow(Data.Cluster), size = n_samples, replace = FALSE)
  sampled_data <- Data.Cluster[sample_idx, ]
  sampled_cluster_data <- cluster_data[sample_idx, ]
  
  # Cluster subsample for each k and compute ARI and silhouette widths 
  for (j in seq_along(k_range)) {
    k <- k_range[j]
   
    # PAM
    pam_sub_result <- pam(sampled_cluster_data, k = k)
    pam_sub_clusters <- pam_sub_result$clustering
    full_clusters_sub_pam <- full_clusters_pam[[j]][sample_idx]
    ari_results_pam[i, j] <- adjustedRandIndex(full_clusters_sub_pam, pam_sub_clusters)
    
    sil_pam <- silhouette(pam_sub_clusters, dist(sampled_cluster_data, method = "euclidean"))
    sil_results_pam[i,j] <- mean(sil_pam[, 3], na.rm = TRUE)
  
    # # KAMILA -- dont necessarily need to calculate, as the model chooses optimal k when run 
    # kamila_sub_result <- kamila(
    #   conVar = sampled_data[, continuous_features],
    #   catFactor = sampled_data[, categorical_features],
    #   numClust = k,
    #   numInit = 10
    # )
    # kamila_sub_clusters <- kamila_sub_result$finalMemb
    # full_clusters_sub_kamila <- full_clusters_kam[[j]][sample_idx]
    # ari_results_kam[i, j] <- adjustedRandIndex(full_clusters_sub_kamila, kamila_sub_clusters)
  }
}
for (i in 1:n_iterations) {
  # Subsample 80% of the data
  sample_idx <- sample(1:nrow(Data.Cluster), size = n_samples, replace = FALSE)
  sampled_data <- Data.Cluster[sample_idx, ]
  sampled_cluster_data <- cluster_data[sample_idx, ]
  for (j in seq_along(k_range_ms)) {
    k <- k_range_ms[j]
    
    # Modha-Spangler
    ms_sub_result <- gmsClust(
      conData = sampled_data[, continuous_features],
      catData = sampled_data[, categorical_features],
      nclust = k,
      clustFun = pammix, conDist = L1Dist, catDist = matchingDist
    )
    ms_sub_clusters <- ms_sub_result$results$cluster
    full_clusters_sub_ms <- full_clusters_ms[[j]][sample_idx]
    ari_results_ms[i, j] <- adjustedRandIndex(full_clusters_sub_ms, ms_sub_clusters)
    
    dist_ms <- daisy(sampled_data[, c(continuous_features, categorical_features)], metric = "gower")
    sil_ms <- silhouette(ms_sub_clusters, dist_ms)
    sil_results_ms[i,j] <- mean(sil_ms[, 3], na.rm = TRUE)
  }
}


# Summarize ARI resutls 
ari_sil_summary <- data.frame(
  k = k_range,
  mean_ari_pam = colMeans(ari_results_pam, na.rm = TRUE),
  mean_sil_pam = colMeans(sil_results_pam, na.rm = TRUE),
  mean_ari_ms = colMeans(ari_results_ms, na.rm = TRUE),
  mean_sil_ms = colMeans(sil_results_ms, na.rm = TRUE)
)


print(ari_sil_summary)
write.csv(ari_sil_summary, file="outputs/mixed_results/PAM_MS_ari_sil_summary.csv", row.names = FALSE)

## Determine optimal number of clusters (k) for each clustering method
combined_score_pam <- ari_sil_summary$mean_ari_pam + ari_sil_summary$mean_sil_pam  # Simple sum (could adjust weighting)
optimal_k_idx_pam <- which.max(combined_score_pam)
numclust_pam <- k_range[optimal_k_idx_pam]

combined_score_ms <- ari_sil_summary$mean_ari_ms + ari_sil_summary$mean_sil_ms  # Simple sum (could adjust weighting)
optimal_k_idx_ms <- which.max(combined_score_ms)
numclust_ms <- k_range[optimal_k_idx_ms]

cat("\nSelected k with highest combined mean ARI and silhouette for PAM:", numclust_pam, "\n")
cat("Mean ARI for selected k:", round(ari_sil_summary$mean_ari_pam[optimal_k_idx_pam], 3), "\n")
cat("Mean Silhouette for selected k:", round(ari_sil_summary$mean_sil_pam[optimal_k_idx_pam], 3), "\n")

cat("\nSelected k with highest combined mean ARI and silhouette for MS:", numclust_ms, "\n")
cat("Mean ARI for selected k:", round(ari_sil_summary$mean_ari_ms[optimal_k_idx_ms], 3), "\n")
cat("Mean Silhouette for selected k:", round(ari_sil_summary$mean_sil_ms[optimal_k_idx_ms], 3), "\n")

numclust_pam = ari_summary$k[which.max(ari_summary$mean_ari_pam)]
numclust_kamila = ari_summary$k[which.max(ari_summary$mean_ari_kamila)]
numclust_ms = ari_summary$k[which.max(ari_summary$mean_ari_ms)]

# Run clustering algorithms with optimal number of clusters 
## PAM
pam_fit <- pam(gower_dist, k = numclust_pam, diss = TRUE)

## Modha-Spangler PAM
ms_fit <- gmsClust(conData = conVars, catData = catVars, nclust = numclust_ms,
                   clustFun = pammix, conDist = L1Dist, catDist = matchingDist)

## KAMILA
conVar <- as.data.frame(conVars)
catFactor <- as.data.frame(catVars)
#Data.Cluster <- Data.Cluster %>% ungroup()
kam_fit <- kamila(conVar, catFactor, numClust = 2:10, numInit = 50, maxIter = 50,
                  calcNumClust = "ps", predStrThresh = PredictionThreshold)
    #KAM will choose optimal cluster number from 2:10 based on prediction strength and will choose the smallest k with prediction strength over the threshold

# Add cluster labels to data
Data.Cluster$PAMcluster <- pam_fit$clustering
Data.Cluster$MScluster <- ms_fit$results$cluster
Data.Cluster$KAMcluster <- kam_fit$finalMemb

################################################################################

# Evaluate Cluster Stability 
library(fpc)

# Compute Adjusted Rand Index -- for cluster aggreement 
RAND_PAM_MS <- adjustedRandIndex(Data.Cluster$PAMcluster, Data.Cluster$MScluster)
RAND_PAM_KAM <- adjustedRandIndex(Data.Cluster$PAMcluster, Data.Cluster$KAMcluster)
RAND_KAM_MS <- adjustedRandIndex(Data.Cluster$MScluster, Data.Cluster$KAMcluster)
cluster_agreement <- data.frame(RAND_PAM_MS, RAND_PAM_KAM, RAND_KAM_MS)
write.csv(cluster_agreement, "outputs/mixed_results/mixed_data_cluster_agreement.csv", row.names = FALSE)

# Compute Jaccard stability 
#PAM 
# Functions to calculate stability, e.g., Jaccard Index for each cluster
jaccard_index_cluster <- function(cluster1, cluster2, cluster_id) {
  # Identify points that belong to the specific cluster
  cluster1_points <- which(cluster1 == cluster_id)
  cluster2_points <- which(cluster2 == cluster_id)
  
  # Calculate Jaccard index for the specific cluster
  intersect_len <- length(intersect(cluster1_points, cluster2_points))
  union_len <- length(union(cluster1_points, cluster2_points))
  
  # Return the Jaccard index for this cluster
  return(intersect_len / union_len)
}

## calculate jaccard index between full clusters 
jaccard_index_overall <- function(clustering1, clustering2) {
  ids1 <- unique(clustering1)
  ids2 <- unique(clustering2)
  matched <- 0
  total <- length(clustering1)
  for (i in 1:total) {
    if (clustering1[i] == clustering2[i]) matched <- matched + 1
  }
  return(matched / total)
}

# function for bootstrap resampling and clustering
bootstrap_pam <- function(conData, catData, k, B = 100) {
  # Store cluster results
  cluster_results <- list()
  
  # Initialize empty matrix for Jaccard stability for each cluster
  cluster_stability <- matrix(0, nrow = B, ncol = k)  # B bootstrap samples, k clusters
  
  for (i in 1:B) {
    # Bootstrap resample (with replacement)
    sample_data_con <- conData[sample(1:nrow(conData), replace = TRUE), ]
    sample_data_cat <- catData[sample(1:nrow(catData), replace = TRUE), ]
    
    # Create the distance matrix for PAM
    distMat <- daisy(x = cbind(sample_data_con, sample_data_cat), metric = "gower")
    
    # Run PAM with the specified k
    pam_fit <- pam(distMat, k = k, diss = TRUE)
    
    # Store the clustering results
    cluster_results[[i]] <- pam_fit$clustering
    
    # Calculate the stability for each cluster in this bootstrap sample
    if (i > 1) {
      for (cluster_id in 1:k) {
        cluster_stability[i, cluster_id] <- jaccard_index_cluster(cluster_results[[1]], pam_fit$clustering, cluster_id)
      }
    }
  }
  
  # Return the stability for each cluster across bootstrap samples
  return(list(
    stability = cluster_stability,
    clusterings = cluster_results
  ))
}

# Run for PAM (with k = 2 (numclust_pam))
pam_bootstrap_result <- bootstrap_pam(conVars, catVars, k = numclust_pam, B = 100)
cluster_stability_pam <- pam_bootstrap_result$stability
cluster_results_pam <- pam_bootstrap_result$clusterings

mean_stability_per_cluster <- colMeans(cluster_stability_pam, na.rm = TRUE)

# Calculate summary stats
jaccard_summary_pam <- data.frame(
  Cluster = 1:ncol(cluster_stability_pam),
  Mean_Jaccard = mean_stability_per_cluster,
  SD_Jaccard = apply(cluster_stability_pam, 2, sd, na.rm = TRUE)
)

# Compare results across bootstrap samples
jaccard_vals_pam <- sapply(1:(length(cluster_results_pam) - 1), function(i) {
  jaccard_index_global(cluster_results_pam[[i]], cluster_results_pam[[i + 1]])
})
overall_mean_jaccard_pam <- mean(jaccard_vals_pam)
print(overall_mean_jaccard_pam)


jaccard_summary_pam_extended <- rbind(
  jaccard_summary_pam,
  data.frame(
    Cluster = "Overall_Mean",
    Mean_Jaccard = overall_mean_jaccard_pam,
    SD_Jaccard = NA  # You can leave SD blank or compute it if needed
  )
)
# Save to CSV
write.csv(jaccard_summary_pam_extended, "outputs/mixed_results/jaccard_PAM_with_overall.csv", row.names = FALSE)


#KAM 
library(clue)  # for cluster matching (Hungarian algorithm)

# parameters
B <- 100  # number of bootstrap iterations
k <- length(unique(Data.Cluster$KAMcluster))    # number of clusters
set.seed(123)

# Store original clustering 
kam_orig <- kamila(conVar, catFactor, numClust = k, numInit = 10, maxIter = 50)
orig_clusters <- kam_orig$finalMemb
n <- nrow(conVar)

# Initialize Jaccard matrix 
jaccard_matrix <- matrix(NA, nrow = B, ncol = k)
kam_clusterings <- vector("list", B)  # stores bootstrap clusterings

# Bootstrap Loop 
for (b in 1:B) {
  # Resample indices
  boot_indices <- sample(1:n, replace = TRUE)
  oob_indices <- setdiff(1:n, boot_indices)  # optional: out-of-bag points
  
  # Bootstrap sample
  boot_con <- conVar[boot_indices, ]
  boot_cat <- catFactor[boot_indices, ]
  
  # Run KAMILA on bootstrap sample
  kam_boot <- kamila(boot_con, boot_cat, numClust = k, numInit = 10, maxIter = 50)
  boot_clusters <- kam_boot$finalMemb
  kam_clusterings[[b]] <- boot_clusters  # store current clustering
  
  # Get cluster labels for only the bootstrapped rows in original data
  true_labels <- orig_clusters[boot_indices]
  
  # Compute contingency table
  contingency <- table(true_labels, boot_clusters)
  
  # Match clusters using Hungarian algorithm
  matched <- clue::solve_LSAP(contingency, maximum = TRUE)
  
  # Compute Jaccard index for each cluster
  for (i in 1:k) {
    orig_members <- which(true_labels == i)
    boot_members <- which(boot_clusters == matched[i])
    intersect_size <- length(intersect(orig_members, boot_members))
    union_size <- length(union(orig_members, boot_members))
    jaccard_matrix[b, i] <- intersect_size / union_size
  }
}

# Summary statistics
colnames(jaccard_matrix) <- paste0("Cluster_", 1:k)
jaccard_summary <- data.frame(
  Cluster = 1:k,
  Mean_Jaccard = colMeans(jaccard_matrix, na.rm = TRUE),
  SD_Jaccard = apply(jaccard_matrix, 2, sd, na.rm = TRUE)
)

print(jaccard_summary) # measures how consistently the same points are grouped together (cluster stability) -- irrespective of cluster label 

# Function to compute Jaccard index between two clusterings
jaccard_index <- function(cl1, cl2) {
  n <- length(cl1)
  agree_matrix1 <- outer(cl1, cl1, FUN = "==")
  agree_matrix2 <- outer(cl2, cl2, FUN = "==")
  intersection <- sum(agree_matrix1 & agree_matrix2) - n  # remove diagonal
  union <- sum(agree_matrix1 | agree_matrix2) - n
  return(intersection / union)
}

# Calculate pairwise Jaccard indices between consecutive clusterings
jaccard_vals_kam <- sapply(1:(B - 1), function(i) {
  jaccard_index(kam_clusterings[[i]], kam_clusterings[[i + 1]])
})

# Compute overall mean
overall_mean_jaccard_kam <- mean(jaccard_vals_kam, na.rm = TRUE)

jaccard_summary_extended <- rbind(
  jaccard_summary,
  data.frame(
    Cluster = "Overall_Mean",
    Mean_Jaccard = overall_mean_jaccard_kam,
    SD_Jaccard = NA
  )
)

# Save
write.csv(jaccard_summary_extended, "outputs/mixed_results/jaccard_KAM_with_overall.csv", row.names = FALSE)


#MS CLuster Stability 
bootstrap_ms_pam <- function(conData, catData, k=7, B = 100, conWeight = 0.5) {
  # Store cluster results
  cluster_results <- list()
  
  # Initialize empty matrix for Jaccard stability for each cluster
  cluster_stability <- matrix(0, nrow = B, ncol = k)  # B bootstrap samples, k clusters
  
  for (i in 1:B) {
    # Bootstrap resample (with replacement)
    sample_data_con <- conData[sample(1:nrow(conData), replace = TRUE), ]
    sample_data_cat <- catData[sample(1:nrow(catData), replace = TRUE), ]
    
    # Run MS PAM
    ms_fit <- pammix(conData = sample_data_con, catData = sample_data_cat, nclust = k, conWeight = conWeight)
    
    # Store the clustering labels
    cluster_results[[i]] <- ms_fit$cluster
    
    # Calculate the stability for each cluster
    for (cluster_id in 1:k) {
      cluster_stability[i, cluster_id] <- jaccard_index_cluster(cluster_results[[1]], ms_fit$cluster, cluster_id)
    }
  }
  
  # Return the stability for each cluster across bootstrap samples
  return(list(
    stability = cluster_stability,
    results = cluster_results
  ))
}

# Run MS-PAM bootstrap 
ms_boot_output <- bootstrap_ms_pam(conVars, catVars, k = numclust_ms, B = 100, conWeight = 0.5)

# Extract components 
cluster_stability_ms <- ms_boot_output$stability
cluster_results_ms <- ms_boot_output$results

# calculate the mean Jaccard stability for each cluster
mean_stability_per_cluster <- colMeans(cluster_stability_ms)

mean_stability_per_cluster

jaccard_summary_ms <- data.frame(
  Cluster = 1:ncol(cluster_stability_ms),
  Mean_Jaccard = mean_stability_per_cluster,
  SD_Jaccard = apply(cluster_stability_ms, 2, sd, na.rm = TRUE)
)
print(jaccard_summary_ms)

# Run for MS =
cluster_results_ms <- bootstrap_ms_pam(conVars, catVars, k = numclust_ms, B = 100, conWeight = 0.5)

# Calculate stability -- Jaccard Index (uses same function as PAM)
jaccard_vals_ms <- sapply(1:(B - 1), function(i) {
  jaccard_index_overall(cluster_results_ms$results[[i]], cluster_results_ms$results[[i + 1]])
})

overall_mean_jaccard_ms <- mean(jaccard_vals_ms)  # Average Jaccard index for stability
print(overall_mean_jaccard_ms)

jaccard_summary_ms_extended <- rbind(
  jaccard_summary_ms,
  data.frame(
    Cluster = "Overall_Mean",
    Mean_Jaccard = overall_mean_jaccard_ms,
    SD_Jaccard = NA  # You can leave SD blank or compute it if needed
  )
)
# Save to CSV
write.csv(jaccard_summary_ms_extended, "outputs/mixed_results/jaccard_MS_with_overall.csv", row.names = FALSE)


#Compare Jaccard index between methods -- do the same points stay together, and are the clusters labelled the same across runs? 
df_compare <- data.frame(
  Method = c("PAM", "MS-PAM", "KAMILA"),
  Mean_Jaccard = c(mean(jaccard_vals_pam),
                   mean(jaccard_vals_ms),
                   mean(jaccard_vals_kam))
)

ggplot(df_compare, aes(x = Method, y = Mean_Jaccard)) +
  geom_col(fill = "steelblue") +
  ylim(0, 1) +
  labs(title = "Clustering Method Stability Comparison", y = "Mean Jaccard Index")

################################################################################
## Summarize clusters 
library(tableone)
library(pander)

variable_list <- c("sex", "severity_contra", "kl_contra", "knee_oa_location", "age", "PC1_front", "PC1_sag", "PC2_sag", "PC3_sag", "speed")
factor_variables <- c("sex", "severity_contra", "knee_oa_location", "kl_contra")

## Add kl_contra to dataset for patients with scores (not all have kl scores)
Data.Cluster <- Data.Cluster %>%
  left_join(mixed_data %>%
              dplyr::select("subject", "kl_contra"), by = "subject")

# Overall table
tab0 <- CreateTableOne(vars = variable_list, data = Data.Cluster, factorVars = factor_variables, test = FALSE)
print(tab0)

# PAM clusters
PAM_clusters <- CreateTableOne(vars = variable_list, strata = "PAMcluster", data = Data.Cluster, factorVars = factor_variables, test = FALSE)
print(PAM_clusters)

PAM_clusters_matrix <- print(PAM_clusters,
                     includeNA = TRUE,
                     showAllLevels = TRUE
)

write.csv(PAM_clusters_matrix, file = "outputs/mixed_results/PAM_clusters_summary.csv")

# Modha-Spangler clusters
MS_clusters <- CreateTableOne(vars = variable_list, strata = "MScluster", data = Data.Cluster, factorVars = factor_variables, test = FALSE)
print(MS_clusters)

MS_clusters_matrix <- print(PAM_clusters,
                             includeNA = TRUE,
                             showAllLevels = TRUE
)

write.csv(MS_clusters_matrix, file = "outputs/mixed_results/MS_clusters_summary.csv")

# KAMILA clusters
KAM_clusters <- CreateTableOne(vars = variable_list, strata = "KAMcluster", data = Data.Cluster, factorVars = factor_variables, test = FALSE)
print(KAM_clusters)

KAM_clusters_matrix <- print(KAM_clusters,
                             includeNA = TRUE,
                             showAllLevels = TRUE
)

write.csv(KAM_clusters_matrix, file = "outputs/mixed_results/KAM_clusters_summary.csv")


## Outlier detection ########################################
# Outlier detection function
IsOutlier <- function(x) {
  outlier <- (x < quantile(x, 0.25, na.rm = TRUE) - 3 * IQR(x, na.rm = TRUE)) |
             (x > quantile(x, 0.75, na.rm = TRUE) + 3 * IQR(x, na.rm = TRUE))
  as.numeric(outlier)
}

# Identify outliers for continuous variables
for (i in names(conVars)) {
  Data.Cluster[[paste0("Outlier_", i)]] <- IsOutlier(Data.Cluster[[i]])
}

# Create Exclude column
Data.Cluster$Exclude <- rowSums(Data.Cluster[, grepl("Outlier", names(Data.Cluster))], na.rm = TRUE)

# Check for outliers
has_outliers <- any(Data.Cluster$Exclude != 0)

if (has_outliers) {
  # Separate outliers
  outliers <- Data.Cluster[Data.Cluster$Exclude != 0, ]
  Data.Cluster.NoOutliers <- Data.Cluster[Data.Cluster$Exclude == 0, ]
  
  # Re-run clustering on non-outliers
  # Create conVars and catVars
  conVars.NoOutliers <- Data.Cluster.NoOutliers %>%
    dplyr::select(age, PC1_front, PC1_sag, PC2_sag,PC3_sag, speed) %>%
    mutate(across(everything(), as.numeric)) %>%
    scale() %>%
    as.data.frame()
  
  catVars.NoOutliers <- Data.Cluster.NoOutliers %>%
    dplyr::select(sex, severity_contra, knee_oa_location) %>%
    mutate(across(everything(), as.factor))
  
  # re-Run Modha-Spangler
  ms_fit <- gmsClust(conData = conVars.NoOutliers, catData = catVars.NoOutliers, nclust = numclust_ms,
                     clustFun = pammix, conDist = L1Dist, catDist = matchingDist)
  
  # re-Run KAMILA
  PredictionThreshold <- 0.5
  conVars.NoOutliers <- as.data.frame(conVars.NoOutliers)
  catFactor.NoOutliers <- as.data.frame(catVars.NoOutliers)
  #Data.Cluster.NoOutliers <- Data.Cluster.NoOutliers %>% ungroup()
  kam_fit <- kamila(conVars.NoOutliers, catFactor.NoOutliers, numClust = 2:10, numInit = 50, maxIter = 50,
                    calcNumClust = "ps", predStrThresh = PredictionThreshold)
  
  # re-run PAM for comparison
  gower_dist <- daisy(cbind(conVars.NoOutliers, catVars.NoOutliers), metric = "gower")
  pam_fit <- pam(gower_dist, k = numclust_pam, diss = TRUE)
  
  # Store results
  cat("Outliers detected:", nrow(outliers), "rows removed.\n")
  Data.Cluster <- Data.Cluster.NoOutliers  # Update Data.Cluster for downstream use
  conVars <- conVars.NoOutliers
  catVars <- catVars.NoOutliers
} else {
  cat("No outliers detected. Proceed with original Data.Cluster.\n")
}


#Statistical Analyses on clusters ####################################################

library(dplyr)
library(dunn.test)
library(broom)

# Select dataset
if (exists("Data.Cluster.NoOutliers")) {
  analysis_data <- Data.Cluster.NoOutliers
  cat("Using outlier-filtered data: Data.Cluster.NoOutliers (", nrow(analysis_data), "rows).\n")
} else {
  analysis_data <- Data.Cluster
  cat("Using original data: Data.Cluster (", nrow(analysis_data), "rows).\n")
}

# Ensure cluster columns are factors
analysis_data <- analysis_data %>%
  mutate(
    PAMcluster = as.factor(PAMcluster),
    KAMcluster = as.factor(KAMcluster),
    MScluster = as.factor(MScluster)
  )

# Feature importance function
feature_importance <- function(data, cluster_col, cluster_name, continuous_vars, categorical_vars, test_results) {
  cat(sprintf("\n=== Feature Importance for %s ===\n", cluster_name))
  
  # Extract p-values from test results
  p_values <- data.frame(
    variable = character(),
    p_value = numeric(),
    test = character(),
    stringsAsFactors = FALSE
  )
  
  # Track skipped variables
  skipped_vars <- data.frame(
    variable = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )
  
  # Continuous variables (Kruskal-Wallis or ANOVA)
  for (var in continuous_vars) {
    if (var %in% names(test_results) && !is.null(test_results[[var]]) && !is.null(test_results[[var]]$p_value)) {
      p_value <- test_results[[var]]$p_value
      test_type <- test_results[[var]]$test
      p_values <- rbind(p_values, data.frame(
        variable = var,
        p_value = p_value,
        test = test_type
      ))
    } else {
      skipped_vars <- rbind(skipped_vars, data.frame(
        variable = var,
        reason = "Test skipped (e.g., small clusters or all NAs)"
      ))
    }
  }
  
  # Categorical variables (Chi-squared or Fisher’s)
  for (var in categorical_vars) {
    if (var %in% names(test_results) && !is.null(test_results[[var]]) && !is.null(test_results[[var]]$p_value)) {
      p_value <- test_results[[var]]$p_value
      test_type <- test_results[[var]]$test
      p_values <- rbind(p_values, data.frame(
        variable = var,
        p_value = p_value,
        test = test_type
      ))
    } else {
      skipped_vars <- rbind(skipped_vars, data.frame(
        variable = var,
        reason = "Test skipped (e.g., small clusters or all NAs)"
      ))
    }
  }
  
  # Sort by p-value (lower = more important)
  p_values <- p_values %>%
    arrange(p_value)
  
  # Save p-value summary
  cat("\nFeature Importance (Statistical Tests, sorted by p-value):\n")
  print(p_values)
  write.csv(p_values, file = sprintf("outputs/mixed_results/feature_importance_%s.csv", cluster_name), row.names = FALSE)
  
  # skipped variables
  if (nrow(skipped_vars) > 0) {
    cat("\nSkipped Variables (no p-value available):\n")
    print(skipped_vars)
   }
  
  ##Centroid Analysis -- could also try permuatation importance 
  
  # Centroid analysis for continuous variables
  centroid_importance <- data.frame(
    variable = character(),
    centroid_sd = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (var in continuous_vars) {
    # Compute cluster means
    cluster_means <- data %>%
      group_by(.data[[cluster_col]]) %>%
      summarise(mean_value = mean(.data[[var]], na.rm = TRUE)) %>%
      pull(mean_value)
    # Standard deviation of means (spread across clusters)
    centroid_sd <- sd(cluster_means, na.rm = TRUE)
    centroid_importance <- rbind(centroid_importance, data.frame(
      variable = var,
      centroid_sd = centroid_sd
    ))
  }
  
  # Sort by centroid_sd (higher = more important)
  centroid_importance <- centroid_importance %>%
    arrange(desc(centroid_sd))
  
  # Save centroid importance
  cat("\nFeature Importance (Centroid Spread, sorted by SD):\n")
  print(centroid_importance)
  write.csv(centroid_importance, file = sprintf("outputs/mixed_results/centroid_importance_%s.csv", cluster_name), row.names = FALSE)
  
  # Plot centroid spread
  barplot(centroid_importance$centroid_sd, 
          names.arg = centroid_importance$variable,
          main = sprintf("Centroid Spread for %s", cluster_name),
          xlab = "Variable", ylab = "SD of Cluster Means",
          col = "lightblue")
  
  ## Categorical variable importance 
  categorical_importance <- data.frame(
    variable = character(),
    mean_entropy = numeric(),
    stringsAsFactors = FALSE
  )
  # entropy function (Shannon entropy)
  entropy_func <- function(p) {
    p <- p[p > 0]  # remove zeros to avoid log(0)
    -sum(p * log2(p))
  }
  
  categorical_importance <- data.frame(variable=character(), mean_entropy=numeric(), stringsAsFactors=FALSE)
  
  for (var in categorical_vars) {
    # Get proportion of each category within each cluster
    prop_table <- data %>%
      group_by(.data[[cluster_col]], .data[[var]]) %>%
      summarise(count = n(), .groups = "drop") %>%
      group_by(.data[[cluster_col]]) %>%
      mutate(prop = count / sum(count)) %>%
      ungroup()
    
    # Calculate entropy per cluster for this variable
    entropy_by_cluster <- prop_table %>%
      group_by(.data[[cluster_col]]) %>%
      summarise(entropy = entropy_func(prop), .groups = "drop")
    
    # Average entropy across clusters for this variable
    mean_ent <- mean(entropy_by_cluster$entropy, na.rm = TRUE)
    
    # Store result
    categorical_importance <- rbind(categorical_importance, data.frame(
      variable = var,
      mean_entropy = mean_ent
    ))
  }
  
  # Sort by ascending entropy (lower = more important)
  categorical_importance <- categorical_importance %>%
    arrange(mean_entropy)
  
  # Save categorical importance
  cat("\nCategorical Feature Importance:\n")
  print(categorical_importance)
  write.csv(categorical_importance, file = sprintf("outputs/mixed_results/categorical_importance_%s.csv", cluster_name), row.names = FALSE)
  
  # Plot categorical importance spread
  barplot(categorical_importance$mean_entropy, 
          names.arg = categorical_importance$variable,
          main = sprintf("Categorical Importance for %s", cluster_name),
          xlab = "Variable", ylab = "Mean Entropy",
          col = "lightblue")
  
}


# Continuous variable test function
test_continuous <- function(var, cluster_col, data, cluster_name) {
  cat(sprintf("\n=== Analysis for %s (%s) ===\n", var, cluster_name))
  # Check normality
  shapiro_test <- data %>%
    group_by(.data[[cluster_col]]) %>%
    summarise(p_value = shapiro.test(.data[[var]])$p.value)
  cat(sprintf("\nShapiro-Wilk test for %s:\n", var))
  print(shapiro_test)
  
  # ANOVA or Kruskal-Wallis
  test_result <- NULL
  p_value <- NA
  test_type <- ""
  if (any(shapiro_test$p_value < 0.05, na.rm = TRUE)) {
    test_result <- kruskal.test(reformulate(cluster_col, var), data = data)
    cat(sprintf("\nKruskal-Wallis test for %s:\n", var))
    print(test_result)
    p_value <- test_result$p.value
    test_type <- "Kruskal-Wallis"
    cat("Kruskal p-value:", p_value, "\n")
    # Save Kruskal-Wallis results
    kruskal_df <- data.frame(
      statistic = test_result$statistic,
      p_value = test_result$p.value,
      df = test_result$parameter,
      method = test_result$method
    )
    write.csv(kruskal_df, file = sprintf("outputs/mixed_results/%s_%s_kruskal.csv", var, cluster_name), row.names = FALSE)
  } else {
    test_result <- aov(reformulate(cluster_col, var), data = data)
    cat(sprintf("\nANOVA for %s:\n", var))
    anova_summary <- summary(test_result)
    print(anova_summary)
    p_values <- anova_summary[[1]][["Pr(>F)"]]
    p_value <- p_values[1]
    test_type <- "ANOVA"
    cat("ANOVA p-values:", p_values, "\n")
    cat("Selected ANOVA p-value:", p_value, "\n")
    # Save ANOVA results
    anova_df <- tidy(test_result)
    write.csv(anova_df, file = sprintf("outputs/mixed_results/%s_%s_anova.csv", var, cluster_name), row.names = FALSE)
  }
  
  # Pairwise tests if significant
  if (is.numeric(p_value) && !is.na(p_value) && p_value < 0.05) {
    if (any(shapiro_test$p_value < 0.05, na.rm = TRUE)) {
      cat(sprintf("\nDunn's test for %s (Bonferroni corrected):\n", var))
      dunn_result <- dunn.test(data[[var]], data[[cluster_col]], 
                               method = "bonferroni", kw = FALSE)
      print(dunn_result)
      # Save Dunn's test results
      dunn_df <- data.frame(
        comparison = gsub("-", "vs", dunn_result$comparisons),  # Replace 2-1 with 2vs1
        Z = dunn_result$Z,
        P_unadjusted = dunn_result$P,
        P_adjusted = dunn_result$P.adjusted
      )
      write.csv(dunn_df, file = sprintf("outputs/mixed_results/%s_%s_dunn.csv", var, cluster_name), row.names = FALSE)
    } else {
      cat(sprintf("\nTukey HSD test for %s:\n", var))
      tukey_result <- TukeyHSD(test_result)
      print(tukey_result)
      # Save Tukey's HSD results
      tukey_df <- tidy(tukey_result)
      tukey_df$contrast <- gsub("-", "vs", tukey_df$contrast)  # Replace 2-1 with 2vs1
      write.csv(tukey_df, file = sprintf("outputs/mixed_results/%s_%s_tukey.csv", var, cluster_name), row.names = FALSE)
    }
  } else {
    cat(sprintf("\nNo significant differences for %s (p = %s).\n", var, 
                if (is.na(p_value)) "NA" else round(p_value, 3)))
  }
  
  # Boxplot
  boxplot(reformulate(cluster_col, var), data = data, 
          main = sprintf("%s Across %s Clusters", var, cluster_name),
          xlab = "Cluster", ylab = var)
}

# Categorical variable test function
test_categorical <- function(var, cluster_col, data, cluster_name) {
  cat(sprintf("\n=== Analysis for %s (%s) ===\n", var, cluster_name))
  
  # Use unique subjects (should already be subject-level, but just in case)
  unique_data <- data %>%
    distinct(subject, .keep_all = TRUE) %>%
    dplyr::select(subject, !!sym(var), !!sym(cluster_col))
  
  # Extract vectors
  var_vector <- pull(unique_data, !!sym(var))
  cluster_vector <- pull(unique_data, !!sym(cluster_col))
  
  # Contingency table
  contingency_table <- table(var_vector, cluster_vector)
  cat(sprintf("\nContingency table for %s:\n", var))
  print(contingency_table)
  # Save contingency table
  write.csv(as.data.frame.matrix(contingency_table), 
            file = sprintf("outputs/mixed_results/%s_%s_contingency.csv", var, cluster_name), 
            row.names = TRUE)
  
  # Chi-squared test
  chi_test <- chisq.test(contingency_table)
  cat(sprintf("\nChi-squared test for %s:\n", var))
  print(chi_test)
  # Save chi-squared results
  chi_df <- tidy(chi_test)
  write.csv(chi_df, file = sprintf("outputs/mixed_results/%s_%s_chi.csv", var, cluster_name), row.names = FALSE)
  
  # Check expected counts
  p_value <- chi_test$p.value
  test_type <- "Chi-squared"
  if (any(chi_test$expected < 5)) {
    cat("Warning: Some expected counts < 5, running Fisher's exact test.\n")
    fisher_test <- fisher.test(contingency_table, simulate.p.value = TRUE)
    cat(sprintf("Fisher's exact test for %s:\n", var))
    print(fisher_test)
    p_value <- fisher_test$p.value
    test_type <- "Fisher’s Exact"
    fisher_df <- data.frame(
      p_value = fisher_test$p.value,
      method = "Fisher's Exact Test (simulated)"
    )
    write.csv(fisher_df, file = sprintf("outputs/mixed_results/%s_%s_fisher.csv", var, cluster_name), row.names = FALSE)
  }
  # Return test results for feature importance
  return(list(p_value = p_value, test = test_type))
}

# Define variables
continuous_vars <- c("age", "PC1_front", "PC1_sag", "PC2_sag", "PC3_sag", "speed")
categorical_vars <- c("sex", "severity_contra", "knee_oa_location")
cluster_methods <- c("PAMcluster", "KAMcluster", "MScluster")

# Loop over clustering methods
for (cluster_col in cluster_methods) {
  cluster_name <- sub("cluster", "", cluster_col)
  cat(sprintf("\n==================== Statistical Analysis for %s ====================\n", cluster_name))
  
  # Store test results for feature importance
  test_results <- list()
  
  # Continuous variables
  for (var in continuous_vars) {
    result <- test_continuous(var, cluster_col, analysis_data, cluster_name)
    if (!is.null(result)) {
      test_results[[var]] <- result
    }
  }
  
  # Categorical variables
  for (var in categorical_vars) {
    result <- test_categorical(var, cluster_col, analysis_data, cluster_name)
    if (!is.null(result)) {
      test_results[[var]] <- result
    }
  }
  
  # Run feature importance function
  feature_importance(analysis_data, cluster_col, cluster_name, continuous_vars, categorical_vars, test_results)
}


#####################################################################

#Plot Waveforms by Cluster 
library(viridis)

# Merge cluster assignment to waveform data 
dat_clust_plot <- df_combined %>%
  left_join(analysis_data[, c("subject", "PAMcluster", "KAMcluster", "MScluster")], join_by("subject"))

component_names <- c(
  `X`="Sagittal  Plane",
  `Y`="Frontal Plane",
  `Z`="Transverse Plane")

# PAM Plots
plot_WFbyClust_PAM_sag <- dat_clust_plot %>%
  filter(signal_names == "KNEE_ANGLE") %>%
  filter(signal_components == "X") %>%
  ggplot(aes(x = item, y = value, color = PAMcluster, fill = PAMcluster)) +  # Map both color and fill
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
    legend.position = 'right'
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Mean Waveforms Colored by PAM Cluster",
       x = "Gait Cycle (%)",
       y = "Mean Waveform Value",
       color = "Cluster",
       fill = "Cluster") +  
  scale_color_viridis(discrete = TRUE, option = "D") +  # Viridis for lines
  scale_fill_viridis(discrete = TRUE, option = "D")

plot_WFbyClust_PAM_front <- dat_clust_plot %>%
  filter(signal_names == "KNEE_ANGLE") %>%
  filter(signal_components == "Y") %>%
  filter(item %in% (1:60)) %>%
  ggplot(aes(x = item, y = value, color = PAMcluster, fill = PAMcluster)) +  # Map both color and fill
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
    legend.position = 'right'
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Mean Waveforms Colored by PAM Cluster",
       x = "Gait Cycle (%)",
       y = "Mean Waveform Value",
       color = "Cluster",
       fill = "Cluster") +  
  scale_color_viridis(discrete = TRUE, option = "D") +  # Viridis for lines
  scale_fill_viridis(discrete = TRUE, option = "D")

plot_WFbyClust_PAM_sag
plot_WFbyClust_PAM_front

# KAM Plot
plot_WFbyClust_KAM_sag <- dat_clust_plot %>%
  filter(signal_names == "KNEE_ANGLE") %>%
  filter(signal_components == "X") %>%
  ggplot(aes(x = item, y = value, color = KAMcluster, fill = KAMcluster)) +  # Map both color and fill
  geom_hline(yintercept = 0, color = "black", linewidth = 0.25) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  # stat_summary(
  #   fun.data = "mean_sdl",
  #   fun.args = list(mult = 1),
  #   geom = "ribbon",
  #   alpha = 0.25) +  
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
  labs(title = "Mean Waveforms Colored by PAM Cluster",
       x = "Gait Cycle (%)",
       y = "Mean Waveform Value",
       color = "Cluster",
       fill = "Cluster") +  
  scale_color_viridis(discrete = TRUE, option = "D") +  # Viridis for lines
  scale_fill_viridis(discrete = TRUE, option = "D")

plot_WFbyClust_KAM_sag

plot_WFbyClust_KAM_front <- dat_clust_plot %>%
  filter(signal_names == "KNEE_ANGLE") %>%
  filter(signal_components == "Y") %>%
  filter(item %in% (1:60)) %>%
  ggplot(aes(x = item, y = value, color = KAMcluster, fill = KAMcluster)) +  # Map both color and fill
  geom_hline(yintercept = 0, color = "black", linewidth = 0.25) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  # stat_summary(
  #   fun.data = "mean_sdl",
  #   fun.args = list(mult = 1),
  #   geom = "ribbon",
  #   alpha = 0.25) +  
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
  labs(title = "Mean Waveforms Colored by PAM Cluster",
       x = "Gait Cycle (%)",
       y = "Mean Waveform Value",
       color = "Cluster",
       fill = "Cluster") +  
  scale_color_viridis(discrete = TRUE, option = "D") +  # Viridis for lines
  scale_fill_viridis(discrete = TRUE, option = "D")

plot_WFbyClust_KAM_front

# MS Plot
plot_WFbyClust_MS_sag <- dat_clust_plot %>%
  filter(signal_names == "KNEE_ANGLE") %>%
  filter(signal_components == "X") %>%
  ggplot(aes(x = item, y = value, color = MScluster, fill = MScluster)) +  # Map both color and fill
  geom_hline(yintercept = 0, color = "black", linewidth = 0.25) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  # stat_summary(
  #   fun.data = "mean_sdl",
  #   fun.args = list(mult = 1),
  #   geom = "ribbon",
  #   alpha = 0.25) +  
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
  labs(title = "Mean Waveforms Colored by PAM Cluster",
       x = "Gait Cycle (%)",
       y = "Mean Waveform Value",
       color = "Cluster",
       fill = "Cluster") +  
  scale_color_viridis(discrete = TRUE, option = "D") +  # Viridis for lines
  scale_fill_viridis(discrete = TRUE, option = "D")

plot_WFbyClust_MS_sag

plot_WFbyClust_MS_front <- dat_clust_plot %>%
  filter(signal_names == "KNEE_ANGLE") %>%
  filter(signal_components == "Y") %>%
  filter(item %in% (1:60)) %>%
  ggplot(aes(x = item, y = value, color = MScluster, fill = MScluster)) +  # Map both color and fill
  geom_hline(yintercept = 0, color = "black", linewidth = 0.25) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  # stat_summary(
  #   fun.data = "mean_sdl",
  #   fun.args = list(mult = 1),
  #   geom = "ribbon",
  #   alpha = 0.25) +  
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
  labs(title = "Mean Waveforms Colored by PAM Cluster",
       x = "Gait Cycle (%)",
       y = "Mean Waveform Value",
       color = "Cluster",
       fill = "Cluster") +  
  scale_color_viridis(discrete = TRUE, option = "D") +  # Viridis for lines
  scale_fill_viridis(discrete = TRUE, option = "D")

plot_WFbyClust_MS_front


###############################################################################
# Compare Clusters across methods (Cross-method cluster stability) 
## Co-assignment matrix correlation 

pam_labels <- Data.Cluster$PAMcluster
kam_labels <- Data.Cluster$KAMcluster
ms_labels <- Data.Cluster$MScluster

# Function to compute co-clustering matrix
get_coclustering_matrix <- function(cluster_labels) {
  outer(cluster_labels, cluster_labels, FUN = "==") * 1  # convert TRUE/FALSE to 1/0
}

pam_mat    <- get_coclustering_matrix(pam_labels)
kam_mat <- get_coclustering_matrix(kam_labels)
ms_mat     <- get_coclustering_matrix(ms_labels)

# Flatten to a vector (avoid redundancy)
flatten <- function(mat) {
  mat[upper.tri(mat)]
}

# Jaccard similarity between two co-clustering vectors
jaccard_similarity <- function(a, b) {
  intersection <- sum(a & b)
  union <- sum(a | b)
  intersection / union
}

# Extract upper triangle vectors
pam_vec    <- flatten(pam_mat)
kam_vec <- flatten(kam_mat)
ms_vec     <- flatten(ms_mat)

# Pairwise comparisons
jaccard_pam_kam <- jaccard_similarity(pam_vec, kam_vec)
jaccard_pam_ms     <- jaccard_similarity(pam_vec, ms_vec)
jaccard_kam_ms  <- jaccard_similarity(kam_vec, ms_vec)

# View results
data.frame(
  Comparison = c("PAM vs KAMILA", "PAM vs MS", "KAMILA vs MS"),
  Jaccard = c(jaccard_pam_kam, jaccard_pam_ms, jaccard_kam_ms)
)


## Visualization
library(ggalluvial)
library(tidyr)
library(dplyr)
library(ggplot2)

# Ensure consistent factor levels for cluster labels
pam_f     <- factor(pam_labels)
kam_f  <- factor(kam_labels)
ms_f      <- factor(ms_labels)

# Combine into wide format
df_clusters <- data.frame(
  ID     = 1:length(pam_labels),
  PAM    = pam_f,
  KAMILA = kam_f,
  MS     = ms_f
)

# Reshape to long format for ggalluvial
df_long <- pivot_longer(df_clusters, 
                        cols = c(PAM, KAMILA, MS), 
                        names_to = "Method", 
                        values_to = "Cluster")

# Plot alluvial flow
ggplot(df_long,
       aes(x = Method, stratum = Cluster, alluvium = ID,
           fill = Cluster)) +
  geom_flow(stat = "alluvium", lode.guidance = "forward", alpha = 0.6) +
  geom_stratum(width = 0.3) +
  theme_minimal() +
  ggtitle("Cluster Membership Flow Across Methods") +
  theme(axis.text.x = element_text(size = 12),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")











