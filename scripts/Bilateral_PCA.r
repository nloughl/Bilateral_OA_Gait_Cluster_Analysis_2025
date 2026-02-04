# PCA For Clustering, knee angles (frontal and sagittal only)
# Code modified from /Volumes/hmrl-data/r-code/PCA_in_R_helpers/PCA_in_R_Example.R (Jereme Outerley's)

# July 20: Previous analyses including all (lateral and medial) OA locations were conducted. See saved html outputs dated before July 20, 2025. Will now redo analysis excluding lateral OA (as it has distinct gait affects from medial and may be dictating clustering over other features)
# Aug. 13: Will filter out isolated patellofemoral to see if it alters results (as rest of cases are medial). Results in Bilateral_HCPC_clValid_noisopatellAug13.html

set.seed(123)  # for reproducibility

library(here)
library(fst)
library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork) # needed for plots
library(factoextra) # needed for fviz_eig (nice scree plots)

# Load the data (combined 2023/2025 data set --one leg, walk, first visit only, at least one severe leg)
load("data/df_2023_2025_combined.RData")


## Prepare data for PCA 
varlist_long<-c("subject", "signal_side", "item", "value", "severity_contra", "severity", "sex", "signal_names", "signal_components", "knee_oa_location", "kl_score", "kl_contra")
varlist_wide<-c("subject", "signal_side", "severity","severity_contra", "sex", "signal_names", "signal_components", "knee_oa_location", "kl_score", "kl_contra")
itemlist_stance <- as.character(1:1:60)       # Stance phase: 60 points
itemlist_full_cycle <- as.character(1:1:101)  # Full gait cycle

col_head_order_stance<- c(varlist_wide,itemlist_stance)
col_head_order_full <- c(varlist_wide,itemlist_full_cycle)

#rename variables
df_for_train_mean <- df_combined %>%  
  rename("knee_oa_location" = "oa_location") %>%
  dplyr::select(all_of(varlist_long))#%>%
# mutate(knee=paste0(subject, signal_side))

######## Filter OA location ########
df_for_train_mean <- df_for_train_mean %>%
  filter(!knee_oa_location %in% c("lateral", "lateral, patellofemoral", "medial, lateral")) %>%
  filter(!knee_oa_location %in% c("isolated patellofemoral"))


length(unique(df_for_train_mean$subject))
####################################

# Split into frontal and sagittal datasets (frontal = stance phase, sag = entire gait cycle)
df_frontal <- df_for_train_mean %>%
  filter(signal_names == "KNEE_ANGLE", signal_components == "Y") %>%
  filter(item %in% itemlist_stance)

df_sagittal <- df_for_train_mean %>%
  filter(signal_names == "KNEE_ANGLE", signal_components == "X") %>%
  filter(item %in% itemlist_full_cycle)


# pivot to make columns x percent of gait cycle
df_frontal_wide <- df_frontal %>%
  pivot_wider(values_from = "value", names_from = "item") %>%
  ungroup()

df_sagittal_wide <- df_sagittal %>%
  pivot_wider(values_from = "value", names_from = "item") %>%
  ungroup()

# ensure correct order of columns
df_frontal_wide <- df_frontal_wide[, col_head_order_stance]
df_sagittal_wide <- df_sagittal_wide[, col_head_order_full]

# Ensure factors are properly formatted for PCA to work properly 
df_frontal_wide <- df_frontal_wide %>%
  mutate(
    severity = factor(severity),
    severity_contra = factor(severity_contra),
    signal_names = factor(signal_names),
    signal_components = factor(signal_components),
    sex = factor(sex),
    kl_score = factor(kl_score),
    kl_contra = factor(kl_score)
  )

df_sagittal_wide <- df_sagittal_wide %>%
  mutate(
    severity = factor(severity),
    severity_contra = factor(severity_contra),
    signal_names = factor(signal_names),
    signal_components = factor(signal_components),
    sex = factor(sex),
    kl_score = factor(kl_score),
    kl_contra = factor(kl_score)
  )

## Principal Component Analysis (Stance phase, frontal and sagittal plane knee angles only)

source(file=file.path(here(), 'scripts', 'helpers', 'pca_zscores.R'), echo=FALSE)
# NOTE pca_zscore is cranky re: factors. Looks for numeric only columns.
# Any other columns put into df_for_PCA should be factors or omitted if possible.
#rank = number of PCs to retain
#scale = TRUE or FALSE (note: data is centred in pca_zscores)
df.pca.frontal <- pca_zscores(df_frontal_wide, scale = TRUE, rank = 3)
df.pca.sagittal <- pca_zscores(df_sagittal_wide, scale = TRUE, rank = 3)


# Merge back to get  severity groups
meta_frontal <- df_frontal_wide %>% dplyr::select(subject, severity, severity_contra, kl_score, kl_contra) %>% distinct()
meta_sagittal <- df_sagittal_wide %>% dplyr::select(subject, severity, severity_contra, kl_score, kl_contra) %>% distinct()

df_analysis_PC_frontal <- merge(meta_frontal, as.data.frame(df.pca.frontal$zscores))
df_analysis_PC_sagittal <- merge(meta_sagittal, as.data.frame(df.pca.sagittal$zscores))


# Should Export fst here as if you rerun the PC directions may change...
# Exports the dataframe to fst file for speed of import later. 
# You need to make a chkpts folder in the root of wherever "here" is or in the "data" folder of the project. type here() in console to find out
write.fst(df_analysis_PC_frontal, here('data', 'chkpts', paste0(Sys.Date(), '_frontal_PC.fst')))
write.fst(df_analysis_PC_sagittal, here('data', 'chkpts', paste0(Sys.Date(), '_sagittal_PC.fst')))


summary(df.pca.frontal$pca_list$KNEE_ANGLE_Y)  #Frontal plane
summary(df.pca.sagittal$pca_list$KNEE_ANGLE_X) #Sagittal
#summary(df.pca$pca_list$KNEE_ANGLE_Z) # transverse 

# Save summaries to .txt files (reference for analysis -- typically only include PCs with 80% culmulative proportion or more) 
pca_summary_text_front <- capture.output(summary(df.pca.frontal$pca_list$KNEE_ANGLE_Y))
writeLines(pca_summary_text_front, "outputs/pca_summary_frontal.txt")
pca_summary_text_sag <- capture.output(summary(df.pca.sagittal$pca_list$KNEE_ANGLE_X))
writeLines(pca_summary_text_sag, "outputs/pca_summary_sagittal.txt")

#fviz_eig(df.pca$pca_list$KNEE_ANGLE_X, ncp = 3)

# Build Single Component Reconstruction
source(file=file.path(here(), 'scripts', 'helpers', 'scr_plot.R'), echo=FALSE)

#Plots mean +/- 1SD of PC scores (reconstructions - different from example high and low PC example waveforms)
PC1_sag_plot <- scr_plot(df.pca.sagittal$pca_list$KNEE_ANGLE_X, "PC1", "Stance Only Knee Angle - Flexion/Extension - PC1")
PC2_sag_plot <- scr_plot(df.pca.sagittal$pca_list$KNEE_ANGLE_X, "PC2", "Stance Only Knee Angle - Flexion/Extension - PC2")
PC3_sag_plot <- scr_plot(df.pca.sagittal$pca_list$KNEE_ANGLE_X, "PC3", "Stance Only Knee Angle - Flexion/Extension - PC3")
PC1_frontal_plot <- scr_plot(df.pca.frontal$pca_list$KNEE_ANGLE_Y, "PC1", "Stance only Knee Angle - Ab-/Aduction - PC1")
PC2_frontal_plot <-scr_plot(df.pca.frontal$pca_list$KNEE_ANGLE_Y, "PC2", "Stance only Knee Angle - Ab-/Aduction - PC2")
PC3_frontal_plot <-scr_plot(df.pca.frontal$pca_list$KNEE_ANGLE_Y, "PC3", "Stance only Knee Angle - Ab-/Aduction - PC3")

PC1_sag_plot
PC2_sag_plot
PC3_sag_plot
PC1_frontal_plot
PC2_frontal_plot
PC3_frontal_plot

## Save  plots with ggsave
ggsave("outputs/PC1_sag_plot.png", plot = PC1_sag_plot, width = 12, height = 8, dpi = 300)
ggsave("outputs/PC2_sag_plot.png", plot = PC2_sag_plot, width = 12, height = 8, dpi = 300)
ggsave("outputs/PC3_sag_plot.png", plot = PC3_sag_plot, width = 12, height = 8, dpi = 300)
ggsave("outputs/PC1_frontal_plot.png", plot = PC1_frontal_plot, width = 12, height = 8, dpi = 300)
ggsave("outputs/PC2_frontal_plot.png", plot = PC2_frontal_plot, width = 12, height = 8, dpi = 300)
ggsave("outputs/PC3_frontal_plot.png", plot = PC3_frontal_plot, width = 12, height = 8, dpi = 300)



####### Loading vector plots #########################################################################

# Define time_points based on the number of rows in the rotation matrix
n_time_points_sag <- nrow(df.pca.sagittal$pca_list$KNEE_ANGLE_X$rotation)  # Should be 101
n_time_points_front <- nrow(df.pca.frontal$pca_list$KNEE_ANGLE_Y$rotation)  # Check this value
time_points_sag <- seq(0, 100, length.out = n_time_points_sag)
time_points_front <- seq(0, 100, length.out = n_time_points_front)

# Function to create loading plots with consistent styling
plot_loading <- function(pca_obj, pc_name, plane, time_points) {
  if (nrow(pca_obj$rotation) != length(time_points)) {
    stop(paste("Mismatch: rotation matrix has", nrow(pca_obj$rotation), "rows, but time_points has", length(time_points), "elements."))
  }
  loadings <- as.data.frame(pca_obj$rotation[, pc_name])
  colnames(loadings) <- "Loading"
  loadings <- loadings %>%
    mutate(Time = time_points) %>%
    rowid_to_column(var = "Index")
  
  p1 <- scr_plot(pca_obj, pc_name, paste0("KNEE_ANGLE_", plane))  # Assuming scr_plot works with PC name
  p2 <- ggplot(loadings, aes(x = Time, y = Loading)) +
    geom_line(linewidth = 1, color = "#1f77b4") +  # Distinctive color
    geom_hline(yintercept = 0, color = "black", linewidth = 1) +
    geom_hline(yintercept = c(0.3, -0.3), linetype = "dashed", color = "gray", linewidth = 0.5) +  # Threshold for significance
    theme_minimal() +
    theme(
      axis.line = element_line(linewidth = 1, colour = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      legend.position = "none"
    ) +
    labs(title = paste(plane, "Plane", pc_name, "Loadings"),
         x = "Gait Cycle (%)", y = "Loading") +
    scale_x_continuous(expand = c(0, 0), breaks = seq(0, 100, by = 25)) +
    coord_cartesian(ylim = c(-1, 1))  # Fixed y-axis for consistency
  
  loading_plot <- p1 / p2
  ggsave(paste0("outputs/loading_vector_plot_", tolower(plane), "_", pc_name, ".png"),
         plot = loading_plot, width = 12, height = 8, dpi = 300)
  return(loading_plot)
}

# Generate plots for frontal plane
plot_loading(df.pca.frontal$pca_list$KNEE_ANGLE_Y, "PC1", "Y", time_points_front)
plot_loading(df.pca.frontal$pca_list$KNEE_ANGLE_Y, "PC2", "Y", time_points_front)

# Generate plots for sagittal plane
plot_loading(df.pca.sagittal$pca_list$KNEE_ANGLE_X, "PC1", "X", time_points_sag)
plot_loading(df.pca.sagittal$pca_list$KNEE_ANGLE_X, "PC2", "X", time_points_sag)
plot_loading(df.pca.sagittal$pca_list$KNEE_ANGLE_X, "PC3", "X", time_points_sag)

# Print and inspect loadings for key PCs
print("Frontal Plane PC1 Loadings:")
print(df.pca.frontal$pca_list$KNEE_ANGLE_Y$rotation[, "PC1"])
print("Sagittal Plane PC1 Loadings:")
print(df.pca.sagittal$pca_list$KNEE_ANGLE_X$rotation[, "PC1"])
print("Sagittal Plane PC2 Loadings:")
print(df.pca.sagittal$pca_list$KNEE_ANGLE_X$rotation[, "PC2"])
print("Sagittal Plane PC3 Loadings:")
print(df.pca.sagittal$pca_list$KNEE_ANGLE_X$rotation[, "PC3"])


##############################################################################

#All knees (Frontal plane)
KneeFrontalPlanePC_stance<- df_analysis_PC_frontal %>%
  filter(signal_names == 'KNEE_ANGLE') %>% 
  filter(signal_components == 'Y') #%>%  #Filter for frontal plane (Y component)
  #filter(severity != "Replaced") %>% #Exclude patients with replacements 
  #filter(knee_oa_location != "Lateral") %>% 
  #filter(knee_oa_location != "Lateral_Patellofemoral") #%>%

#All knees (Sag plane)
KneeSagPlanePC_full<- df_analysis_PC_sagittal %>%
  filter(signal_names == 'KNEE_ANGLE') %>% 
  filter(signal_components == 'X') #%>%  #Filter for sagittal plane (X component)
  #filter(severity != "Replaced") %>% #Exclude patients with replacements 
  #filter(knee_oa_location != "Lateral") %>% 
  #filter(knee_oa_location != "Lateral_Patellofemoral") #%>%

# Plot PC1 by severity_contra score 
boxplot(PC1 ~ severity_contra, data =KneeFrontalPlanePC_stance, main="Stance Only Knee Frontal Plane", xlab = "Contralateral Severity")
boxplot(PC1 ~ severity_contra, data =KneeSagPlanePC_full, main="Knee Sagittal Plane", xlab = "Contralateral Severity")

# Statistical Analysis comparing PC1 by severity contra 
library(broom)
anova_PC1_contra_frontal<-aov(PC1 ~ severity_contra, data =KneeFrontalPlanePC_stance)
summary(anova_PC1_contra_frontal)
TukeyHSD(anova_PC1_contra_frontal)
anova_tidy_PC1_front <- tidy(anova_PC1_contra_frontal)
write.csv(anova_tidy_PC1_front, "outputs/PCA_anova_sevcontra_front_PC1.csv", row.names = FALSE)

anova_PC1_contra_sag<-aov(PC1 ~ severity_contra, data =KneeSagPlanePC_full)
summary(anova_PC1_contra_sag)
TukeyHSD(anova_PC1_contra_sag)
anova_tidy_PC1_sag <- tidy(anova_PC1_contra_sag)
write.csv(anova_tidy_PC1_sag, "outputs/PCA_anova_sevcontra_sag_PC1.csv", row.names = FALSE)

# Statistical Analysis comparing PC1 by kl contra 
library(broom)
anova_PC1_kl_contra_frontal<-aov(PC1 ~ kl_contra, data =KneeFrontalPlanePC_stance)
summary(anova_PC1_kl_contra_frontal)
TukeyHSD(anova_PC1_kl_contra_frontal)
anova_tidy_PC1_kl_front <- tidy(anova_PC1_kl_contra_frontal)
write.csv(anova_tidy_PC1_kl_front, "outputs/PCA_anova_klcontra_front_PC1.csv", row.names = FALSE)

anova_PC1_kl_contra_sag<-aov(PC1 ~ kl_contra, data =KneeSagPlanePC_full)
summary(anova_PC1_kl_contra_sag)
TukeyHSD(anova_PC1_kl_contra_sag)
anova_tidy_PC1_kl_sag <- tidy(anova_PC1_kl_contra_sag)
write.csv(anova_tidy_PC1_kl_sag, "outputs/PCA_anova_klcontra_sag_PC1.csv", row.names = FALSE)

# Plot PC1 by sex
boxplot(PC1 ~ sex, data =KneeFrontalPlanePC_stance, main="Stance Only Knee Frontal Plane", xlab = "Sex")
boxplot(PC1 ~ sex, data =KneeSagPlanePC_full, main="Knee Sagittal Plane", xlab = "Sex")

#Summary statistics of PCs by severity_contra groups 
library(broom)
KneeFrontalPlanePC_stance %>%
  group_by(severity_contra) %>%
  summarise(mean = mean(PC1), SD = sd(PC1), n = n())
KneeSagPlanePC_full %>%
  group_by(severity_contra) %>%
  summarise(mean = mean(PC1), SD = sd(PC1), n = n())

# Plot PC2 by severity contra score
boxplot(PC2 ~ severity_contra, data =KneeFrontalPlanePC_stance, xlab = "Contralateral Severity") 
anova_frontal_PC2<-aov(PC2 ~ severity_contra, data =KneeFrontalPlanePC_stance)
summary(anova_frontal_PC2)
anova_tidy_PC2_front <- tidy(anova_frontal_PC2)
write.csv(anova_tidy_PC2_front, "outputs/PCA_anova_front_PC2.csv", row.names = FALSE)

boxplot(PC2 ~ severity_contra, data =KneeSagPlanePC_full, xlab = "Contralateral Severity") 
anova_sag_PC2<-aov(PC2 ~ severity_contra, data =KneeSagPlanePC_full)
summary(anova_sag_PC2)
anova_tidy_PC2_sag <- tidy(anova_sag_PC2)
write.csv(anova_tidy_PC2_sag, "outputs/PCA_anova_sag_PC2.csv", row.names = FALSE)

# Plot PC3 by KL_contra score
boxplot(PC3 ~ severity_contra, data =KneeFrontalPlanePC_stance, xlab = "Contralateral Severity")             
anova_frontal_PC3<-aov(PC3 ~ severity_contra, data =KneeFrontalPlanePC_stance)
summary(anova_frontal_PC3)


table(KneeFrontalPlanePC_stance$knee_oa_location)
table(KneeFrontalPlanePC_stance$severity)

table(KneeFrontalPlanePC_stance$knee_oa_location)
table(KneeFrontalPlanePC_stance$severity)

#Frontal Plane (ab/duction) -- ALL PCs 
KneeFrontalPlanePCAll<- df_analysis_PC_frontal  %>% filter(signal_names == 'KNEE_ANGLE') %>% filter(signal_components == 'Y') 
save(KneeFrontalPlanePCAll, file=("data/PCA_frontal_plane_results.RData"))
boxplot(PC1 ~ severity_contra, data =KneeFrontalPlanePCAll, xlab = "Contralateral Severity")
anova1<-aov(PC1 ~ factor(severity_contra), data =KneeFrontalPlanePCAll)
summary(anova1)
TukeyHSD(anova1)


#Sag Plane (flex/exten) -- ALL PCs 
KneeSagPlanePCAll<- df_analysis_PC_sagittal  %>% filter(signal_names == 'KNEE_ANGLE') %>% filter(signal_components == 'X') 
save(KneeSagPlanePCAll, file=("data/PCA_sagittal_plane_results.RData"))
boxplot(PC1 ~ severity_contra, data =KneeSagPlanePCAll, xlab = "Contralateral Severity")
anova2<-aov(PC1 ~ factor(severity_contra), data =KneeSagPlanePCAll)
summary(anova2)
TukeyHSD(anova2)

table(KneeFrontalPlanePCAll$knee_oa_location)
table(KneeFrontalPlanePCAll$severity_contra)
table(KneeFrontalPlanePCAll$kl_contra)
table(interaction(KneeFrontalPlanePCAll$knee_oa_location,KneeFrontalPlanePCAll$severity))


table(KneeSagPlanePCAll$knee_oa_location)
table(KneeSagPlanePCAll$severity_contra)
table(KneeSagPlanePCAll$kl_contra)
table(interaction(KneeSagPlanePCAll$knee_oa_location,KneeSagPlanePCAll$severity_contra))

## Further investigate PC1 results 

# print PC1 rotations
print(df.pca.frontal$pca_list$KNEE_ANGLE_Y$rotation[,1])  # Loadings for PC1 (frontal plane)
print(df.pca.sagittal$pca_list$KNEE_ANGLE_X$rotation[,1])  # Loadings for PC1 (sag plane)

#plot PC1 results against waveform time points 
  ## Positive peaks show where high PC1 subjects differ most from low PC1
PC_loadings_frontal <- plot(df.pca.frontal$pca_list$KNEE_ANGLE_Y$rotation[,1], type = "l", xlab = "Waveform Time Point", ylab = "PC1 Loading",
     main = "Frontal Plane PC1 Loadings Across Waveform")
PC_loadings_sag <- plot(df.pca.sagittal$pca_list$KNEE_ANGLE_X$rotation[,1], type = "l", xlab = "Waveform Time Point", ylab = "PC1 Loading",
     main = "Sagittal Plane PC1 Loadings Across Waveform")
ggsave("outputs/loadings_plot_front_PC1.png", plot = PC_loadings_frontal, width = 12, height = 8, dpi = 300)
ggsave("outputs/loadings_plot_sag_PC1.png", plot = PC_loadings_sag, width = 12, height = 8, dpi = 300)

  
# Colour by severity_contra
library(ggplot2)
ggplot(KneeFrontalPlanePC_stance, aes(x = PC1, y = PC2, color = factor(severity_contra))) +
  geom_point() +
  labs(title = "Frontal Plane PCA Score Plot by Contralateral Severity", color = "severity_contra")

ggplot(KneeSagPlanePC_full, aes(x = PC1, y = PC2, color = factor(severity_contra))) +
  geom_point() +
  labs(title = "Sagittal Plane PCA Score Plot by Contralateral Severity", color = "severity_contra")

