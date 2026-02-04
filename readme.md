# Bilateral Gait Cluster Analysis

**Author:** Naomi Loughlin
**Created:** 2025-05-22
**Affiliation:** University of Waterloo

## Project Description

This project applies unsupervised clustering to bilateral knee gait waveforms from individuals with knee osteoarthritis (OA). The goal is to identify distinct gait pattern subgroups within a medial knee OA cohort using dimensionality reduction and cluster analysis, and then validate the resulting clusters using Statistical Parametric Mapping (SPM).

The analysis focuses on knee angle waveforms in the **frontal** (stance phase, 0-60% gait cycle) and **sagittal** (full gait cycle) planes. Individuals with lateral OA and isolated patellofemoral OA are excluded to maintain a homogeneous medial OA sample.

## Clustering Methods

Multiple clustering approaches were evaluated to determine the most appropriate method for this mixed-type dataset (continuous gait features + categorical demographics):

- **PAM** (Partitioning Around Medoids) with Gower distance
- **KAMILA** (KAy-means for MIxed LArge data)
- **PAM-MS** (Modha-Spangler extension of PAM for mixed data)
- **Hierarchical Clustering** (Ward's minimum variance, Manhattan distance)

All methods were assessed via subsampling stability (100 iterations, 80% fraction), silhouette width, Adjusted Rand Index (ARI), and Jaccard bootstrap stability. See `Bilateral_PAM_KAM_MS_with_MDS.Rmd` for the full comparison of PAM, KAMILA, and PAM-MS methods.

**Final method selected: Hierarchical Clustering (Ward's method)** with Manhattan distance on MDS-reduced coordinates (`scripts/Bilateral_HC_Ward.r`). This method was chosen based on superior stability and interpretability metrics.

## Analysis Pipeline

The pipeline is executed in order via `scripts/~main.R`:

```
1. Data Preprocessing
   ├── Bilateral_v3d_data2023_preprocessing.R
   ├── Bilateral_RC_MWF_merge_2023.R
   ├── Bilateral_RC_data2025_preprocessing.R
   ├── Bilateral_v3d_data2025_preprocessing.R
   ├── Bilateral_merge_RC_v3d_2025.R
   ├── Bilateral_one_leg_2023.R
   ├── Bilateral_one_leg_2025.R
   └── Bilateral_merge_2023_2025.R

2. Principal Component Analysis
   └── Bilateral_PCA.r
       (3 PCs per plane, helper: pca_zscores.R)

3. Multidimensional Scaling
   └── Bilateral_MDS.r
       (Manhattan distance, classical MDS to 2D)

4. Hierarchical Clustering (Final Method)
   └── Bilateral_HC_Ward.r
       (Ward's method, subsampling validation, Jaccard stability,
        cluster characterization, statistical tests, SPM data export)

5. Statistical Parametric Mapping (Python)
   ├── spm_analysis_v3.py        (between-cluster & OA vs healthy comparisons)
   └── paired_spm_analysis_v2.py (ipsilateral vs contralateral within clusters)
```

### Additional Analysis Scripts (Method Comparison)

- `Bilateral_PAM_KAM_MS_with_MDS.Rmd` -- Comprehensive comparison of PAM, KAMILA, and PAM-MS clustering on MDS coordinates. Includes optimal k selection, subsampling stability, Jaccard bootstrap, cross-method agreement, feature importance, and waveform visualization.
- `scripts/Bilateral_HCPC_clValid_Aug13.Rmd` -- Cluster validation using the `clValid` package across multiple distance metrics and clustering algorithms.

## Directory Structure

```
.
├── scripts/
│   ├── ~main.R                              # Runs the full R pipeline in order
│   ├── Bilateral_v3d_data2023_preprocessing.R
│   ├── Bilateral_RC_MWF_merge_2023.R
│   ├── Bilateral_RC_data2025_preprocessing.R
│   ├── Bilateral_v3d_data2025_preprocessing.R
│   ├── Bilateral_merge_RC_v3d_2025.R
│   ├── Bilateral_one_leg_2023.R
│   ├── Bilateral_one_leg_2025.R
│   ├── Bilateral_merge_2023_2025.R
│   ├── Bilateral_PCA.r
│   ├── Bilateral_MDS.r
│   ├── Bilateral_HC_Ward.r                  # Final clustering method
│   ├── Bilateral_PAM_KAM_MS.r              # PAM/KAMILA/PAM-MS script version
│   ├── Bilateral_HCPC_clValid_Aug13.Rmd     # Cluster validation
│   ├── spm_analysis_v3.py                   # SPM: cluster & healthy comparisons
│   ├── paired_spm_analysis_v2.py            # SPM: bilateral paired analysis
│   └── helpers/
│       ├── pca_zscores.R                    # PCA computation helper
│       ├── scr_plot.R                       # Single component reconstruction plots
│       └── shift_legend.R                   # Plot legend positioning helper
├── Bilateral_PAM_KAM_MS_with_MDS.Rmd       # Mixed clustering methods comparison
├── data/                                    # Input data (not tracked in git)
├── outputs/                                 # Generated results (not tracked in git)
├── Bilateral_Gait_Cluster_Analysis_2025.Rproj
└── readme.md
```

## Requirements

### R

- R >= 4.0
- Key packages: `cluster`, `factoextra`, `FactoMineR`, `fpc`, `clValid`, `clusterSim`, `mclust`, `kamila`, `dplyr`, `ggplot2`, `viridis`, `tableone`, `broom`, `here`, `fst`

### Python

- Python >= 3.8
- Packages: `numpy`, `pandas`, `spm1d`, `matplotlib`

## How to Run

### R Pipeline

1. Open `Bilateral_Gait_Cluster_Analysis_2025.Rproj` in RStudio.
2. Place the required raw data files in the `data/` directory.
3. Run `scripts/~main.R` to execute the full preprocessing, PCA, MDS, and clustering pipeline.
4. For the mixed clustering method comparison, knit `Bilateral_PAM_KAM_MS_with_MDS.Rmd`.

### SPM Analysis (Python)

After the R pipeline produces the SPM export CSVs in `data/spm_export/`:

```bash
cd scripts
python spm_analysis_v3.py
python paired_spm_analysis_v2.py
```

## Data

Raw gait data files are not included in this repository due to size and privacy constraints. The pipeline expects the following input files in `data/`:

- V3D gait exports (`.fst` format) for 2023 and 2025 cohorts
- RC (Range of Change) clinical data
- Severity and KL grade classifications

All intermediate `.RData` files (PCA results, MDS coordinates, cluster assignments) are regenerated by the pipeline.

## Outputs

All outputs are written to the `outputs/` directory (gitignored). These include:

- PCA reconstruction plots and loading vectors (PNG)
- MDS scatter plots colored by severity/KL grade (PNG)
- Cluster validation metrics (silhouette, ARI, Jaccard) (CSV)
- Cluster summary tables and statistical test results (CSV)
- Waveform comparison plots by cluster (PNG)
- SPM inference plots with significant regions highlighted (PNG)
- Effect size plots across the gait cycle (PNG)
