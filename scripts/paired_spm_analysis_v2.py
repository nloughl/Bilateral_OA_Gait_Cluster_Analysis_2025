#!/usr/bin/env python3
"""
Paired SPM Analysis: Ipsilateral vs Contralateral within Clusters
WITH EFFECT SIZE CALCULATIONS
"""

import numpy as np
import pandas as pd
import spm1d
import matplotlib.pyplot as plt
from pathlib import Path

DATA_DIR = Path('../data/spm_export/bilateral')
OUTPUT_DIR = Path('../outputs/bilateral_analysis')
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

ALPHA = 0.05
DPI = 300


def load_paired_data(base_name):
    """Load ipsi and contra data, ensure pairing."""
    
    ipsi_file = DATA_DIR / f'{base_name}_ipsi.csv'
    contra_file = DATA_DIR / f'{base_name}_contra.csv'
    
    ipsi_data = pd.read_csv(ipsi_file)
    contra_data = pd.read_csv(contra_file)
    
    # CRITICAL: Verify same subjects in same order
    if not ipsi_data['subject'].equals(contra_data['subject']):
        print("WARNING: Subjects not matched - attempting to align...")
        # Merge and align
        merged = ipsi_data.merge(contra_data, on=['subject', 'cluster'], 
                                suffixes=('_ipsi', '_contra'))
        print(f"  Aligned: {len(merged)} paired subjects")
        
        # Extract waveforms
        ipsi_cols = [c for c in merged.columns if c.endswith('_ipsi')]
        contra_cols = [c for c in merged.columns if c.endswith('_contra')]
        
        # Remove suffix from column names and sort
        ipsi_waves = merged[ipsi_cols].values
        contra_waves = merged[contra_cols].values
        
        cluster_labels = merged['cluster'].values
        subjects = merged['subject'].values
        
    else:
        # Already aligned
        ipsi_waves = ipsi_data.iloc[:, 2:].values  # Skip subject, cluster
        contra_waves = contra_data.iloc[:, 2:].values
        cluster_labels = ipsi_data['cluster'].values
        subjects = ipsi_data['subject'].values
    
    return ipsi_waves, contra_waves, cluster_labels, subjects


def calculate_cohens_d_paired(ipsi, contra):
    """
    Calculate Cohen's d for paired samples at each time point.
    
    For paired samples:
    d = mean(differences) / std(differences)
    
    Parameters:
    -----------
    ipsi : ndarray
        Ipsilateral waveforms (n_subjects × n_timepoints)
    contra : ndarray
        Contralateral waveforms (n_subjects × n_timepoints)
    
    Returns:
    --------
    cohens_d : ndarray
        Effect size at each time point (n_timepoints,)
    """
    # Calculate differences for each subject
    differences = ipsi - contra
    
    # Cohen's d for paired samples
    mean_diff = differences.mean(axis=0)
    std_diff = differences.std(axis=0, ddof=1)  # Sample std deviation
    
    # Avoid division by zero
    cohens_d = np.divide(mean_diff, std_diff, 
                         out=np.zeros_like(mean_diff), 
                         where=std_diff!=0)
    
    return cohens_d


def paired_spm_by_cluster(ipsi_all, contra_all, clusters, cluster_val, cluster_name):
    """Perform paired SPM for one cluster WITH EFFECT SIZE."""
    
    # Select subjects from this cluster
    mask = clusters == cluster_val
    ipsi = ipsi_all[mask]
    contra = contra_all[mask]
    
    print(f"\n{'='*70}")
    print(f"Paired SPM: {cluster_name} - Ipsilateral vs Contralateral")
    print(f"{'='*70}")
    print(f"  Paired subjects: {ipsi.shape[0]}")
    
    # Paired t-test
    spm_test = spm1d.stats.ttest_paired(ipsi, contra)
    spmi = spm_test.inference(alpha=ALPHA, two_tailed=True, interp=True)
    
    # ===== ADDED: Calculate effect size =====
    cohens_d = calculate_cohens_d_paired(ipsi, contra)
    
    # Effect size statistics
    max_d = np.max(np.abs(cohens_d))
    mean_d = np.mean(np.abs(cohens_d))
    
    print(f"  Effect sizes:")
    print(f"    Maximum |Cohen's d|: {max_d:.3f}")
    print(f"    Mean |Cohen's d|: {mean_d:.3f}")
    # ========================================
    
    if spmi.h0reject:
        print(f"  ✓ SIGNIFICANT asymmetry detected")
        print(f"  Number of regions: {len(spmi.clusters)}")
        for i, cluster in enumerate(spmi.clusters, 1):
            # Calculate mean effect size in this significant region
            start_idx = int(cluster.endpoints[0])
            end_idx = int(cluster.endpoints[1])
            region_d = np.mean(np.abs(cohens_d[start_idx:end_idx+1]))
            
            print(f"    Region {i}: p={cluster.P:.4f}, "
                  f"{cluster.endpoints[0]:.1f}%-{cluster.endpoints[1]:.1f}%, "
                  f"mean |d|={region_d:.3f}")
    else:
        print(f"  ✗ No significant asymmetry")
    
    return spm_test, spmi, ipsi, contra, cohens_d


def plot_paired_comparison(ipsi, contra, spmi, cohens_d, cluster_name, signal_name, output_prefix):
    """Plot ipsi vs contra with SPM results AND EFFECT SIZE."""
    
    # ===== MODIFIED: Now 3 panels instead of 2 =====
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 14))
    
    time = np.linspace(0, 100, ipsi.shape[1])
    
    # === PANEL 1: Mean waveforms ===
    mean_ipsi = ipsi.mean(axis=0)
    se_ipsi = ipsi.std(axis=0) / np.sqrt(ipsi.shape[0])
    mean_contra = contra.mean(axis=0)
    se_contra = contra.std(axis=0) / np.sqrt(contra.shape[0])
    
    ax1.plot(time, mean_ipsi, 'b-', linewidth=2.5, label='Ipsilateral', zorder=3)
    ax1.fill_between(time, mean_ipsi - 1.96*se_ipsi, mean_ipsi + 1.96*se_ipsi,
                     color='b', alpha=0.2, zorder=2)
    
    ax1.plot(time, mean_contra, 'r--', linewidth=2.5, label='Contralateral', zorder=3)
    ax1.fill_between(time, mean_contra - 1.96*se_contra, mean_contra + 1.96*se_contra,
                     color='r', alpha=0.2, zorder=2)
    
    # Highlight significant regions
    if spmi.h0reject:
        for cluster in spmi.clusters:
            ax1.axvspan(cluster.endpoints[0], cluster.endpoints[1],
                       alpha=0.3, color='yellow', zorder=1)
    
    ax1.set_ylabel('Angle (degrees)', fontsize=12, fontweight='bold')
    ax1.set_title(f'{cluster_name}: {signal_name}\nIpsilateral vs Contralateral',
                 fontsize=13, fontweight='bold')
    ax1.legend(fontsize=11, loc='best')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim([0, 100])
    
    # === PANEL 2: SPM inference ===
    ax2 = plt.subplot(3, 1, 2)
    spmi.plot()
    spmi.plot_threshold_label()
    if spmi.h0reject:
        spmi.plot_p_values()
    
    ax2.set_xlabel('Gait Cycle (%)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('SPM{t} statistic', fontsize=12)
    ax2.set_title('Statistical Comparison (Paired t-test)', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.axhline(y=0, color='k', linestyle='-', linewidth=0.5, alpha=0.5)
    
    # ===== PANEL 3: EFFECT SIZE (NEW) =====
    ax3.plot(time, cohens_d, 'k-', linewidth=2, label="Cohen's d")
    ax3.axhline(y=0, color='gray', linestyle='-', linewidth=1, alpha=0.5)
    
    # Reference lines for effect size magnitude
    ax3.axhline(y=0.2, color='green', linestyle='--', alpha=0.4, linewidth=1, label='Small (|d|=0.2)')
    ax3.axhline(y=-0.2, color='green', linestyle='--', alpha=0.4, linewidth=1)
    ax3.axhline(y=0.5, color='orange', linestyle='--', alpha=0.4, linewidth=1, label='Medium (|d|=0.5)')
    ax3.axhline(y=-0.5, color='orange', linestyle='--', alpha=0.4, linewidth=1)
    ax3.axhline(y=0.8, color='red', linestyle='--', alpha=0.4, linewidth=1, label='Large (|d|=0.8)')
    ax3.axhline(y=-0.8, color='red', linestyle='--', alpha=0.4, linewidth=1)
    
    # Highlight significant regions
    if spmi.h0reject:
        for cluster in spmi.clusters:
            ax3.axvspan(cluster.endpoints[0], cluster.endpoints[1],
                       alpha=0.3, color='yellow', zorder=1)
    
    ax3.set_xlabel('Gait Cycle (%)', fontsize=12, fontweight='bold')
    ax3.set_ylabel("Cohen's d", fontsize=12, fontweight='bold')
    ax3.set_title('Effect Size (Ipsilateral - Contralateral)', fontsize=12, fontweight='bold')
    ax3.legend(loc='best', fontsize=9)
    ax3.grid(True, alpha=0.3)
    ax3.set_xlim([0, 100])
    # ========================================
    
    plt.tight_layout()
    
    output_path = OUTPUT_DIR / f'{output_prefix}_{cluster_name.replace(" ", "_")}_paired.png'
    plt.savefig(output_path, dpi=DPI, bbox_inches='tight')
    print(f"  Saved: {output_path.name}")
    plt.close()


def calculate_asymmetry_index(ipsi, contra):
    """Calculate asymmetry index at each time point."""
    # AI = |ipsi - contra| / mean × 100
    mean_vals = (ipsi + contra) / 2
    # Avoid division by zero
    ai = np.divide(np.abs(ipsi - contra), mean_vals, 
                   out=np.zeros_like(mean_vals), 
                   where=mean_vals!=0) * 100
    return ai


def plot_asymmetry_comparison(ipsi1, contra1, ipsi2, contra2, signal_name, output_prefix):
    """Compare asymmetry between clusters."""
    
    ai1 = calculate_asymmetry_index(ipsi1, contra1)
    ai2 = calculate_asymmetry_index(ipsi2, contra2)
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    time = np.linspace(0, 100, ai1.shape[1])
    
    # Mean AI for each cluster
    mean_ai1 = ai1.mean(axis=0)
    se_ai1 = ai1.std(axis=0) / np.sqrt(ai1.shape[0])
    mean_ai2 = ai2.mean(axis=0)
    se_ai2 = ai2.std(axis=0) / np.sqrt(ai2.shape[0])
    
    ax.plot(time, mean_ai1, 'darkblue', linewidth=2.5, label='Cluster 1', zorder=3)
    ax.fill_between(time, mean_ai1 - 1.96*se_ai1, mean_ai1 + 1.96*se_ai1,
                    color='darkblue', alpha=0.2, zorder=2)
    
    ax.plot(time, mean_ai2, 'cyan', linewidth=2.5, label='Cluster 2', zorder=3)
    ax.fill_between(time, mean_ai2 - 1.96*se_ai2, mean_ai2 + 1.96*se_ai2,
                    color='cyan', alpha=0.2, zorder=2)
    
    ax.set_xlabel('Gait Cycle (%)', fontsize=13, fontweight='bold')
    ax.set_ylabel('Asymmetry Index (%)', fontsize=13, fontweight='bold')
    ax.set_title(f'{signal_name}\nAsymmetry Index (Ipsi vs Contra) by Cluster',
                fontsize=14, fontweight='bold')
    ax.legend(fontsize=12, loc='best')
    ax.grid(True, alpha=0.3)
    ax.set_xlim([0, 100])
    ax.axhline(y=10, color='orange', linestyle='--', alpha=0.5, label='10% threshold')
    
    plt.tight_layout()
    
    output_path = OUTPUT_DIR / f'{output_prefix}_asymmetry_index.png'
    plt.savefig(output_path, dpi=DPI, bbox_inches='tight')
    print(f"  Saved: {output_path.name}")
    plt.close()
    
    # Return mean AIs for summary
    return mean_ai1.mean(), mean_ai2.mean()


def plot_effect_size_comparison(cohens_d1, cohens_d2, signal_name, output_prefix):
    """
    Compare effect sizes between clusters.
    NEW FUNCTION to visualize effect size differences.
    """
    fig, ax = plt.subplots(figsize=(12, 6))
    
    time = np.linspace(0, 100, len(cohens_d1))
    
    ax.plot(time, cohens_d1, 'darkblue', linewidth=2.5, label='Cluster 1', zorder=3)
    ax.plot(time, cohens_d2, 'cyan', linewidth=2.5, label='Cluster 2', zorder=3)
    
    ax.axhline(y=0, color='gray', linestyle='-', linewidth=1, alpha=0.5)
    
    # Reference lines
    for val, color, label in [(0.2, 'green', 'Small'), 
                               (0.5, 'orange', 'Medium'), 
                               (0.8, 'red', 'Large')]:
        ax.axhline(y=val, color=color, linestyle='--', alpha=0.3, linewidth=1)
        ax.axhline(y=-val, color=color, linestyle='--', alpha=0.3, linewidth=1)
    
    ax.set_xlabel('Gait Cycle (%)', fontsize=13, fontweight='bold')
    ax.set_ylabel("Cohen's d", fontsize=13, fontweight='bold')
    ax.set_title(f"{signal_name}\nEffect Size Comparison: Ipsilateral vs Contralateral",
                fontsize=14, fontweight='bold')
    ax.legend(fontsize=12, loc='best')
    ax.grid(True, alpha=0.3)
    ax.set_xlim([0, 100])
    
    # Add text annotations for reference
    ax.text(102, 0.2, 'Small', fontsize=9, color='green', va='center')
    ax.text(102, 0.5, 'Medium', fontsize=9, color='orange', va='center')
    ax.text(102, 0.8, 'Large', fontsize=9, color='red', va='center')
    
    plt.tight_layout()
    
    output_path = OUTPUT_DIR / f'{output_prefix}_effect_size_comparison.png'
    plt.savefig(output_path, dpi=DPI, bbox_inches='tight')
    print(f"  Saved: {output_path.name}")
    plt.close()


def analyze_signal_bilateral(base_name, signal_name):
    """Complete bilateral analysis for one signal WITH EFFECT SIZES."""
    
    print(f"\n{'#'*80}")
    print(f"# BILATERAL ANALYSIS: {signal_name}")
    print(f"{'#'*80}")
    
    # Load data
    ipsi_all, contra_all, clusters, subjects = load_paired_data(base_name)
    
    # ===== MODIFIED: Now returns cohens_d =====
    # Cluster 1 paired analysis
    spm1, spmi1, ipsi1, contra1, cohens_d1 = paired_spm_by_cluster(
        ipsi_all, contra_all, clusters, 1, 'Cluster 1'
    )
    plot_paired_comparison(ipsi1, contra1, spmi1, cohens_d1, 'Cluster 1', 
                          signal_name, base_name)
    
    # Cluster 2 paired analysis
    spm2, spmi2, ipsi2, contra2, cohens_d2 = paired_spm_by_cluster(
        ipsi_all, contra_all, clusters, 2, 'Cluster 2'
    )
    plot_paired_comparison(ipsi2, contra2, spmi2, cohens_d2, 'Cluster 2',
                          signal_name, base_name)
    # ========================================
    
    # Asymmetry index comparison
    print("\nCalculating asymmetry indices...")
    mean_ai1, mean_ai2 = plot_asymmetry_comparison(
        ipsi1, contra1, ipsi2, contra2, signal_name, base_name
    )
    
    # ===== ADDED: Effect size comparison plot =====
    print("Creating effect size comparison plot...")
    plot_effect_size_comparison(cohens_d1, cohens_d2, signal_name, base_name)
    # ==============================================
    
    # Compile results
    results = {
        'Signal': signal_name,
        'N_Cluster1': ipsi1.shape[0],
        'N_Cluster2': ipsi2.shape[0],
        'Cluster1_Asymmetric': spmi1.h0reject,
        'Cluster2_Asymmetric': spmi2.h0reject,
        'Cluster1_Mean_AI': round(mean_ai1, 2),
        'Cluster2_Mean_AI': round(mean_ai2, 2),
        # ===== ADDED: Effect size metrics =====
        'Cluster1_Max_EffectSize': round(np.max(np.abs(cohens_d1)), 3),
        'Cluster1_Mean_EffectSize': round(np.mean(np.abs(cohens_d1)), 3),
        'Cluster2_Max_EffectSize': round(np.max(np.abs(cohens_d2)), 3),
        'Cluster2_Mean_EffectSize': round(np.mean(np.abs(cohens_d2)), 3)
        # ======================================
    }
    
    if spmi1.h0reject:
        results['C1_NumRegions'] = len(spmi1.clusters)
        results['C1_MinP'] = min([c.P for c in spmi1.clusters])
        # ===== ADDED: Mean effect size in significant regions =====
        region_effects = []
        for cluster in spmi1.clusters:
            start_idx = int(cluster.endpoints[0])
            end_idx = int(cluster.endpoints[1])
            region_d = np.mean(np.abs(cohens_d1[start_idx:end_idx+1]))
            region_effects.append(region_d)
        results['C1_Mean_SigRegion_EffectSize'] = round(np.mean(region_effects), 3)
        # =========================================================
    
    if spmi2.h0reject:
        results['C2_NumRegions'] = len(spmi2.clusters)
        results['C2_MinP'] = min([c.P for c in spmi2.clusters])
        # ===== ADDED: Mean effect size in significant regions =====
        region_effects = []
        for cluster in spmi2.clusters:
            start_idx = int(cluster.endpoints[0])
            end_idx = int(cluster.endpoints[1])
            region_d = np.mean(np.abs(cohens_d2[start_idx:end_idx+1]))
            region_effects.append(region_d)
        results['C2_Mean_SigRegion_EffectSize'] = round(np.mean(region_effects), 3)
        # =========================================================
    
    return results


def main():
    """Main bilateral analysis WITH EFFECT SIZES."""
    
    print("="*80)
    print("BILATERAL ANALYSIS: Ipsilateral vs Contralateral by Cluster")
    print("WITH EFFECT SIZE CALCULATIONS")
    print("="*80)
    
    signals = [
        ('knee_frontal', 'Knee Frontal Angle'),
        ('knee_sagittal', 'Knee Sagittal Angle'),
        ('hip_frontal', 'Hip Frontal Angle'),
        ('hip_sagittal', 'Hip Sagittal Angle'),
    ]
    
    all_results = []
    for base_name, signal_name in signals:
        try:
            results = analyze_signal_bilateral(base_name, signal_name)
            all_results.append(results)
        except Exception as e:
            print(f"\n✗ Error analyzing {signal_name}: {e}")
            import traceback
            traceback.print_exc()
    
    # Summary
    if all_results:
        summary_df = pd.DataFrame(all_results)
        
        print("\n" + "="*80)
        print("SUMMARY WITH EFFECT SIZES")
        print("="*80 + "\n")
        print(summary_df.to_string(index=False))
        
        summary_path = OUTPUT_DIR / 'bilateral_analysis_summary_with_effect_sizes.csv'
        summary_df.to_csv(summary_path, index=False)
        print(f"\n✓ Summary saved: {summary_path}")
        
        # ===== ADDED: Effect size interpretation guide =====
        print("\n" + "="*80)
        print("EFFECT SIZE INTERPRETATION")
        print("="*80)
        print("\nCohen's d guidelines:")
        print("  Small effect:  |d| = 0.2")
        print("  Medium effect: |d| = 0.5")
        print("  Large effect:  |d| = 0.8")
        print("\nCluster 1 asymmetry:")
        for _, row in summary_df.iterrows():
            if row['Cluster1_Asymmetric']:
                print(f"  {row['Signal']}: Max |d|={row['Cluster1_Max_EffectSize']}, "
                      f"Mean |d|={row['Cluster1_Mean_EffectSize']}")
                if 'C1_Mean_SigRegion_EffectSize' in row:
                    print(f"    Significant regions mean |d|={row['C1_Mean_SigRegion_EffectSize']}")
        
        print("\nCluster 2 asymmetry:")
        has_c2_asymmetry = False
        for _, row in summary_df.iterrows():
            if row['Cluster2_Asymmetric']:
                has_c2_asymmetry = True
                print(f"  {row['Signal']}: Max |d|={row['Cluster2_Max_EffectSize']}, "
                      f"Mean |d|={row['Cluster2_Mean_EffectSize']}")
                if 'C2_Mean_SigRegion_EffectSize' in row:
                    print(f"    Significant regions mean |d|={row['C2_Mean_SigRegion_EffectSize']}")
        
        if not has_c2_asymmetry:
            print("  No significant asymmetry detected")
        # ===================================================
    
    print("\n" + "="*80)
    print("BILATERAL ANALYSIS COMPLETE")
    print("="*80 + "\n")


if __name__ == "__main__":
    main()