#!/usr/bin/env python3
"""
Comprehensive SPM Analysis with Healthy Controls
Performs multiple comparisons:
1. OA Cluster 1 vs Cluster 2 (identify cluster characteristics)
2. OA Cluster 1 vs Healthy (clinical relevance)
3. OA Cluster 2 vs Healthy (clinical relevance)  
4. All OA vs Healthy (overall pathology)
"""

import numpy as np
import pandas as pd
import spm1d
import matplotlib.pyplot as plt
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Configuration
DATA_DIR = Path('../data/spm_export')
OUTPUT_DIR = Path('../outputs/spm_results')
OUTPUT_DIR.mkdir(exist_ok=True, parents=True)

ALPHA = 0.05
TWO_TAILED = True
DPI = 300

# Color schemes
COLORS = {
    '1': "#11326E",  # dark Blue
    '2': "#359aa1",  # cyab
    'healthy': "#797f7f",  # grey
    'OA': "#056ded"  # blue
}

COMPARISON_NAMES = {
    'cluster1_vs_cluster2': 'OA Cluster 1 vs OA Cluster 2',
    'cluster1_vs_healthy': 'OA Cluster 1 vs Healthy Controls',
    'cluster2_vs_healthy': 'OA Cluster 2 vs Healthy Controls',
    'all_oa_vs_healthy': 'All OA vs Healthy Controls'
}


def load_comparison_data(filename):
    """Load and prepare comparison data from CSV."""
    filepath = DATA_DIR / filename
    
    if not filepath.exists():
        raise FileNotFoundError(f"Cannot find {filepath}")
    
    data = pd.read_csv(filepath)
    
    # Get unique groups
    groups = data['group'].unique()
    if len(groups) != 2:
        raise ValueError(f"Expected 2 groups, found {len(groups)}: {groups}")
    
    # Extract waveforms (columns 3+ are timepoints)
    group1_waves = data[data['group'] == groups[0]].iloc[:, 3:].values
    group2_waves = data[data['group'] == groups[1]].iloc[:, 3:].values
    
    metadata = {
        'filename': filename,
        'group1': str(groups[0]),
        'group2': str(groups[1]),
        'n_group1': group1_waves.shape[0],
        'n_group2': group2_waves.shape[0],
        'n_timepoints': group1_waves.shape[1]
    }
    
    return group1_waves, group2_waves, metadata


def perform_spm_ttest(group1_waves, group2_waves, alpha=0.05, two_tailed=True):
    """Perform SPM two-sample t-test."""
    spm_test = spm1d.stats.ttest2(group1_waves, group2_waves)
    spmi = spm_test.inference(alpha=alpha, two_tailed=two_tailed, interp=True)
    return spm_test, spmi


def get_comparison_labels(metadata, comparison_type):
    """Get appropriate labels for the comparison."""
    g1, g2 = metadata['group1'], metadata['group2']
    
    if comparison_type == 'cluster1_vs_cluster2':
        label1 = f'OA Cluster {g1}'
        label2 = f'OA Cluster {g2}'
        color1, color2 = COLORS.get(g1, 'blue'), COLORS.get(g2, 'orange')
    elif 'healthy' in [g1, g2]:
        if g1 == 'healthy':
            label1, label2 = 'Healthy', f'OA Cluster {g2}'
            color1, color2 = COLORS['healthy'], COLORS.get(g2, 'red')
        else:
            label1, label2 = f'OA Cluster {g1}', 'Healthy'
            color1, color2 = COLORS.get(g1, 'red'), COLORS['healthy']
    elif g1 == 'OA' or g2 == 'OA':
        if g1 == 'OA':
            label1, label2 = 'All OA', 'Healthy'
            color1, color2 = COLORS['OA'], COLORS['healthy']
        else:
            label1, label2 = 'Healthy', 'All OA'
            color1, color2 = COLORS['healthy'], COLORS['OA']
    else:
        label1, label2 = g1, g2
        color1, color2 = 'blue', 'orange'
    
    return label1, label2, color1, color2


def print_spm_results(spmi, comparison_name, metadata):
    """Print formatted SPM results."""
    print(f"\n{'='*80}")
    print(f"SPM RESULTS: {comparison_name}")
    print(f"{'='*80}")
    print(f"Groups: {metadata['group1']} (n={metadata['n_group1']}) vs "
          f"{metadata['group2']} (n={metadata['n_group2']})")
    
    if spmi.h0reject:
        print(f"\n✓ SIGNIFICANT DIFFERENCES DETECTED (α = {ALPHA})")
        print(f"  Number of significant regions: {len(spmi.clusters)}")
        
        for i, cluster in enumerate(spmi.clusters, 1):
            duration = cluster.endpoints[1] - cluster.endpoints[0]
            print(f"\n  Region {i}:")
            print(f"    P-value: {cluster.P:.4f}")
            print(f"    Gait cycle: {cluster.endpoints[0]:.1f}% - {cluster.endpoints[1]:.1f}%")
            print(f"    Duration: {duration:.1f}% of gait cycle")
            
            # Interpret gait phases
            phase = interpret_gait_phase(cluster.endpoints[0], cluster.endpoints[1])
            if phase:
                print(f"    Phase: {phase}")
    else:
        print(f"\n✗ No significant differences detected (α = {ALPHA})")
    
    print(f"{'='*80}\n")


def interpret_gait_phase(start, end):
    """Interpret which gait phase the significant region falls in."""
    mid = (start + end) / 2
    
    if start < 2:
        return "Initial contact/Loading response"
    elif 2 <= mid < 12:
        return "Loading response"
    elif 12 <= mid < 31:
        return "Mid-stance"
    elif 31 <= mid < 50:
        return "Terminal stance"
    elif 50 <= mid < 60:
        return "Pre-swing"
    elif 60 <= mid < 73:
        return "Initial swing"
    elif 73 <= mid < 87:
        return "Mid-swing"
    elif 87 <= mid <= 100:
        return "Terminal swing"
    else:
        return None


def plot_spm_inference(spmi, comparison_name, signal_name, metadata, output_prefix, 
                       label1, label2):
    """Create SPM inference plot."""
    fig = plt.figure(figsize=(14, 6))
    
    spmi.plot()
    spmi.plot_threshold_label()
    spmi.plot_p_values()
    
    plt.xlabel('Gait Cycle (%)', fontsize=13, fontweight='bold')
    plt.ylabel('SPM{t} statistic', fontsize=13, fontweight='bold')
    plt.title(f'{signal_name}\n{comparison_name}\n({label1} vs {label2})', 
              fontsize=14, fontweight='bold', pad=20)
    plt.grid(True, alpha=0.3)
    plt.axhline(y=0, color='k', linestyle='-', linewidth=0.5, alpha=0.5)
    
    # Add sample sizes
    plt.text(0.02, 0.98, f'n₁={metadata["n_group1"]}, n₂={metadata["n_group2"]}',
             transform=plt.gca().transAxes, fontsize=10, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    
    output_path = OUTPUT_DIR / f'{output_prefix}_spm_inference.png'
    plt.savefig(output_path, dpi=DPI, bbox_inches='tight')
    print(f"  Saved: {output_path.name}")
    plt.close()


def plot_mean_waveforms(group1_waves, group2_waves, spmi, comparison_name, 
                       signal_name, metadata, output_prefix, label1, label2, 
                       color1, color2, ylabel='Angle (degrees)'):
    """Plot mean waveforms with confidence intervals."""
    fig, ax = plt.subplots(figsize=(14, 7))
    
    time = np.linspace(0, 100, group1_waves.shape[1])
    mean1 = group1_waves.mean(axis=0)
    se1 = group1_waves.std(axis=0) / np.sqrt(group1_waves.shape[0])
    mean2 = group2_waves.mean(axis=0)
    se2 = group2_waves.std(axis=0) / np.sqrt(group2_waves.shape[0])
    
    # Plot group 1
    ax.plot(time, mean1, color=color1, linewidth=3, label=label1, zorder=3)
    ax.fill_between(time, mean1 - 1.96*se1, mean1 + 1.96*se1, 
                     color=color1, alpha=0.2, zorder=2)
    
    # Plot group 2
    ax.plot(time, mean2, color=color2, linewidth=3, label=label2, zorder=3)
    ax.fill_between(time, mean2 - 1.96*se2, mean2 + 1.96*se2, 
                     color=color2, alpha=0.2, zorder=2)
    
    # Highlight significant regions
    if spmi.h0reject:
        for i, cluster in enumerate(spmi.clusters):
            start_pct = cluster.endpoints[0]
            end_pct = cluster.endpoints[1]
            
            ax.axvspan(start_pct, end_pct, alpha=0.25, color='gold', zorder=1)
            
            # Add p-value annotation
            mid_pct = (start_pct + end_pct) / 2
            y_pos = ax.get_ylim()[1] * 0.96
            ax.text(mid_pct, y_pos, f'p={cluster.P:.3f}', 
                   ha='center', va='top', fontsize=11, fontweight='bold',
                   bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7, 
                            edgecolor='orange', linewidth=2))
    
    # Add gait phase markers
    ax.axvline(x=60, color='gray', linestyle=':', alpha=0.5, linewidth=1.5)
    y_bottom = ax.get_ylim()[0]
    ax.text(30, y_bottom, 'Stance', ha='center', va='bottom', 
           fontsize=11, color='gray', fontweight='bold', alpha=0.7)
    ax.text(80, y_bottom, 'Swing', ha='center', va='bottom', 
           fontsize=11, color='gray', fontweight='bold', alpha=0.7)
    
    # Formatting
    ax.set_xlabel('Gait Cycle (%)', fontsize=14, fontweight='bold')
    ax.set_ylabel(ylabel, fontsize=14, fontweight='bold')
    ax.set_title(f'{signal_name}\n{comparison_name}\nMean Waveforms with 95% CI', 
                fontsize=15, fontweight='bold', pad=20)
    ax.legend(fontsize=12, loc='best', framealpha=0.9)
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.set_xlim([0, 100])
    
    plt.tight_layout()
    
    output_path = OUTPUT_DIR / f'{output_prefix}_mean_waveforms.png'
    plt.savefig(output_path, dpi=DPI, bbox_inches='tight')
    print(f"  Saved: {output_path.name}")
    plt.close()


def calculate_effect_size(group1_waves, group2_waves, output_prefix, 
                         comparison_name, signal_name):
    """Calculate and plot Cohen's d effect size."""
    mean1 = group1_waves.mean(axis=0)
    mean2 = group2_waves.mean(axis=0)
    std1 = group1_waves.std(axis=0)
    std2 = group2_waves.std(axis=0)
    
    n1 = group1_waves.shape[0]
    n2 = group2_waves.shape[0]
    pooled_std = np.sqrt(((n1 - 1) * std1**2 + (n2 - 1) * std2**2) / (n1 + n2 - 2))
    
    cohens_d = (mean1 - mean2) / pooled_std
    
    fig, ax = plt.subplots(figsize=(14, 6))
    time = np.linspace(0, 100, len(cohens_d))
    
    ax.plot(time, cohens_d, color='black', linewidth=2.5, label="Cohen's d")
    ax.fill_between(time, 0, cohens_d, where=(cohens_d > 0), 
                     alpha=0.3, color='blue', label='Group 1 > Group 2')
    ax.fill_between(time, 0, cohens_d, where=(cohens_d < 0), 
                     alpha=0.3, color='red', label='Group 2 > Group 1')
    
    ax.axhline(y=0, color='gray', linestyle='-', linewidth=1.5, alpha=0.5)
    
    # Effect size reference lines
    for magnitude, value, color in [('Small', 0.2, 'green'), 
                                    ('Medium', 0.5, 'orange'), 
                                    ('Large', 0.8, 'red')]:
        ax.axhline(y=value, color=color, linestyle='--', alpha=0.4, linewidth=1)
        ax.axhline(y=-value, color=color, linestyle='--', alpha=0.4, linewidth=1)
        if magnitude == 'Large':
            ax.text(102, value, magnitude, fontsize=9, color=color, 
                   va='center', fontweight='bold')
    
    ax.set_xlabel('Gait Cycle (%)', fontsize=13, fontweight='bold')
    ax.set_ylabel("Cohen's d Effect Size", fontsize=13, fontweight='bold')
    ax.set_title(f"{signal_name}\n{comparison_name}\nEffect Size Across Gait Cycle", 
                fontsize=14, fontweight='bold', pad=15)
    ax.legend(loc='best', fontsize=11, framealpha=0.9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim([0, 100])
    
    plt.tight_layout()
    
    output_path = OUTPUT_DIR / f'{output_prefix}_effect_size.png'
    plt.savefig(output_path, dpi=DPI, bbox_inches='tight')
    print(f"  Saved: {output_path.name}")
    plt.close()
    
    return cohens_d


def analyze_comparison(filename, signal_name, comparison_type, ylabel='Angle (degrees)'):
    """Complete SPM analysis for a single comparison."""
    print(f"\n{'#'*90}")
    print(f"# {COMPARISON_NAMES.get(comparison_type, comparison_type)}")
    print(f"# Signal: {signal_name}")
    print(f"{'#'*90}")
    
    output_prefix = filename.replace('.csv', '')
    
    try:
        # Load data
        group1_waves, group2_waves, metadata = load_comparison_data(filename)
        
        print(f"\nLoaded: {filename}")
        print(f"  Group 1 ({metadata['group1']}): {metadata['n_group1']} waveforms")
        print(f"  Group 2 ({metadata['group2']}): {metadata['n_group2']} waveforms")
        print(f"  Time points: {metadata['n_timepoints']}")
        
        # Get appropriate labels and colors
        label1, label2, color1, color2 = get_comparison_labels(metadata, comparison_type)
        comparison_name = COMPARISON_NAMES.get(comparison_type, comparison_type)
        
        # Perform SPM test
        spm_test, spmi = perform_spm_ttest(group1_waves, group2_waves, 
                                           alpha=ALPHA, two_tailed=TWO_TAILED)
        
        # Print results
        print_spm_results(spmi, comparison_name, metadata)
        
        # Create plots
        print("\nGenerating plots...")
        plot_spm_inference(spmi, comparison_name, signal_name, metadata, 
                          output_prefix, label1, label2)
        plot_mean_waveforms(group1_waves, group2_waves, spmi, comparison_name,
                           signal_name, metadata, output_prefix, label1, label2,
                           color1, color2, ylabel)
        cohens_d = calculate_effect_size(group1_waves, group2_waves, 
                                         output_prefix, comparison_name, signal_name)
        
        # Compile results
        results = {
            'Signal': signal_name,
            'Comparison': comparison_name,
            'Comparison_Type': comparison_type,
            'Group1': metadata['group1'],
            'Group2': metadata['group2'],
            'N_Group1': metadata['n_group1'],
            'N_Group2': metadata['n_group2'],
            'Significant': spmi.h0reject,
            'N_Significant_Regions': len(spmi.clusters) if spmi.h0reject else 0,
            'Max_Effect_Size': np.max(np.abs(cohens_d)),
            'Mean_Effect_Size': np.mean(np.abs(cohens_d))
        }
        
        if spmi.h0reject:
            for i, cluster in enumerate(spmi.clusters, 1):
                results[f'Region{i}_P'] = cluster.P
                results[f'Region{i}_Start'] = cluster.endpoints[0]
                results[f'Region{i}_End'] = cluster.endpoints[1]
                results[f'Region{i}_Duration'] = cluster.endpoints[1] - cluster.endpoints[0]
        
        print(f"\n✓ Analysis complete")
        return results
        
    except Exception as e:
        print(f"\n✗ Error: {str(e)}")
        return {
            'Signal': signal_name,
            'Comparison': comparison_type,
            'Error': str(e)
        }


def main():
    """Main analysis workflow."""
    print("="*90)
    print("COMPREHENSIVE SPM ANALYSIS WITH HEALTHY CONTROLS")
    print("="*90)
    
    # Define all signals and comparisons
    signals = {
        'Knee Frontal Angle': 'knee_frontal_angle',
        'Knee Sagittal Angle': 'knee_sagittal_angle'
    }
    
    comparison_types = [
        'cluster1_vs_cluster2',
        'cluster1_vs_healthy',
        'cluster2_vs_healthy',
        'all_oa_vs_healthy'
    ]
    
    # Collect all results
    all_results = []
    
    # Analyze each signal × comparison combination
    for signal_name, signal_file in signals.items():
        for comp_type in comparison_types:
            filename = f'spm_{signal_file}_{comp_type}.csv'
            
            if not (DATA_DIR / filename).exists():
                print(f"\n⊗ Skipping {filename} - file not found")
                continue
            
            result = analyze_comparison(filename, signal_name, comp_type)
            all_results.append(result)
    
    # Create comprehensive summary
    print("\n" + "="*90)
    print("CREATING SUMMARY REPORTS")
    print("="*90)
    
    results_df = pd.DataFrame(all_results)
    
    # Save detailed results
    results_df.to_csv(OUTPUT_DIR / 'spm_comprehensive_results.csv', index=False)
    print(f"\n✓ Detailed results saved: spm_comprehensive_results.csv")
    
    # Create summary tables
    if 'Error' not in results_df.columns or results_df['Error'].isna().all():
        # Summary by comparison type
        summary_by_comparison = results_df.groupby('Comparison_Type').agg({
            'Significant': 'sum',
            'Signal': 'count'
        }).rename(columns={'Signal': 'Total_Signals'})
        summary_by_comparison['Proportion_Significant'] = (
            summary_by_comparison['Significant'] / summary_by_comparison['Total_Signals']
        )
        summary_by_comparison.to_csv(OUTPUT_DIR / 'summary_by_comparison.csv')
        print(f"✓ Summary by comparison saved")
        
        # Summary by signal
        summary_by_signal = results_df.groupby('Signal').agg({
            'Significant': 'sum',
            'Comparison': 'count'
        }).rename(columns={'Comparison': 'Total_Comparisons'})
        summary_by_signal['Proportion_Significant'] = (
            summary_by_signal['Significant'] / summary_by_signal['Total_Comparisons']
        )
        summary_by_signal.to_csv(OUTPUT_DIR / 'summary_by_signal.csv')
        print(f"✓ Summary by signal saved")
        
        # Print final summary to console
        print("\n" + "="*90)
        print("FINAL SUMMARY")
        print("="*90)
        print("\nResults by Comparison Type:")
        print(summary_by_comparison.to_string())
        print("\nResults by Signal:")
        print(summary_by_signal.to_string())
        
        # Key findings
        print("\n" + "="*90)
        print("KEY FINDINGS")
        print("="*90)
        
        sig_results = results_df[results_df['Significant'] == True]
        
        if len(sig_results) > 0:
            print(f"\n✓ {len(sig_results)} significant comparisons found")
            
            # Cluster differences
            cluster_comp = sig_results[sig_results['Comparison_Type'] == 'cluster1_vs_cluster2']
            if len(cluster_comp) > 0:
                print(f"\n  OA Cluster 1 vs Cluster 2:")
                print(f"    - {len(cluster_comp)} signals show significant differences")
                print(f"    - Signals: {', '.join(cluster_comp['Signal'].unique())}")
            
            # Clinical relevance
            for cluster_num in ['1', '2']:
                cluster_health = sig_results[sig_results['Comparison_Type'] == f'cluster{cluster_num}_vs_healthy']
                if len(cluster_health) > 0:
                    print(f"\n  OA Cluster {cluster_num} vs Healthy:")
                    print(f"    - {len(cluster_health)} signals differ from healthy")
                    print(f"    - Signals: {', '.join(cluster_health['Signal'].unique())}")
            
            # Overall OA effect
            oa_health = sig_results[sig_results['Comparison_Type'] == 'all_oa_vs_healthy']
            if len(oa_health) > 0:
                print(f"\n  All OA vs Healthy:")
                print(f"    - {len(oa_health)} signals differ from healthy")
                print(f"    - Signals: {', '.join(oa_health['Signal'].unique())}")
        else:
            print("\n✗ No significant differences detected in any comparison")
    
    print("\n" + "="*90)
    print(f"ANALYSIS COMPLETE")
    print(f"All outputs saved to: {OUTPUT_DIR}")
    print("="*90 + "\n")


if __name__ == "__main__":
    main()