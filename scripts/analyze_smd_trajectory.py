#!/usr/bin/env python3
"""
Analyze SMD trajectory to check if ligand is moving properly toward the active site.
This script examines both the pull force and distance data to assess SMD quality.

Usage:
    python scripts/analyze_smd_trajectory.py --pose pose_1 --output analysis_report.txt
    python scripts/analyze_smd_trajectory.py --pose pose_1 --plot  # Creates plots
"""

import argparse
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict

def read_xvg_file(filepath: str) -> Tuple[np.ndarray, np.ndarray]:
    """Read GROMACS XVG file and return time and data arrays."""
    times = []
    values = []
    
    with open(filepath, 'r') as f:
        for line in f:
            # Skip comments and headers
            if line.startswith('#') or line.startswith('@') or not line.strip():
                continue
            
            parts = line.split()
            if len(parts) >= 2:
                try:
                    times.append(float(parts[0]))
                    values.append(float(parts[1]))
                except ValueError:
                    continue
    
    return np.array(times), np.array(values)

def analyze_smd_quality(times: np.ndarray, forces: np.ndarray, distances: np.ndarray) -> Dict:
    """Analyze SMD trajectory quality and ligand movement."""
    analysis = {}
    
    # Basic statistics
    analysis['simulation_time'] = float(times[-1] - times[0])  # ps
    analysis['initial_distance'] = float(distances[0])  # nm
    analysis['final_distance'] = float(distances[-1])  # nm
    analysis['distance_change'] = float(distances[-1] - distances[0])  # nm
    
    # Force statistics
    analysis['mean_force'] = float(np.mean(forces))  # kJ/mol/nm
    analysis['max_force'] = float(np.max(forces))
    analysis['min_force'] = float(np.min(forces))
    analysis['force_std'] = float(np.std(forces))
    
    # Movement analysis
    expected_distance_change = -0.001 * analysis['simulation_time'] / 1000  # nm (rate * time)
    analysis['expected_distance_change'] = expected_distance_change
    analysis['movement_efficiency'] = abs(analysis['distance_change'] / expected_distance_change) if expected_distance_change != 0 else 0
    
    # Check if ligand is moving in the right direction
    analysis['moving_toward_active_site'] = analysis['distance_change'] < 0
    
    # Force trend analysis (should generally become more negative as ligand approaches)
    force_trend = np.polyfit(times, forces, 1)[0]  # Linear slope
    analysis['force_trend'] = float(force_trend)
    analysis['force_trend_correct'] = force_trend < 0  # Should become more negative
    
    # Smoothness check (large force spikes indicate problems)
    force_diff = np.diff(forces)
    analysis['force_smoothness'] = float(np.std(force_diff))
    analysis['large_spikes'] = int(np.sum(np.abs(force_diff) > 3 * np.std(force_diff)))
    
    # Distance smoothness (should be relatively smooth decrease)
    distance_diff = np.diff(distances)
    analysis['distance_smoothness'] = float(np.std(distance_diff))
    
    return analysis

def create_smd_plots(times: np.ndarray, forces: np.ndarray, distances: np.ndarray, 
                    pose: str, output_dir: str):
    """Create diagnostic plots for SMD analysis."""
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
    fig.suptitle(f'SMD Analysis for {pose}', fontsize=16, fontweight='bold')
    
    # Plot 1: Distance vs Time
    ax1.plot(times, distances, 'b-', linewidth=2, label='Distance')
    ax1.set_xlabel('Time (ps)')
    ax1.set_ylabel('Distance (nm)')
    ax1.set_title('Ligand-Active Site Distance vs Time')
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    # Add trend line
    z = np.polyfit(times, distances, 1)
    p = np.poly1d(z)
    ax1.plot(times, p(times), 'r--', alpha=0.8, label=f'Trend: {z[0]:.6f} nm/ps')
    ax1.legend()
    
    # Plot 2: Force vs Time
    ax2.plot(times, forces, 'r-', linewidth=2, label='Pull Force')
    ax2.set_xlabel('Time (ps)')
    ax2.set_ylabel('Force (kJ/mol/nm)')
    ax2.set_title('Pull Force vs Time')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    
    # Add horizontal line at zero
    ax2.axhline(y=0, color='k', linestyle='--', alpha=0.5)
    
    # Plot 3: Force vs Distance
    ax3.scatter(distances, forces, c=times, cmap='viridis', alpha=0.7, s=20)
    ax3.set_xlabel('Distance (nm)')
    ax3.set_ylabel('Force (kJ/mol/nm)')
    ax3.set_title('Force vs Distance (colored by time)')
    ax3.grid(True, alpha=0.3)
    cbar = plt.colorbar(ax3.scatter(distances, forces, c=times, cmap='viridis', alpha=0.7, s=20), ax=ax3)
    cbar.set_label('Time (ps)')
    
    # Plot 4: Force distribution
    ax4.hist(forces, bins=30, alpha=0.7, color='green', edgecolor='black')
    ax4.axvline(np.mean(forces), color='red', linestyle='--', linewidth=2, label=f'Mean: {np.mean(forces):.2f}')
    ax4.set_xlabel('Force (kJ/mol/nm)')
    ax4.set_ylabel('Frequency')
    ax4.set_title('Force Distribution')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot
    plot_file = os.path.join(output_dir, f'{pose}_smd_analysis.png')
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"Saved SMD analysis plot: {plot_file}")
    plt.close()

def generate_report(analysis: Dict, pose: str) -> str:
    """Generate a text report of SMD analysis."""
    
    report = f"""
SMD Analysis Report for {pose}
{'='*50}

Simulation Parameters:
- Total simulation time: {analysis['simulation_time']:.1f} ps
- Expected pull rate: -0.001 nm/ps (toward active site)

Distance Analysis:
- Initial distance: {analysis['initial_distance']:.3f} nm
- Final distance: {analysis['final_distance']:.3f} nm
- Distance change: {analysis['distance_change']:.3f} nm
- Expected change: {analysis['expected_distance_change']:.3f} nm
- Movement efficiency: {analysis['movement_efficiency']:.1%}

Force Analysis:
- Mean force: {analysis['mean_force']:.2f} kJ/mol/nm
- Force range: {analysis['min_force']:.2f} to {analysis['max_force']:.2f} kJ/mol/nm
- Force standard deviation: {analysis['force_std']:.2f} kJ/mol/nm
- Force trend slope: {analysis['force_trend']:.6f} kJ/mol/nm/ps

Quality Assessment:
- Moving toward active site: {'✓' if analysis['moving_toward_active_site'] else '✗'}
- Force trend correct: {'✓' if analysis['force_trend_correct'] else '✗'}
- Force smoothness: {analysis['force_smoothness']:.2f} (lower is better)
- Large force spikes: {analysis['large_spikes']}
- Distance smoothness: {analysis['distance_smoothness']:.6f} (lower is better)

Interpretation:
"""
    
    # Add interpretation
    if analysis['moving_toward_active_site']:
        if analysis['movement_efficiency'] > 0.8:
            report += "✓ GOOD: Ligand is moving efficiently toward the active site.\n"
        elif analysis['movement_efficiency'] > 0.5:
            report += "⚠ MODERATE: Ligand is moving toward active site but with some resistance.\n"
        else:
            report += "⚠ POOR: Ligand movement is significantly impeded.\n"
    else:
        report += "✗ ERROR: Ligand is moving AWAY from the active site! Check setup.\n"
    
    if analysis['mean_force'] < -5:
        report += "✓ Strong attractive forces toward active site.\n"
    elif analysis['mean_force'] < 0:
        report += "⚠ Moderate attractive forces toward active site.\n"
    else:
        report += "⚠ Repulsive forces - ligand may be encountering steric clashes.\n"
    
    if analysis['large_spikes'] == 0:
        report += "✓ Smooth force profile - good simulation stability.\n"
    elif analysis['large_spikes'] < 5:
        report += "⚠ Few force spikes detected - mostly stable.\n"
    else:
        report += "✗ Many force spikes - simulation may be unstable.\n"
    
    return report

def main():
    parser = argparse.ArgumentParser(
        description="Analyze SMD trajectory for ligand movement toward active site",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Analyze specific pose and generate report
    python scripts/analyze_smd_trajectory.py --pose pose_1 --output pose_1_smd_report.txt
    
    # Analyze with plots
    python scripts/analyze_smd_trajectory.py --pose pose_1 --plot --output pose_1_smd_report.txt
    
    # Analyze all poses
    for i in {1..9}; do
        python scripts/analyze_smd_trajectory.py --pose pose_$i --output pose_${i}_smd_report.txt --plot
    done
        """
    )
    
    parser.add_argument("--pose", required=True, help="Pose name (e.g., pose_1)")
    parser.add_argument("--output", help="Output report file (optional)")
    parser.add_argument("--plot", action="store_true", help="Generate diagnostic plots")
    parser.add_argument("--output-dir", default="output/smd_analysis", help="Directory for output files")
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # File paths
    pose_dir = f"output/pose_runs/{args.pose}/smd"
    pullf_file = os.path.join(pose_dir, "pullf.xvg")
    pullx_file = os.path.join(pose_dir, "pullx.xvg")
    
    # Check if files exist
    for file_path in [pullf_file, pullx_file]:
        if not os.path.exists(file_path):
            print(f"Error: File not found: {file_path}", file=sys.stderr)
            print("Make sure SMD simulation has completed for this pose.", file=sys.stderr)
            sys.exit(1)
    
    # Read data
    try:
        times_f, forces = read_xvg_file(pullf_file)
        times_x, distances = read_xvg_file(pullx_file)
        
        # Use force times as reference (should be the same)
        times = times_f
        
        print(f"Loaded SMD data for {args.pose}:")
        print(f"  - {len(times)} data points")
        print(f"  - Simulation time: {times[-1] - times[0]:.1f} ps")
        print(f"  - Distance change: {distances[-1] - distances[0]:.3f} nm")
        
    except Exception as e:
        print(f"Error reading SMD data: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Analyze trajectory
    analysis = analyze_smd_quality(times, forces, distances)
    
    # Generate report
    report = generate_report(analysis, args.pose)
    print(report)
    
    # Save report if requested
    if args.output:
        with open(args.output, 'w') as f:
            f.write(report)
        print(f"\nReport saved to: {args.output}")
    
    # Generate plots if requested
    if args.plot:
        try:
            create_smd_plots(times, forces, distances, args.pose, args.output_dir)
        except ImportError:
            print("Warning: matplotlib not available for plotting")
        except Exception as e:
            print(f"Warning: Could not create plots: {e}")
    
    # Summary assessment
    print(f"\n{'='*50}")
    print(f"SUMMARY for {args.pose}:")
    if analysis['moving_toward_active_site'] and analysis['movement_efficiency'] > 0.5:
        print("✓ SMD appears to be working correctly")
    else:
        print("⚠ SMD may have issues - review the detailed analysis above")
    print(f"{'='*50}")

if __name__ == "__main__":
    main()
