#!/usr/bin/env python3
"""
Plot RMSF from GROMACS analysis to show per-residue flexibility
Identifies flexible and rigid regions of the protein
"""

import matplotlib.pyplot as plt
import numpy as np
import argparse
import os

def read_xvg_file(filename):
    """Read GROMACS .xvg file and extract residue and RMSF data"""
    residues = []
    rmsfs = []
    
    with open(filename, 'r') as f:
        for line in f:
            # Skip comments and headers (lines starting with @ or #)
            if line.startswith('@') or line.startswith('#'):
                continue
            
            # Parse data lines
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    residue = int(float(parts[0]))  # Residue number
                    rmsf = float(parts[1]) * 10     # Convert nm to Å
                    residues.append(residue)
                    rmsfs.append(rmsf)
                except ValueError:
                    continue
    
    return np.array(residues), np.array(rmsfs)

def plot_rmsf(residues, rmsfs, output_file, title="RMSF Analysis"):
    """Create RMSF plot with flexibility regions highlighted"""
    
    plt.figure(figsize=(12, 8))
    
    # Main RMSF plot
    plt.plot(residues, rmsfs, 'b-', linewidth=1.5, alpha=0.8)
    plt.fill_between(residues, rmsfs, alpha=0.3, color='lightblue')
    
    # Highlight flexible regions (RMSF > 2.0 Å)
    flexible_mask = rmsfs > 2.0
    if np.any(flexible_mask):
        plt.fill_between(residues, 0, rmsfs, where=flexible_mask, 
                        color='red', alpha=0.3, label='Flexible regions (>2.0 Å)')
    
    # Highlight very flexible regions (RMSF > 3.0 Å)  
    very_flexible_mask = rmsfs > 3.0
    if np.any(very_flexible_mask):
        plt.fill_between(residues, 0, rmsfs, where=very_flexible_mask,
                        color='darkred', alpha=0.5, label='Very flexible regions (>3.0 Å)')
    
    # Add average RMSF line
    avg_rmsf = np.mean(rmsfs)
    plt.axhline(y=avg_rmsf, color='green', linestyle='--', 
                label=f'Average RMSF: {avg_rmsf:.2f} Å')
    
    # Formatting
    plt.xlabel('Residue Number', fontsize=12)
    plt.ylabel('RMSF (Å)', fontsize=12)
    plt.title(title, fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Set reasonable y-axis limits
    plt.ylim(0, max(np.max(rmsfs) * 1.1, 4.0))
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    return avg_rmsf, np.max(rmsfs)

def main():
    parser = argparse.ArgumentParser(description='Plot RMSF analysis from GROMACS')
    parser.add_argument('--input', default='analysis/rmsf_apo.xvg', help='Input RMSF XVG file')
    parser.add_argument('--output', default='analysis/rmsf_apo_plot.png', help='Output plot file')
    parser.add_argument('--title', default='RMSF Analysis - Apo Protein', help='Plot title')
    
    args = parser.parse_args()
    
    if not os.path.exists(args.input):
        print(f"❌ Error: Input file {args.input} not found!")
        return 1
    
    print(f"🔍 Analyzing RMSF data: {args.input}")
    
    # Read data
    residues, rmsfs = read_xvg_file(args.input)
    
    if len(residues) == 0:
        print("❌ Error: No data found in XVG file!")
        return 1
    
    # Create plot
    avg_rmsf, max_rmsf = plot_rmsf(residues, rmsfs, args.output, args.title)
    
    print(f"✅ RMSF plot saved: {args.output}")
    print(f"\n📊 RMSF Analysis Summary:")
    print(f"   Number of residues: {len(residues)}")
    print(f"   Average RMSF: {avg_rmsf:.2f} Å")
    print(f"   Maximum RMSF: {max_rmsf:.2f} Å")
    
    # Identify flexible regions
    flexible_residues = residues[rmsfs > 2.0]
    very_flexible_residues = residues[rmsfs > 3.0]
    
    if len(flexible_residues) > 0:
        print(f"   Flexible regions (>2.0 Å): {len(flexible_residues)} residues")
        print(f"   Most flexible residue: {residues[np.argmax(rmsfs)]} (RMSF: {max_rmsf:.2f} Å)")
    
    if len(very_flexible_residues) > 0:
        print(f"   ⚠️  Very flexible regions (>3.0 Å): {len(very_flexible_residues)} residues")
    else:
        print(f"   ✅ No extremely flexible regions detected")
    
    return 0

if __name__ == "__main__":
    exit(main()) 