#!/usr/bin/env python3
"""
Plot RMSD from GROMACS analysis to assess protein relaxation
Shows plateau regions indicating structural stabilization
"""

import matplotlib.pyplot as plt
import numpy as np
import sys
import os

def read_xvg_file(filename):
    """Read GROMACS .xvg file and extract time and RMSD data"""
    times = []
    rmsds = []
    
    with open(filename, 'r') as f:
        for line in f:
            # Skip comments and headers (lines starting with @ or #)
            if line.startswith('@') or line.startswith('#'):
                continue
            
            # Parse data lines
            parts = line.strip().split()
            if len(parts) >= 2:
                try:
                    time = float(parts[0])  # Time in ns
                    rmsd = float(parts[1])  # RMSD in nm
                    times.append(time)
                    rmsds.append(rmsd * 10)  # Convert nm to Å for clarity
                except ValueError:
                    continue
    
    return np.array(times), np.array(rmsds)

def analyze_plateau(times, rmsds, window_size=1000):
    """Simple plateau detection using moving average slope"""
    if len(rmsds) < window_size:
        return None, None
    
    # Calculate moving average slope
    slopes = []
    for i in range(window_size, len(rmsds)):
        window_data = rmsds[i-window_size:i]
        window_times = times[i-window_size:i]
        
        # Linear regression slope
        slope = np.polyfit(window_times, window_data, 1)[0]
        slopes.append(abs(slope))
    
    # Find where slope becomes very small (plateau region)
    threshold = np.percentile(slopes, 25)  # Bottom 25% of slopes
    plateau_start_idx = None
    
    for i, slope in enumerate(slopes):
        if slope < threshold:
            plateau_start_idx = window_size + i
            break
    
    if plateau_start_idx:
        plateau_start_time = times[plateau_start_idx]
        return plateau_start_time, threshold
    
    return None, threshold

def plot_rmsd(filename):
    """Create RMSD plot with plateau analysis"""
    if not os.path.exists(filename):
        print(f"❌ RMSD file not found: {filename}")
        print("   Make sure the MD simulation and analysis have completed!")
        return
    
    # Read data
    times, rmsds = read_xvg_file(filename)
    
    if len(times) == 0:
        print(f"❌ No data found in {filename}")
        return
    
    # Analyze plateau
    plateau_time, threshold = analyze_plateau(times, rmsds)
    
    # Create plot
    plt.figure(figsize=(12, 8))
    
    # Main RMSD plot
    plt.subplot(2, 1, 1)
    plt.plot(times, rmsds, 'b-', linewidth=1.5, alpha=0.8)
    plt.xlabel('Time (ns)')
    plt.ylabel('RMSD (Å)')
    plt.title('Protein Backbone RMSD vs Time\n(Assessing Apo Protein Relaxation)', fontsize=14, pad=20)
    plt.grid(True, alpha=0.3)
    
    # Mark plateau region if detected
    if plateau_time:
        plt.axvline(x=plateau_time, color='red', linestyle='--', alpha=0.7, 
                   label=f'Plateau starts ~{plateau_time:.1f} ns')
        plt.fill_betweenx(plt.ylim(), plateau_time, times[-1], alpha=0.2, color='green',
                         label='Plateau region')
        plt.legend()
    
    # Statistics box
    final_rmsd = rmsds[-1]
    max_rmsd = np.max(rmsds)
    avg_rmsd = np.mean(rmsds)
    
    stats_text = f"""Statistics:
    Final RMSD: {final_rmsd:.2f} Å
    Maximum RMSD: {max_rmsd:.2f} Å  
    Average RMSD: {avg_rmsd:.2f} Å
    Simulation Time: {times[-1]:.1f} ns"""
    
    plt.text(0.02, 0.98, stats_text, transform=plt.gca().transAxes, 
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    # Moving average plot (smoothed trend)
    plt.subplot(2, 1, 2)
    window = min(100, len(rmsds)//10)  # Adaptive window size
    if window > 1:
        smoothed = np.convolve(rmsds, np.ones(window)/window, mode='valid')
        smooth_times = times[window-1:]
        plt.plot(smooth_times, smoothed, 'r-', linewidth=2, label=f'Moving average (window={window})')
        plt.plot(times, rmsds, 'b-', alpha=0.3, linewidth=0.5, label='Raw data')
    else:
        plt.plot(times, rmsds, 'b-', linewidth=1.5)
    
    plt.xlabel('Time (ns)')
    plt.ylabel('RMSD (Å)')
    plt.title('Smoothed RMSD Trend (Plateau Assessment)', fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Assessment text
    assessment = ""
    if plateau_time:
        if plateau_time < times[-1] * 0.5:  # Plateau in first half
            assessment = "✅ GOOD: Protein appears well-relaxed with early plateau"
        elif plateau_time < times[-1] * 0.8:  # Plateau in latter part
            assessment = "⚠️  MODERATE: Protein relaxing, plateau in later part"
        else:
            assessment = "❌ POOR: Late plateau - may need longer simulation"
    else:
        if rmsds[-1] - rmsds[len(rmsds)//2] < 0.5:  # Small change in second half
            assessment = "⚠️  UNCLEAR: No clear plateau, but RMSD seems stable"
        else:
            assessment = "❌ POOR: No plateau detected - protein still changing"
    
    plt.figtext(0.02, 0.02, f"Assessment: {assessment}", fontsize=12, 
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    plt.tight_layout()
    
    # Save plot
    output_file = filename.replace('.xvg', '_plot.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✅ RMSD plot saved: {output_file}")
    
    # Show plot
    plt.show()
    
    # Print summary
    print(f"\n📊 RMSD Analysis Summary:")
    print(f"   Simulation time: {times[-1]:.1f} ns")
    print(f"   Final RMSD: {final_rmsd:.2f} Å")
    print(f"   {assessment}")
    
    if plateau_time:
        relaxation_time = plateau_time
        print(f"   Estimated relaxation time: ~{relaxation_time:.1f} ns")
        
        # Recommendation
        if relaxation_time < times[-1] * 0.3:
            print(f"   💡 Recommendation: Protein well-relaxed, good for docking!")
        else:
            print(f"   💡 Recommendation: Consider longer simulation for better relaxation")

def main():
    """Main function - can specify file or use default"""
    rmsd_file = "analysis/rmsd_apo.xvg"
    
    if len(sys.argv) > 1:
        rmsd_file = sys.argv[1]
    
    print(f"🔍 Analyzing RMSD data: {rmsd_file}")
    plot_rmsd(rmsd_file)

if __name__ == "__main__":
    main() 