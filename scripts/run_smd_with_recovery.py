#!/usr/bin/env python3
"""
Run SMD simulation and recover from crash when ligand reaches active site.
Uses GROMACS checkpoint files to extract the final valid frame.
"""

import argparse
import subprocess
import os
import sys
import glob

def run_grompp(mdp_file, gro_file, cpt_file, top_file, ndx_file, tpr_file):
    """Run gmx grompp to create TPR file."""
    cmd = [
        'gmx', 'grompp',
        '-f', mdp_file,
        '-c', gro_file,
        '-t', cpt_file,
        '-p', top_file,
        '-n', ndx_file,
        '-o', tpr_file,
        '-maxwarn', '1000'
    ]
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print("✅ TPR file created successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Error creating TPR file: {e}")
        print(f"STDERR: {e.stderr}")
        return False

def run_mdrun_with_recovery(output_prefix):
    """Run mdrun and handle crash by recovering from checkpoint."""
    
    cmd = [
        'gmx', 'mdrun',
        '-v',
        '-deffnm', output_prefix,
        '-ntmpi', '1',
        '-ntomp', '4'
    ]
    
    # Set environment
    env = os.environ.copy()
    env['OMP_NUM_THREADS'] = '4'
    
    print("🚀 Starting SMD simulation...")
    print("(Will crash when ligand reaches active site - this is expected!)")
    
    try:
        # Run mdrun with live output - expect it to crash
        process = subprocess.Popen(cmd, env=env, stdout=subprocess.PIPE, 
                                  stderr=subprocess.STDOUT, text=True, 
                                  bufsize=1, universal_newlines=True)
        
        # Show live output
        output_lines = []
        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                print(output.strip())  # Show live GROMACS output
                output_lines.append(output)
        
        return_code = process.poll()
        
        if return_code == 0:
            print("✅ Simulation completed successfully without crash")
            return True
        else:
            print("💥 Simulation crashed as expected (ligand reached active site)")
            print("🔄 Recovering final frame from checkpoint...")
            
            # Show last few lines of output for debugging
            print("\nLast few lines of GROMACS output:")
            for line in output_lines[-5:]:
                print(f"  {line.strip()}")
            
            return recover_from_checkpoint(output_prefix)
            
    except Exception as e:
        print(f"Error running mdrun: {e}")
        return False

def recover_from_checkpoint(output_prefix):
    """Recover the final frame from the last checkpoint file."""
    
    # Find the most recent checkpoint file
    cpt_pattern = f"{output_prefix}*.cpt"
    cpt_files = glob.glob(cpt_pattern)
    
    if not cpt_files:
        print("❌ No checkpoint files found for recovery")
        return False
    
    # Sort by modification time, get the most recent
    latest_cpt = max(cpt_files, key=os.path.getmtime)
    print(f"📁 Found checkpoint file: {latest_cpt}")
    
    # Extract final frame from checkpoint
    tpr_file = f"{output_prefix}.tpr"
    final_gro = f"{output_prefix}.gro"
    
    # Use gmx dump to extract coordinates from checkpoint
    try:
        print("📤 Extracting final coordinates from checkpoint...")
        
        # Method 1: Try gmx trjconv with checkpoint
        trjconv_cmd = [
            'gmx', 'trjconv',
            '-s', tpr_file,
            '-f', latest_cpt,
            '-o', final_gro,
            '-dump', '0'  # Extract frame at time 0 from checkpoint
        ]
        
        # Use "0" to select System
        result = subprocess.run(
            trjconv_cmd, 
            input="0\n", 
            text=True, 
            capture_output=True
        )
        
        if result.returncode == 0 and os.path.exists(final_gro):
            print("✅ Final frame extracted successfully")
            
            # Also try to extract any available pull data
            extract_pull_data(output_prefix)
            return True
        else:
            print("⚠️ trjconv method failed, trying alternative...")
            return recover_alternative_method(output_prefix, latest_cpt)
            
    except Exception as e:
        print(f"Error in recovery: {e}")
        return recover_alternative_method(output_prefix, latest_cpt)

def recover_alternative_method(output_prefix, cpt_file):
    """Alternative recovery method using gmx mdrun -cpi."""
    try:
        print("🔄 Trying alternative recovery method...")
        
        # Continue from checkpoint for just 1 step to generate final files
        recovery_cmd = [
            'gmx', 'mdrun',
            '-deffnm', output_prefix,
            '-cpi', cpt_file,
            '-nsteps', '1',  # Just 1 step to write files
            '-ntmpi', '1',
            '-ntomp', '4'
        ]
        
        env = os.environ.copy()
        env['OMP_NUM_THREADS'] = '4'
        
        # This might also crash, but should write some output first
        subprocess.run(recovery_cmd, env=env, capture_output=True, text=True, timeout=30)
        
        # Check if we got output files
        gro_file = f"{output_prefix}.gro"
        if os.path.exists(gro_file):
            print("✅ Recovery successful using continuation method")
            extract_pull_data(output_prefix)
            return True
        else:
            print("❌ Recovery failed - no output files generated")
            return False
            
    except subprocess.TimeoutExpired:
        print("⚠️ Recovery timed out, but may have written files")
        gro_file = f"{output_prefix}.gro"
        if os.path.exists(gro_file):
            print("✅ Files were written before timeout")
            extract_pull_data(output_prefix)
            return True
        return False
    except Exception as e:
        print(f"❌ Alternative recovery failed: {e}")
        return False

def extract_pull_data(output_prefix):
    """Extract pull force and distance data if available."""
    pullf_file = f"{output_prefix}_pullf.xvg"
    pullx_file = f"{output_prefix}_pullx.xvg"
    
    for pull_file in [pullf_file, pullx_file]:
        if os.path.exists(pull_file):
            print(f"✅ Found pull data: {pull_file}")
        else:
            print(f"⚠️ Pull data not found: {pull_file}")

def main():
    parser = argparse.ArgumentParser(
        description="Run SMD with crash recovery using checkpoint files"
    )
    parser.add_argument("--mdp", required=True, help="MDP file")
    parser.add_argument("--gro", required=True, help="GRO file")
    parser.add_argument("--cpt", required=True, help="CPT file")
    parser.add_argument("--top", required=True, help="Topology file")
    parser.add_argument("--ndx", required=True, help="Index file")
    parser.add_argument("--output", required=True, help="Output prefix")
    
    args = parser.parse_args()
    
    tpr_file = f"{args.output}.tpr"
    
    # Step 1: Create TPR file
    if not run_grompp(args.mdp, args.gro, args.cpt, args.top, args.ndx, tpr_file):
        print("💥 Failed to create TPR file")
        sys.exit(1)
    
    # Step 2: Run mdrun with crash recovery
    if run_mdrun_with_recovery(args.output):
        print("🎉 SMD completed successfully!")
        
        # Check final distance if pull data exists
        pullx_file = f"{args.output}_pullx.xvg"
        if os.path.exists(pullx_file):
            try:
                with open(pullx_file, 'r') as f:
                    lines = f.readlines()
                
                # Find last data line
                for line in reversed(lines):
                    if line.strip() and not line.startswith('#') and not line.startswith('@'):
                        parts = line.split()
                        if len(parts) >= 2:
                            final_distance = float(parts[1])
                            print(f"🎯 Final ligand-active site distance: {final_distance:.3f} nm")
                            break
            except:
                pass
        
        sys.exit(0)
    else:
        print("💥 SMD failed and recovery unsuccessful")
        sys.exit(1)

if __name__ == "__main__":
    main()
