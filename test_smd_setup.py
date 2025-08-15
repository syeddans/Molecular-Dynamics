#!/usr/bin/env python3
"""
Test the new SMD setup with ActiveSite group using existing run3 data.
Run a very short simulation to verify the pulling works correctly.
"""

import os
import subprocess
import sys

def run_command(cmd, description=""):
    """Run a command and handle errors."""
    print(f"\n=== {description} ===")
    print(f"Running: {cmd}")
    
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"ERROR: {description} failed")
        print(f"STDOUT: {result.stdout}")
        print(f"STDERR: {result.stderr}")
        return False
    else:
        print(f"SUCCESS: {description} completed")
        if result.stdout:
            print(f"Output: {result.stdout}")
        return True

def main():
    print("Testing new SMD setup with ActiveSite group...")
    
    # Check if we have the necessary files
    required_files = [
        "run3/best_pose/npt/npt.gro",
        "run3/best_pose/npt/npt.cpt", 
        "run3/best_pose/ions/ion.top",
        "mdp/pull_short.mdp",
        "output/autodock_vina/active_site_residues_updated.txt"
    ]
    
    missing_files = []
    for file in required_files:
        if not os.path.exists(file):
            missing_files.append(file)
    
    if missing_files:
        print("ERROR: Missing required files:")
        for file in missing_files:
            print(f"  - {file}")
        return False
    
    # Create test directory
    test_dir = "test_smd"
    os.makedirs(test_dir, exist_ok=True)
    
    # Step 1: Create topology with ActiveSite group
    print("\n" + "="*50)
    print("STEP 1: Creating topology with ActiveSite group")
    print("="*50)
    
    cmd = f"""python scripts/add_active_site_group.py \
        --topology run3/best_pose/ions/ion.top \
        --residues output/autodock_vina/active_site_residues_updated.txt \
        --output {test_dir}/test_smd.top"""
    
    if not run_command(cmd, "Create ActiveSite topology"):
        return False
    
    # Step 2: Create test MDP with very short simulation
    print("\n" + "="*50)
    print("STEP 2: Creating test MDP file (100 steps = 0.1 ps)")
    print("="*50)
    
    # Read original pull_short.mdp and modify for very short test
    with open("mdp/pull_short.mdp", "r") as f:
        mdp_content = f.read()
    
    # Replace nsteps with a very small number for testing
    mdp_content = mdp_content.replace("nsteps                  = 2500000", "nsteps                  = 100")
    mdp_content = mdp_content.replace("nstxout                 = 500", "nstxout                 = 10")
    mdp_content = mdp_content.replace("nstvout                 = 500", "nstvout                 = 10")
    mdp_content = mdp_content.replace("nstxout-compressed      = 500", "nstxout-compressed      = 10")
    
    with open(f"{test_dir}/test_pull.mdp", "w") as f:
        f.write(mdp_content)
    
    print("Created test MDP with 100 steps (0.1 ps)")
    
    # Step 3: Create TPR file
    print("\n" + "="*50)
    print("STEP 3: Creating TPR file")
    print("="*50)
    
    cmd = f"""gmx grompp -f {test_dir}/test_pull.mdp \
        -c run3/best_pose/npt/npt.gro \
        -t run3/best_pose/npt/npt.cpt \
        -p {test_dir}/test_smd.top \
        -o {test_dir}/test_smd.tpr \
        -maxwarn 1000"""
    
    if not run_command(cmd, "Create TPR file"):
        return False
    
    # Step 4: Run very short test simulation
    print("\n" + "="*50)
    print("STEP 4: Running test SMD simulation (0.1 ps)")
    print("="*50)
    
    cmd = f"""cd {test_dir} && gmx mdrun -v -deffnm test_smd -ntmpi 1 -ntomp 2"""
    
    if not run_command(cmd, "Run test SMD simulation"):
        return False
    
    # Step 5: Check results
    print("\n" + "="*50)
    print("STEP 5: Analyzing results")
    print("="*50)
    
    # Check if pull files were created
    pull_files = [f"{test_dir}/test_smd_pullf.xvg", f"{test_dir}/test_smd_pullx.xvg"]
    
    for pull_file in pull_files:
        if os.path.exists(pull_file):
            print(f"✓ Found: {pull_file}")
            
            # Read first few lines to check data
            with open(pull_file, "r") as f:
                lines = f.readlines()
            
            data_lines = [line for line in lines if not line.startswith("#") and not line.startswith("@")]
            if len(data_lines) >= 2:
                first_line = data_lines[0].strip().split()
                last_line = data_lines[-1].strip().split()
                
                if "pullf" in pull_file:
                    print(f"  Force - Start: {first_line[1]} kJ/mol/nm, End: {last_line[1]} kJ/mol/nm")
                else:
                    print(f"  Position - Start: {first_line[1]} nm, End: {last_line[1]} nm")
                    
                    # Check if position changed (indicating pulling is working)
                    try:
                        start_pos = float(first_line[1])
                        end_pos = float(last_line[1])
                        change = end_pos - start_pos
                        print(f"  Position change: {change:.6f} nm")
                        
                        if abs(change) > 0.0001:  # Some movement detected
                            print("  ✓ Ligand position is changing - pulling is working!")
                        else:
                            print("  ⚠ Very small position change - may need longer simulation")
                    except ValueError:
                        print("  ⚠ Could not parse position values")
        else:
            print(f"✗ Missing: {pull_file}")
    
    print("\n" + "="*50)
    print("TEST COMPLETED")
    print("="*50)
    print(f"Test files created in: {test_dir}/")
    print("Check the pull files to verify the SMD setup is working correctly.")
    
    return True

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
