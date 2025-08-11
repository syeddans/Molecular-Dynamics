#!/usr/bin/env python3
"""
Setup AMBER ligand complex for GROMACS simulation.
This script handles only AMBER-specific ligand parameterization and conversion.
GROMACS commands are handled in the Snakefile.
"""

import sys
import argparse
import os
import subprocess

def run_command(cmd, description):
    """Run a command and handle errors."""
    print(f"Running: {description}")
    print(f"Command: {cmd}")
    
    try:
        result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        print(f"✓ {description} completed successfully")
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"❌ {description} failed:")
        print(f"Error: {e.stderr}")
        sys.exit(1)

def copy_file(src, dst, description):
    """Copy file using Linux cp command."""
    run_command(f"cp {src} {dst}", description)

def process_ligand_parameterization(ligand_sdf, output_dir):
    """
    Process ligand parameterization using AMBER/Antechamber GAFF.
    This is the only AMBER-specific part that can't be done with GROMACS.
    """
    print("=== Step 1: Ligand Parameterization with AMBER/Antechamber ===")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Get absolute paths
    abs_ligand_sdf = os.path.abspath(ligand_sdf)
    abs_output_dir = os.path.abspath(output_dir)
    
    # Save current directory
    original_cwd = os.getcwd()
    
    try:
        # Change to output directory to contain Antechamber's temporary files
        os.chdir(abs_output_dir)
        
        # Copy input file using Linux command
        input_file = "input_ligand.sdf"
        copy_file(abs_ligand_sdf, input_file, "Copying ligand input file")
        
        # Generate MOL2 with Gasteiger charges
        mol2_file = "ligand.mol2"
        run_command(f"antechamber -i {input_file} -fi sdf -o {mol2_file} -fo mol2 -c gas -s 2 -nc 0", 
                    "Generating MOL2 with Gasteiger charges")
        
        # Generate GAFF parameters
        gaff_mol2 = "ligand_gaff.mol2"
        run_command(f"antechamber -i {mol2_file} -fi mol2 -o {gaff_mol2} -fo mol2 -c gas -s 2", 
                    "Generating GAFF parameters")
        
        # Generate additional parameters
        frcmod_file = "ligand.frcmod"
        run_command(f"parmchk2 -i {gaff_mol2} -f mol2 -o {frcmod_file}", 
                    "Generating additional force field parameters")
        
        # Clean up Antechamber intermediate files
        cleanup_files = [
            "ANTECHAMBER_*.AC*", "ATOMTYPE.INF", "NEWPDB.PDB", "PREP.INF",
            "leap.log", "tleap.log", "*.amb2gmx"
        ]
        for pattern in cleanup_files:
            # Use subprocess directly for cleanup to avoid verbose output
            subprocess.run(f"rm -rf {pattern}", shell=True, capture_output=True)
        
        print("✓ Ligand parameterization complete!")
        return os.path.join(abs_output_dir, mol2_file), os.path.join(abs_output_dir, gaff_mol2), os.path.join(abs_output_dir, frcmod_file)
        
    finally:
        # Always return to original directory
        os.chdir(original_cwd)

def convert_to_gromacs_format(gaff_mol2, frcmod_file, output_dir):
    """
    Convert AMBER parameters to GROMACS format.
    This is the conversion step that can't be done with GROMACS tools.
    """
    print("=== Step 2: Converting to GROMACS Format ===")
    
    # Get absolute paths
    abs_output_dir = os.path.abspath(output_dir)
    original_cwd = os.getcwd()
    
    try:
        # Change to output directory
        os.chdir(abs_output_dir)
        
        # Get relative paths within output directory
        gaff_mol2_name = os.path.basename(gaff_mol2)
        frcmod_name = os.path.basename(frcmod_file)
        
        # Create AMBER topology files
        prep_file = "ligand.prep"
        ac_file = "ligand.ac"
        
        run_command(f"antechamber -i {gaff_mol2_name} -fi mol2 -o {prep_file} -fo prepi", 
                    "Creating AMBER prep file")
        run_command(f"antechamber -i {gaff_mol2_name} -fi mol2 -o {ac_file} -fo ac", 
                    "Creating AMBER ac file")
        
        # Create tleap input file
        leap_in = "leap.in"
        leap_content = f"""source leaprc.gaff
loadAmberParams {frcmod_name}
loadAmberPrep {prep_file}
mol = loadmol2 {gaff_mol2_name}
saveAmberParm mol ligand.prmtop ligand.rst7
quit
"""
        
        with open(leap_in, 'w') as f:
            f.write(leap_content)
        
        # Run tleap
        run_command(f"tleap -f {leap_in}", "Running tleap to create AMBER topology")
        
        # Convert to GROMACS using parmed (need absolute path to script)
        convert_script = os.path.join(original_cwd, "scripts/convert_to_gromacs.py")
        run_command(f"python {convert_script} --prmtop ligand.prmtop --rst7 ligand.rst7 --output_prefix ligand", 
                    "Converting AMBER files to GROMACS format")
        
        # Clean up additional intermediate files
        cleanup_files = [
            "ANTECHAMBER_*.AC*", "ATOMTYPE.INF", "NEWPDB.PDB", "PREP.INF",
            "leap.log", "tleap.log", "*.amb2gmx"
        ]
        for pattern in cleanup_files:
            # Use subprocess directly for cleanup to avoid verbose output
            subprocess.run(f"rm -rf {pattern}", shell=True, capture_output=True)
        
        print("✓ GROMACS conversion complete!")
        return os.path.join(abs_output_dir, "ligand.gro"), os.path.join(abs_output_dir, "ligand.itp")
        
    finally:
        # Always return to original directory
        os.chdir(original_cwd)


def main():
    parser = argparse.ArgumentParser(description='Setup AMBER ligand for GROMACS simulation')
    parser.add_argument('--ligand', required=True, help='Input ligand SDF file')
    parser.add_argument('--output-dir', default='merged', help='Output directory')
    
    args = parser.parse_args()
    
    print("=== AMBER Ligand Setup ===")
    print("Handling AMBER-specific ligand parameterization")
    print("GROMACS commands will be handled in Snakefile")
    
    # Step 1: Ligand parameterization (AMBER-specific)
    mol2_file, gaff_mol2, frcmod_file = process_ligand_parameterization(args.ligand, args.output_dir)
    
    # Step 2: Convert to GROMACS format (AMBER-specific)
    ligand_gro, ligand_itp = convert_to_gromacs_format(gaff_mol2, frcmod_file, args.output_dir)
    
    print("\n=== SUCCESS ===")
    print(f"Ligand GRO file: {ligand_gro}")
    print(f"Ligand ITP file: {ligand_itp}")
    print("Ready for GROMACS processing in Snakefile!")

if __name__ == "__main__":
    main() 