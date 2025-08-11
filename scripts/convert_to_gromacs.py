#!/usr/bin/env python3
"""
Convert AMBER files to GROMACS format
"""

import argparse
import os
import sys
import subprocess

def main():
    parser = argparse.ArgumentParser(description='Convert AMBER files to GROMACS format')
    parser.add_argument('--prmtop', required=True, help='AMBER topology file (.prmtop)')
    parser.add_argument('--rst7', required=True, help='AMBER restart file (.rst7)')
    parser.add_argument('--output_prefix', required=True, help='Output prefix for GROMACS files')
    
    args = parser.parse_args()
    
    # Ensure output directory exists
    output_dir = os.path.dirname(args.output_prefix)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    print(f"Converting AMBER files to GROMACS format...")
    print(f"Input topology: {args.prmtop}")
    print(f"Input coordinates: {args.rst7}")
    print(f"Output prefix: {args.output_prefix}")
    
    try:
        # Convert using acpype
        cmd = f"acpype -p {args.prmtop} -x {args.rst7} -o gmx"
        print(f"Running: {cmd}")
        
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"Error running acpype: {result.stderr}")
            sys.exit(1)
        
        # Find the generated files - acpype creates a directory with the files
        gro_file = None
        itp_file = None
        
        print("Looking for generated files...")
        
        # Look for the acpype output directory
        acpype_dir = None
        for item in os.listdir('.'):
            if item.endswith('.amb2gmx'):
                acpype_dir = item
                print(f"Found acpype directory: {item}")
                break
        
        if acpype_dir:
            print(f"Looking in directory: {acpype_dir}")
            for file in os.listdir(acpype_dir):
                print(f"Found file: {file}")
                if file.endswith('_GMX.gro'):
                    gro_file = os.path.join(acpype_dir, file)
                    print(f"Found GRO file: {file}")
                elif file.endswith('_GMX.top'):
                    itp_file = os.path.join(acpype_dir, file)
                    print(f"Found TOP file: {file}")
        
        if gro_file and itp_file:
            # Move files to desired location
            os.rename(gro_file, f"{args.output_prefix}.gro")
            os.rename(itp_file, f"{args.output_prefix}.itp")
            
            print(f"✅ Conversion successful!")
            print(f"📁 GRO file: {args.output_prefix}.gro")
            print(f"📁 ITP file: {args.output_prefix}.itp")
        else:
            print("❌ Generated files not found")
            print(f"GRO file found: {gro_file}")
            print(f"ITP file found: {itp_file}")
            sys.exit(1)
            
    except Exception as e:
        print(f"❌ Error during conversion: {e}")
        sys.exit(1)
 
if __name__ == "__main__":
    main() 