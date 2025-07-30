#!/usr/bin/env python3
"""
Handle position restraint files for single and multi-chain proteins
"""

import argparse
import os
import glob

def main():
    parser = argparse.ArgumentParser(description='Handle position restraint files')
    parser.add_argument('--protein_dir', required=True, help='Protein directory containing posre files')
    
    args = parser.parse_args()
    
    # Handle position restraint files (single vs multi-chain)
    posre_files = glob.glob(os.path.join(args.protein_dir, "posre*.itp"))
    
    if len(posre_files) == 0:
        print("Warning: No position restraint files found!")
    elif len(posre_files) == 1:
        # Single chain - file is already named correctly or rename it
        existing_file = posre_files[0]
        target_file = os.path.join(args.protein_dir, "posre.itp")
        if existing_file != target_file:
            os.rename(existing_file, target_file)
        print("✅ Single-chain protein processed")
    else:
        # Multi-chain - combine all posre files (GROMACS already has correct atom numbers)
        print(f"✅ Multi-chain protein: found {len(posre_files)} chains")
        combined_file = os.path.join(args.protein_dir, "posre.itp")
        
        with open(combined_file, 'w') as combined:
            combined.write("; Combined position restraints for multi-chain protein\\n")
            combined.write("[ position_restraints ]\\n")
            combined.write("; atom  type      fx      fy      fz\\n")
            
            for i, posre_file in enumerate(sorted(posre_files)):
                print(f"  Adding chain {i+1}: {os.path.basename(posre_file)}")
                with open(posre_file, 'r') as f:
                    lines = f.readlines()
                    
                # Extract only the position restraint entries (GROMACS atom numbers are already correct)
                in_posre_section = False
                for line in lines:
                    line_stripped = line.strip()
                    if line_stripped.startswith("[ position_restraints ]"):
                        in_posre_section = True
                        continue
                    elif line_stripped.startswith("[") and in_posre_section:
                        in_posre_section = False
                    elif in_posre_section and line_stripped and not line_stripped.startswith(";"):
                        # Copy the line as-is (atom numbers are already globally correct)
                        combined.write(line)
        
        print(f"✅ Combined {len(posre_files)} position restraint files into posre.itp")

if __name__ == "__main__":
    main() 