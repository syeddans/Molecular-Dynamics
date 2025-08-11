#!/usr/bin/env python3
"""
Script to create complex topology by merging protein and ligand topologies.
Based on GROMACS tutorial guidelines for including ligand parameters.
"""

import argparse
import os
import re

def create_complex_topology(protein_top, ligand_itp, output_top):
    """
    Create complex topology by merging protein and ligand topologies.
    
    Args:
        protein_top: Path to protein topology file (clean_apo.top)
        ligand_itp: Path to ligand .itp file (cleaned, without defaults/atomtypes)
        output_top: Path to output complex topology file
    """
    
    # Read the clean_apo.top file
    with open(protein_top, 'r') as f:
        protein_lines = f.readlines()
    
    # Read ligand topology to extract atomtypes
    ligand_atomtypes = []
    with open(ligand_itp, 'r') as f:
        ligand_lines = f.readlines()
    
    # Extract atomtypes from the ligand file
    in_atomtypes = False
    for line in ligand_lines:
        if line.strip() == '[ atomtypes ]':
            in_atomtypes = True
            continue
        elif line.strip().startswith('[') and in_atomtypes:
            # End of atomtypes section
            break
        elif in_atomtypes and line.strip() and not line.strip().startswith(';'):
            ligand_atomtypes.append(line)
    
    # Find the ligand molecule name from the ligand .itp file
    ligand_mol_name = None
    for idx, line in enumerate(ligand_lines):
        if line.strip().startswith('[ moleculetype ]'):
            # Next non-comment, non-empty line should contain the molecule name
            for j in range(idx + 1, min(idx + 6, len(ligand_lines))):
                next_line = ligand_lines[j]
                if next_line.strip() and not next_line.strip().startswith(';'):
                    parts = next_line.split()
                    if len(parts) >= 1:
                        ligand_mol_name = parts[0]
                        break
            break
    
    if not ligand_mol_name:
        raise ValueError("Could not find molecule name in ligand .itp file")
    
    print(f"Found ligand molecule name: {ligand_mol_name}")
    
    # Create output topology by copying the protein topology
    output_lines = []
    
    # Copy protein topology and add atomtypes after forcefield include
    for line in protein_lines:
        output_lines.append(line)
        if 'amber99sb-ildn.ff/forcefield.itp' in line:
            # Add ligand atomtypes immediately after forcefield include
            if ligand_atomtypes:
                output_lines.append(f"\n; Include ligand atomtypes\n")
                output_lines.append(f"[ atomtypes ]\n")
                for atomtype_line in ligand_atomtypes:
                    output_lines.append(atomtype_line)
                output_lines.append(f"\n")
    
    # Always include ligand topology before the [ system ] section, if not already present
    include_line = '#include "clean_ligand.itp"\n'
    has_include_already = any('clean_ligand.itp' in l for l in output_lines)
    if not has_include_already:
        # Find [ system ] section to insert before it
        insert_idx = None
        for i, line in enumerate(output_lines):
            if line.strip().lower().startswith('[ system ]'):
                insert_idx = i
                break
        if insert_idx is None:
            # Fallback: append near the end
            insert_idx = len(output_lines)
        # Insert a small header then the include
        output_lines.insert(insert_idx, "\n")
        output_lines.insert(insert_idx, "; Include ligand topology\n")
        output_lines.insert(insert_idx, include_line)
    
    # Add ligand molecule to the molecules section
    in_molecules = False
    molecules_added = False
    
    for i, line in enumerate(output_lines):
        if line.strip() == '[ molecules ]':
            in_molecules = True
        elif in_molecules and line.strip().startswith('['):
            # End of molecules section, add ligand before this
            output_lines.insert(i, f"{ligand_mol_name}                 1\n")
            molecules_added = True
            break
        elif in_molecules and line.strip() and not line.strip().startswith(';') and not molecules_added:
            # This is a molecule entry, add ligand after the last one if next line starts a new section
            if i + 1 < len(output_lines) and output_lines[i + 1].strip().startswith('['):
                output_lines.insert(i + 1, f"{ligand_mol_name}                 1\n")
                molecules_added = True
                break
            elif i + 1 >= len(output_lines):
                # End of file, add ligand
                output_lines.append(f"{ligand_mol_name}                 1\n")
                molecules_added = True
                break
    
    # If we never found a molecules section, append one
    if not molecules_added and not any(l.strip() == '[ molecules ]' for l in output_lines):
        output_lines.append("\n[ molecules ]\n")
        output_lines.append("; Compound        #mols\n")
        output_lines.append(f"{ligand_mol_name}                 1\n")
    
    # Write the complex topology
    with open(output_top, 'w') as f:
        f.writelines(output_lines)
    
    print(f"Created complex topology: {output_top}")
    print(f"Added ligand molecule: {ligand_mol_name}")


def main():
    parser = argparse.ArgumentParser(description='Create complex or apo topology files')
    parser.add_argument('--protein-top', required=True, help='Protein topology file')
    parser.add_argument('--ligand-itp', help='Ligand .itp file (cleaned)')
    parser.add_argument('--output', required=True, help='Output topology file')
    
    args = parser.parse_args()
    create_complex_topology(args.protein_top, args.ligand_itp, args.output)

if __name__ == "__main__":
    main() 