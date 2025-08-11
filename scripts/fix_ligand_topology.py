#!/usr/bin/env python3
"""
Script to fix ligand topology by removing defaults and atomtypes sections.
"""

import argparse
import os

def fix_ligand_topology(ligand_itp, output_itp, add_posres=False):
    """
    Fix ligand topology by removing defaults and atomtypes sections.
    
    Args:
        ligand_itp: Path to original ligand .itp file
        output_itp: Path to output fixed ligand .itp file
        add_posres: Whether to add position restraints to the ligand
    """
    
    # Read ligand topology
    with open(ligand_itp, 'r') as f:
        ligand_lines = f.readlines()
    
    output_lines = []
    in_defaults = False
    in_atomtypes = False
    in_molecules = False
    in_system = False

    for line in ligand_lines:
        # Skip defaults section (conflicts with main force field)
        if line.strip() == '[ defaults ]':
            in_defaults = True
            continue
        elif line.strip().startswith('[') and in_defaults:
            in_defaults = False
        elif in_defaults:
            continue

        # Skip atomtypes section (will be included in main topology)
        if line.strip() == '[ atomtypes ]':
            in_atomtypes = True
            continue
        elif line.strip().startswith('[') and in_atomtypes:
            in_atomtypes = False
        elif in_atomtypes:
            continue

        # Skip molecules section
        if line.strip() == '[ molecules ]':
            in_molecules = True
            continue
        elif line.strip().startswith('[') and in_molecules:
            in_molecules = False
        elif in_molecules:
            continue

        # Skip system section
        if line.strip() == '[ system ]':
            in_system = True
            continue
        elif line.strip().startswith('[') and in_system:
            in_system = False
        elif in_system:
            continue

        # Include all other lines
        output_lines.append(line)
    
    # Add position restraints if requested
    if add_posres:
        # Count the number of atoms in the ligand by looking at the [ atoms ] section
        atom_count = 0
        in_atoms_section = False
        for line in output_lines:
            if line.strip() == '[ atoms ]':
                in_atoms_section = True
                continue
            elif line.strip().startswith('[') and in_atoms_section:
                # End of atoms section
                break
            elif in_atoms_section and line.strip() and not line.strip().startswith(';'):
                # This is an atom line
                atom_count += 1
        
        print(f"Found {atom_count} atoms in ligand")
        
        # Add position restraints section
        output_lines.append(f"\n; Position restraints for ligand\n")
        output_lines.append(f"#ifdef POSRES\n")
        output_lines.append(f"[ position_restraints ]\n")
        output_lines.append(f";  i funct       fcx        fcy        fcz\n")
        
        for i in range(1, atom_count + 1):
            output_lines.append(f"{i:6d}    1   10000   10000   10000\n")
        
        output_lines.append(f"#endif\n")
    
    # Write the fixed ligand topology
    with open(output_itp, 'w') as f:
        f.writelines(output_lines)
    
    print(f"Created fixed ligand topology: {output_itp}")
    print(f"Removed defaults section to avoid conflicts with main force field")
    print(f"Removed atomtypes section (will be included in main topology)")
    if add_posres:
        print(f"Added position restraints for {atom_count} atoms")

def main():
    parser = argparse.ArgumentParser(description='Fix ligand topology by removing defaults and atomtypes sections')
    parser.add_argument('--input', required=True, help='Input ligand .itp file')
    parser.add_argument('--output', required=True, help='Output fixed ligand .itp file')
    parser.add_argument('--add-posres', action='store_true', help='Add position restraints to ligand')
    
    args = parser.parse_args()
    
    fix_ligand_topology(args.input, args.output, args.add_posres)

if __name__ == "__main__":
    main()