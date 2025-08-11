#!/usr/bin/env python3
"""
Directly combine protein and ligand coordinates without using gmx insert-molecules.

This preserves the exact ligand positioning at the tunnel entrance.
"""

import argparse
import os

def read_gro_file(gro_file):
    """Read a GRO file and return title, atoms, and box line."""
    with open(gro_file, 'r') as f:
        lines = f.readlines()
    
    title = lines[0].strip()
    natoms = int(lines[1].strip())
    
    atoms = []
    for i in range(2, 2 + natoms):
        atoms.append(lines[i])
    
    box_line = lines[2 + natoms].strip() if len(lines) > 2 + natoms else "0.0 0.0 0.0"
    
    return title, atoms, box_line

def write_gro_file(title, atoms, box_line, output_file):
    """Write a combined GRO file."""
    with open(output_file, 'w') as f:
        f.write(f"{title}\n")
        f.write(f"{len(atoms)}\n")
        
        for atom in atoms:
            f.write(atom)
        
        f.write(f"{box_line}\n")

def renumber_atoms(atoms, start_id=1):
    """Renumber atom IDs sequentially."""
    renumbered = []
    atom_id = start_id
    
    for line in atoms:
        if len(line.strip()) == 0:
            continue
        
        # Parse and rebuild the line with new atom ID
        residue_id = line[0:5]
        residue_name = line[5:10]
        atom_name = line[10:15]
        # Skip old atom ID (15:20)
        coords = line[20:44]  # x, y, z coordinates
        velocities = line[44:] if len(line) > 44 else "\n"
        
        new_line = f"{residue_id}{residue_name}{atom_name}{atom_id:5d}{coords}{velocities}"
        renumbered.append(new_line)
        atom_id += 1
    
    return renumbered

def main():
    parser = argparse.ArgumentParser(description="Combine protein and ligand GRO files directly")
    parser.add_argument("--protein", required=True, help="Protein GRO file")
    parser.add_argument("--ligand", required=True, help="Ligand GRO file (positioned at tunnel)")
    parser.add_argument("--output", required=True, help="Output combined GRO file")
    
    args = parser.parse_args()
    
    # Read protein structure
    protein_title, protein_atoms, protein_box = read_gro_file(args.protein)
    
    # Read ligand structure
    ligand_title, ligand_atoms, ligand_box = read_gro_file(args.ligand)
    
    print(f"Protein atoms: {len(protein_atoms)}")
    print(f"Ligand atoms: {len(ligand_atoms)}")
    
    # Combine atoms
    combined_atoms = protein_atoms + ligand_atoms
    
    # Renumber all atoms sequentially
    combined_atoms = renumber_atoms(combined_atoms)
    
    # Use protein box (larger)
    combined_title = f"Protein-ligand complex (tunnel positioned)"
    
    print(f"Combined atoms: {len(combined_atoms)}")
    
    # Write combined structure
    output_dir = os.path.dirname(args.output)
    if output_dir:  # Only create directory if path contains a directory
        os.makedirs(output_dir, exist_ok=True)
    write_gro_file(combined_title, combined_atoms, protein_box, args.output)
    
    print(f"Combined structure written to: {args.output}")

if __name__ == "__main__":
    main() 