#!/usr/bin/env python3
"""
Convert residue numbers from original PDB to GROMACS-processed PDB using offset
"""

import argparse
import os
from Bio.PDB import PDBParser

def map_residue_numbers(original_pdb, gromacs_pdb, active_site_residues_file, output_file):
    """
    Map residue numbers from original PDB to GROMACS-processed PDB using offset
    
    Args:
        original_pdb: Path to original PDB file
        gromacs_pdb: Path to GROMACS-processed PDB file  
        active_site_residues_file: Path to file with active site residue numbers (original numbering)
        output_file: Output file for converted residue numbers
    """
    
    # Parse PDBs
    parser = PDBParser(QUIET=True)
    original_structure = parser.get_structure("original", original_pdb)
    gromacs_structure = parser.get_structure("gromacs", gromacs_pdb)
    
    # Read active site residues from file
    with open(active_site_residues_file, 'r') as f:
        content = f.read().strip()
        if not content:
            print("Warning: Active site residues file is empty")
            return []
        
        original_residues = [int(r.strip()) for r in content.split(',') if r.strip()]
    
    print(f"Original active site residues: {original_residues}")
    
    # Get residue sequences from both PDBs
    original_sequence = []
    gromacs_sequence = []
    
    # Get original PDB sequence (chain A only)
    for chain in original_structure[0]:
        if chain.id == 'A':
            for res in chain:
                if res.id[0] == " ":  # standard residue
                    original_sequence.append((res.get_id()[1], res.resname))
            break
    
    # Get GROMACS PDB sequence
    for chain in gromacs_structure[0]:
        for res in chain:
            if res.id[0] == " ":  # standard residue
                gromacs_sequence.append((res.get_id()[1], res.resname))
        break
    
    print(f"Original PDB has {len(original_sequence)} residues")
    print(f"GROMACS PDB has {len(gromacs_sequence)} residues")
    
    # Find the offset by matching sequences
    print("Finding offset between original and GROMACS PDB...")
    
    # Find a good starting point for matching (first few residues)
    offset = None
    for i in range(min(10, len(original_sequence), len(gromacs_sequence))):
        orig_resnum, orig_name = original_sequence[i]
        gromacs_resnum, gromacs_name = gromacs_sequence[i]
        
        if orig_name == gromacs_name:
            offset = orig_resnum - gromacs_resnum
            print(f"Found offset: {offset} (residue {orig_name} {orig_resnum} -> {gromacs_resnum})")
            break
    
    if offset is None:
        print("Error: Could not determine offset between PDBs")
        return []
    
    # Convert active site residues using offset
    converted_residues = []
    for orig_resnum in original_residues:
        converted_resnum = orig_resnum - offset
        
        # Verify the conversion is correct
        orig_res = None
        gromacs_res = None
        
        # Find original residue
        for resnum, resname in original_sequence:
            if resnum == orig_resnum:
                orig_res = resname
                break
        
        # Find GROMACS residue
        for resnum, resname in gromacs_sequence:
            if resnum == converted_resnum:
                gromacs_res = resname
                break
        
        if orig_res and gromacs_res and orig_res == gromacs_res:
            converted_residues.append(converted_resnum)
            print(f"✓ Mapped {orig_res} {orig_resnum} -> {converted_resnum}")
        else:
            print(f"✗ Could not verify mapping for residue {orig_resnum} -> {converted_resnum}")
            print(f"  Original: {orig_res}, GROMACS: {gromacs_res}")
    
    # Sort and remove duplicates
    converted_residues = sorted(list(set(converted_residues)))
    
    # Write converted residues to file
    with open(output_file, 'w') as f:
        f.write(','.join(str(r) for r in converted_residues))
    
    print(f"Converted active site residues: {converted_residues}")
    print(f"Saved to: {output_file}")
    
    return converted_residues

def main():
    parser = argparse.ArgumentParser(description='Convert residue numbers from original to GROMACS PDB')
    parser.add_argument('--original_pdb', required=True, help='Original PDB file')
    parser.add_argument('--gromacs_pdb', required=True, help='GROMACS-processed PDB file')
    parser.add_argument('--active_site_residues', required=True, help='File with original active site residue numbers')
    parser.add_argument('--output', required=True, help='Output file for converted residue numbers')
    
    args = parser.parse_args()
    
    try:
        converted_residues = map_residue_numbers(
            args.original_pdb,
            args.gromacs_pdb, 
            args.active_site_residues,
            args.output
        )
        
        if converted_residues:
            print(f"✅ Successfully converted {len(converted_residues)} active site residues")
        else:
            print("❌ No residues were converted")
            
    except Exception as e:
        print(f"Error: {e}")
        exit(1)

if __name__ == "__main__":
    main() 