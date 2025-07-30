#!/usr/bin/env python3
"""
Clean protein by removing ligand and using PDBFixer for robust cleaning
"""

import argparse
import os
import pymol
from pymol import cmd
from pdbfixer import PDBFixer
from openmm.app import PDBFile

def main():
    parser = argparse.ArgumentParser(description='Clean protein by removing ligand')
    parser.add_argument('--input', required=True, help='Input holo PDB file')
    parser.add_argument('--output', required=True, help='Output apo PDB file')
    parser.add_argument('--ligand_resn', default='EST', help='Ligand residue name to remove')
    
    args = parser.parse_args()
    
    # Ensure directory exists
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    # Initialize PyMOL in headless mode
    pymol.pymol_argv = ['pymol', '-c']
    pymol.finish_launching()
    
    # Step 1: Remove specific ligand and keep all protein chains with PyMOL
    print(f"Removing ligand {args.ligand_resn} and keeping all protein chains...")
    cmd.load(args.input, "structure")
    
    # Remove the specific ligand
    cmd.remove(f"resn {args.ligand_resn}")
    
    # Keep all protein chains - just show what chains are available
    chains = cmd.get_chains("structure")
    print(f"Available chains: {chains}")
    print(f"Keeping all protein chains: {', '.join(chains)}")
    
    # Save multi-chain protein without ligand
    temp_file = os.path.join(os.path.dirname(args.output), "temp_no_ligand.pdb")
    cmd.save(temp_file, "structure")
    cmd.delete("all")
    pymol.cmd.quit()
    
    # Step 2: Use PDBFixer to clean up the structure (automatic and robust)
    print("Cleaning structure with PDBFixer...")
    fixer = PDBFixer(filename=temp_file)
    
    # Remove all heterogens (waters, ions, other ligands)
    fixer.removeHeterogens(keepWater=False)
    
    # Fix missing residues and atoms (correct order)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    
    print(f"PDBFixer fixed {len(fixer.missingAtoms)} missing atoms")
    
    # Save the cleaned structure
    with open(args.output, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    
    # Clean up temp file
    os.remove(temp_file)
    print(f"✅ Clean multi-chain apo protein saved: {args.output}")

if __name__ == "__main__":
    main() 