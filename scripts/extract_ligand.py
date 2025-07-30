#!/usr/bin/env python3
"""
Extract clean protein from holo structure by removing ligands, waters, ions, etc.
Keeps only protein chains for MD simulation preparation
"""

import sys
import argparse
import pymol
from pymol import cmd

def extract_clean_protein(input_pdb, output_pdb, ligand_resn="EST", keep_chain=None):
    """Extract clean protein by removing everything except protein chains"""
    
    # Initialize PyMOL in headless mode
    pymol.pymol_argv = ['pymol', '-c']
    pymol.finish_launching()
    
    try:
        # Load the structure
        print(f"Loading structure: {input_pdb}")
        cmd.load(input_pdb, "structure")
        
        # Get information about the structure
        chains = cmd.get_chains("structure")
        print(f"Available chains: {chains}")
        
        # Remove the specified ligand
        if ligand_resn:
            cmd.remove(f"resn {ligand_resn}")
            print(f"Removed ligand: {ligand_resn}")
        
        # Remove all heterogens (waters, ions, other ligands, etc.)
        # Keep only standard amino acids
        cmd.remove("not polymer.protein")
        print("Removed all non-protein molecules (waters, ions, other ligands)")
        
        # If specific chain requested, keep only that chain
        if keep_chain:
            if keep_chain in chains:
                cmd.remove(f"not chain {keep_chain}")
                print(f"Kept only chain {keep_chain}")
            else:
                print(f"Warning: Chain {keep_chain} not found, keeping all protein chains")
        
        # Remove any remaining non-standard residues
        cmd.remove("not (resn ALA+ARG+ASN+ASP+CYS+GLN+GLU+GLY+HIS+ILE+LEU+LYS+MET+PHE+PRO+SER+THR+TRP+TYR+VAL)")
        
        # Get final chain information
        final_chains = cmd.get_chains("structure")
        atom_count = cmd.count_atoms("structure")
        
        print(f"Final structure:")
        print(f"  Chains: {final_chains}")
        print(f"  Atoms: {atom_count}")
        
        # Save the clean protein
        cmd.save(output_pdb, "structure")
        print(f"Clean protein saved to: {output_pdb}")
        
    except Exception as e:
        print(f"Error processing structure: {e}")
        sys.exit(1)
    
    finally:
        # Clean up
        cmd.delete("all")
        cmd.quit()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract clean protein from PDB file (remove ligands, waters, ions)")
    parser.add_argument("--input", required=True, help="Input PDB file")
    parser.add_argument("--output", required=True, help="Output clean protein PDB file")
    parser.add_argument("--ligand_resn", default="EST", help="Specific ligand residue name to remove (default: EST)")
    parser.add_argument("--keep_chain", help="Keep only specified chain (e.g., A, B, C)")
    
    args = parser.parse_args()
    
    extract_clean_protein(args.input, args.output, args.ligand_resn, args.keep_chain) 