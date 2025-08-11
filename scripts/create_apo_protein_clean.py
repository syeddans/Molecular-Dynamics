#!/usr/bin/env python3
"""
Clean protein by removing ligand and using PDBFixer for robust cleaning
Find active site residues in GROMACS-processed single-chain structure
"""

import argparse
import os
import pymol
from pymol import cmd
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from Bio.PDB import PDBParser, NeighborSearch

def find_active_site_residues(pdb_file, ligand_resname, cutoff=5.0, output_file=None):
    """
    Find active site residues within cutoff distance of the ligand.
    GROMACS doesn't preserve chain information, so we work without chains.
    
    Args:
        pdb_file (str): Input PDB file path (GROMACS-processed)
        ligand_resname (str): Ligand residue name
        cutoff (float): Cutoff distance in Angstroms (default: 10.0)
        output_file (str): Output file for active site residues
    
    Returns:
        list: List of residue numbers
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdb_file)
    model = structure[0]
    
    protein_atoms = []
    ligand_atoms = []
    
    # Collect all protein and ligand atoms (no chain filtering needed)
    for chain in model:
        print(f"Chain {chain.id}: {len(list(chain))} residues")
        for res in chain:
            print(f"  Residue: {res.resname} {res.get_id()} (resname.strip()='{res.resname.strip()}')")
            if res.resname.strip() == ligand_resname:
                ligand_atoms.extend(res.get_atoms())
                print(f"Found ligand residue: {res.resname} {res.get_id()}")
            elif res.id[0] == " ":  # standard residue
                protein_atoms.extend(res.get_atoms())
            else:
                print(f"  Not protein or ligand: {res.resname} {res.get_id()}")
                if res.resname.strip() == "EST":
                    print(f"    EST found but not matching ligand_resname='{ligand_resname}'")
                    print(f"    Comparison: '{res.resname.strip()}' == '{ligand_resname}' = {res.resname.strip() == ligand_resname}")
    
    print(f"Total protein atoms: {len(protein_atoms)}")
    print(f"Total ligand atoms: {len(ligand_atoms)}")
    print(f"Looking for ligand: '{ligand_resname}'")
    
    if not ligand_atoms:
        print(f"Warning: No ligand atoms found for {ligand_resname}")
        print("Available residues:")
        for chain in model:
            for res in chain:
                print(f"  {res.resname} {res.get_id()}")
        return []
    
    if not protein_atoms:
        print(f"Warning: No protein atoms found")
        return []
    
    # Debug: Show ligand and protein coordinate ranges
    ligand_coords = [atom.coord for atom in ligand_atoms]
    protein_coords = [atom.coord for atom in protein_atoms[:10]]  # First 10 protein atoms
    
    print(f"Ligand coordinate range: X={min(c[0] for c in ligand_coords):.1f}-{max(c[0] for c in ligand_coords):.1f}, "
          f"Y={min(c[1] for c in ligand_coords):.1f}-{max(c[1] for c in ligand_coords):.1f}, "
          f"Z={min(c[2] for c in ligand_coords):.1f}-{max(c[2] for c in ligand_coords):.1f}")
    print(f"Protein coordinate range: X={min(c[0] for c in protein_coords):.1f}-{max(c[0] for c in protein_coords):.1f}, "
          f"Y={min(c[1] for c in protein_coords):.1f}-{max(c[1] for c in protein_coords):.1f}, "
          f"Z={min(c[2] for c in protein_coords):.1f}-{max(c[2] for c in protein_coords):.1f}")
    print(f"Using cutoff distance: {cutoff} Å")
    
    # Build neighbor search from all protein atoms
    ns = NeighborSearch(list(protein_atoms))
    
    # Find all residues with atoms within cutoff of any ligand atom
    active_site_residues = set()
    for atom in ligand_atoms:
        nearby_atoms = ns.search(atom.coord, cutoff)
        for neighbor in nearby_atoms:
            res = neighbor.get_parent()
            if res.id[0] == " ":  # standard residue
                active_site_residues.add(res.get_id()[1])
    
    # Sort residues
    active_site_residues = sorted(active_site_residues)
    
    # Save to file if output_file is provided
    if output_file:
        with open(output_file, 'w') as f:
            if active_site_residues:
                f.write(','.join(str(r) for r in active_site_residues))
                print(f"Active site residues saved to: {output_file}")
                print(f"Residues written: {active_site_residues}")
            else:
                f.write("")
                print(f"Warning: No active site residues found, creating empty file: {output_file}")
    
    print(f"Found {len(active_site_residues)} active site residues: {active_site_residues}")
    return active_site_residues

def main():
    parser = argparse.ArgumentParser(description='Clean protein by removing ligand and keeping only chain A')
    parser.add_argument('--input', required=True, help='Input holo PDB file')
    parser.add_argument('--output', required=True, help='Output apo PDB file')
    parser.add_argument('--ligand_resn', default='EST', help='Ligand residue name to remove')
    parser.add_argument('--chain', default='A', help='Chain to keep (default: A)')
    parser.add_argument('--active_site_output', help='Output file for active site residues')
    parser.add_argument('--cutoff', type=float, default=5.0, help='Cutoff distance for active site (default: 10.0 Å)')
    
    args = parser.parse_args()
    
    # Ensure directory exists
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    
    # Step 1: Load original PDB and remove second chain
    print(f"=== Step 1: Loading original PDB and removing second chain ===")
    temp_single_chain = os.path.join(os.path.dirname(args.output), "temp_single_chain.pdb")
    
    # Initialize PyMOL in headless mode
    pymol.pymol_argv = ['pymol', '-c']
    pymol.finish_launching()
    
    print(f"Loading original PDB: {args.input}")
    cmd.load(args.input, "original_structure")
    
    chains = cmd.get_chains("original_structure")
    print(f"Available chains: {chains}")
    
    other_chains = [c for c in chains if c != args.chain]
    print(f"Removing chains: {other_chains}")
    
    for chain in other_chains:
        print(f"Removing chain {chain} and its ligand {args.ligand_resn}")
        cmd.remove(f"chain {chain}")
    
    # Step 2: Save single-chain PDB with ligand
    print(f"=== Step 2: Saving single-chain structure with ligand ===")
    print(f"Saving single-chain structure with ligand to: {temp_single_chain}")
    cmd.save(temp_single_chain, "original_structure")
    
    # Step 3: Find active site residues from the single-chain PDB
    if args.active_site_output:
        print(f"=== Step 3: Finding active site residues around ligand {args.ligand_resn} ===")
        
        active_site_residues = find_active_site_residues(
            temp_single_chain, 
            args.ligand_resn, 
            args.cutoff,
            args.active_site_output
        )
        
        print(f"Active site detection completed")
    
    # Step 4: Remove the ligand (keep PyMOL open)
    print(f"=== Step 4: Removing ligand {args.ligand_resn} ===")
    cmd.remove(f"resn {args.ligand_resn}")
    print(f"Removed ligand {args.ligand_resn}")
    
    # Step 5: Save ligand-free PDB
    print(f"=== Step 5: Saving ligand-free structure ===")
    temp_no_ligand = os.path.join(os.path.dirname(args.output), "temp_no_ligand.pdb")
    cmd.save(temp_no_ligand, "original_structure")
    print(f"Saved ligand-free structure to: {temp_no_ligand}")
    
    # Close PyMOL
    cmd.delete("all")
    pymol.cmd.quit()
    
    # Step 4: Clean the ligand-free PDB with PDBFixer
    print(f"=== Step 4: Cleaning ligand-free PDB with PDBFixer ===")
    
    fixer = PDBFixer(filename=temp_no_ligand)
    
    # Remove all heterogens (waters, ions, other ligands)
    fixer.removeHeterogens(keepWater=False)
    
    # Fix missing residues and atoms (correct order)
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    
    print(f"PDBFixer fixed {len(fixer.missingAtoms)} missing atoms")
    
    # Save the final cleaned structure
    with open(args.output, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)
    
    # Clean up temp files
    os.remove(temp_single_chain)
    os.remove(temp_no_ligand)
    


if __name__ == "__main__":
    main() 