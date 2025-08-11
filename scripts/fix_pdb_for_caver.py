#!/usr/bin/env python3
"""
Fix PDB file for CAVER analysis by adding missing atoms and hydrogens.
"""

import argparse
import sys
from pdbfixer import PDBFixer
from openmm.app import PDBFile

def fix_pdb_for_caver(input_pdb, output_pdb):
    """
    Fix PDB file for CAVER analysis.
    
    Args:
        input_pdb (str): Input PDB file path
        output_pdb (str): Output PDB file path
    """
    try:
        # Create PDBFixer object
        fixer = PDBFixer(filename=input_pdb)
        
        # Find and add missing residues
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        
        # Add missing hydrogens with pH 7.0
        fixer.addMissingHydrogens(7.0)
        
        # Write the fixed structure
        PDBFile.writeFile(fixer.topology, fixer.positions, open(output_pdb, "w"))
        
        print(f"Successfully fixed PDB file: {input_pdb} -> {output_pdb}")
        
    except Exception as e:
        print(f"Error fixing PDB file: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Fix PDB file for CAVER analysis")
    parser.add_argument("--input", required=True, help="Input PDB file")
    parser.add_argument("--output", required=True, help="Output PDB file")
    
    args = parser.parse_args()
    
    fix_pdb_for_caver(args.input, args.output)

if __name__ == "__main__":
    main() 