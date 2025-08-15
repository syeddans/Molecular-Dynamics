#!/usr/bin/env python3
"""
Create GROMACS index file with ActiveSite group for SMD simulations.
This creates an .ndx file that can be used with gmx grompp.
"""

import argparse
import os

def read_active_site_residues(residues_file):
    """Read active site residue numbers from file."""
    with open(residues_file, 'r') as f:
        content = f.read().strip()
        if not content:
            raise ValueError("Active site residues file is empty")
        
        residues = [int(r.strip()) for r in content.split(',') if r.strip()]
        return residues

def create_index_file(gro_file, topology_file, active_site_residues, output_dir):
    """Create GROMACS index file with standard groups plus ActiveSite group."""
    import subprocess
    import tempfile
    
    # Read the GRO file to get atom indices for each residue
    atom_indices = []
    
    with open(gro_file, 'r') as f:
        lines = f.readlines()
        
        # Skip header and get number of atoms
        if len(lines) < 3:
            raise ValueError("Invalid GRO file format")
        
        try:
            n_atoms = int(lines[1].strip())
        except ValueError:
            raise ValueError("Could not parse number of atoms from GRO file")
        
        # Parse atoms (lines 2 to n_atoms+1)
        for i in range(2, min(2 + n_atoms, len(lines))):
            line = lines[i]
            if len(line) < 20:  # Basic GRO format check
                continue
                
            # GRO format: resid(5) + resname(5) + atomname(5) + atomid(5) + ...
            try:
                res_id = int(line[0:5].strip())
                atom_id = int(line[15:20].strip())
                
                # Check if this residue is in our active site list
                if res_id in active_site_residues:
                    atom_indices.append(atom_id)
            except (ValueError, IndexError):
                continue
    
    if not atom_indices:
        raise ValueError(f"No atoms found for active site residues {active_site_residues}")
    
    # Sort atom indices
    atom_indices.sort()
    
    # First, generate standard index groups using gmx make_ndx
    index_file = os.path.join(output_dir, "index.ndx")
    
    # Create temporary input file for gmx make_ndx (just quit to generate default groups)
    with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as temp_input:
        temp_input.write("q\n")  # Just quit to generate default groups
        temp_input_path = temp_input.name
    
    try:
        # Generate standard groups using gmx make_ndx
        subprocess.run([
            'gmx', 'make_ndx',
            '-f', gro_file,
            '-o', index_file
        ], input='q\n', text=True, check=True, capture_output=True)
        
        # Now append our ActiveSite group to the generated index file
        with open(index_file, 'a') as f:
            f.write("\n[ ActiveSite ]\n")
            
            # Write atom indices, 15 per line (GROMACS standard)
            for i, atom_idx in enumerate(atom_indices):
                if i > 0 and i % 15 == 0:
                    f.write("\n")
                f.write(f"{atom_idx:>6}")
            f.write("\n")
        
    finally:
        # Clean up temporary file
        os.unlink(temp_input_path)
    
    print(f"Created ActiveSite group with {len(atom_indices)} atoms from {len(active_site_residues)} residues: {active_site_residues}")
    print(f"Index file written to: {index_file}")
    
    return index_file

def main():
    parser = argparse.ArgumentParser(description="Create GROMACS index file with ActiveSite group")
    parser.add_argument("--gro", required=True, help="Input GRO file for atom indices")
    parser.add_argument("--topology", required=True, help="Input topology file")
    parser.add_argument("--residues", required=True, help="Active site residues file")
    parser.add_argument("--output", required=True, help="Output directory for index file")
    
    args = parser.parse_args()
    
    # Read active site residues
    active_site_residues = read_active_site_residues(args.residues)
    
    # Create index file
    create_index_file(args.gro, args.topology, active_site_residues, args.output)

if __name__ == "__main__":
    main()
