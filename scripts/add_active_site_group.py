#!/usr/bin/env python3
"""
Add ActiveSite group definition to GROMACS topology file.
This allows using ActiveSite as a pull group in SMD simulations.
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

def add_active_site_group(topology_file, active_site_residues, output_file):
    """Add ActiveSite group definition to topology file."""
    
    with open(topology_file, 'r') as f:
        lines = f.readlines()
    
    # Find the end of the file or existing [ system ] section
    insert_index = len(lines)
    for i, line in enumerate(lines):
        if line.strip().startswith('[ system ]'):
            insert_index = i
            break
    
    # Create ActiveSite group definition
    group_lines = [
        "\n",
        "[ ActiveSite ]\n",
        "; Active site residues for SMD pulling\n"
    ]
    
    # Add residue ranges for active site atoms
    # Format: residue_start-residue_end for each active site residue
    for res_id in active_site_residues:
        group_lines.append(f"r_{res_id}\n")
    
    group_lines.append("\n")
    
    # Insert group definition before system section
    new_lines = lines[:insert_index] + group_lines + lines[insert_index:]
    
    # Write updated topology
    with open(output_file, 'w') as f:
        f.writelines(new_lines)
    
    print(f"Added ActiveSite group with {len(active_site_residues)} residues: {active_site_residues}")
    print(f"Updated topology written to: {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Add ActiveSite group to GROMACS topology")
    parser.add_argument("--topology", required=True, help="Input topology file (.top)")
    parser.add_argument("--residues", required=True, help="Active site residues file")
    parser.add_argument("--output", required=True, help="Output topology file")
    
    args = parser.parse_args()
    
    # Read active site residues
    active_site_residues = read_active_site_residues(args.residues)
    
    # Add group definition
    add_active_site_group(args.topology, active_site_residues, args.output)

if __name__ == "__main__":
    main()
