#!/usr/bin/env python3
"""
Translate a pose-oriented ligand to the tunnel entrance position.

The pose-oriented ligand has the correct orientation from AutoDock Vina
but is positioned at the binding site. We need to translate it to the
tunnel entrance position (2nm outside) while preserving the orientation.
"""

import argparse
import numpy as np
import os

def read_gro_coordinates(gro_file):
    """Read coordinates from a GRO file."""
    atoms = []
    with open(gro_file, 'r') as f:
        lines = f.readlines()
    
    # Parse header
    title = lines[0].strip()
    natoms = int(lines[1].strip())
    
    # Parse atoms
    for i in range(2, 2 + natoms):
        line = lines[i]
        # Parse GRO format: residue_id residue_name atom_name atom_id x y z
        try:
            residue_id = int(line[0:5])
            residue_name = line[5:10].strip()
            atom_name = line[10:15].strip()
            atom_id = int(line[15:20])
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
            
            atoms.append({
                'residue_id': residue_id,
                'residue_name': residue_name,
                'atom_name': atom_name,
                'atom_id': atom_id,
                'coords': np.array([x, y, z]),
                'line_format': line
            })
        except (ValueError, IndexError):
            print(f"Warning: Could not parse line {i}: {line.strip()}")
            continue
    
    # Box line
    box_line = lines[2 + natoms].strip() if len(lines) > 2 + natoms else "0.0 0.0 0.0"
    
    return title, atoms, box_line

def read_tunnel_position(position_file):
    """Read tunnel position from the placement file."""
    with open(position_file, 'r') as f:
        line = f.read().strip()
    
    coords = [float(x) for x in line.split()]
    if len(coords) != 3:
        raise ValueError(f"Expected 3 coordinates, got {len(coords)}: {coords}")
    
    return np.array(coords)

def compute_center_of_mass(atoms):
    """Compute center of mass of ligand atoms."""
    coords = np.array([atom['coords'] for atom in atoms])
    return np.mean(coords, axis=0)

def write_gro_file(title, atoms, box_line, output_file):
    """Write atoms to a GRO file."""
    with open(output_file, 'w') as f:
        f.write(f"{title}\n")
        f.write(f"{len(atoms)}\n")
        
        for atom in atoms:
            # Format: residue_id(5) residue_name(5) atom_name(5) atom_id(5) x(8.3) y(8.3) z(8.3)
            f.write(f"{atom['residue_id']:5d}{atom['residue_name']:>5s}{atom['atom_name']:>5s}{atom['atom_id']:5d}"
                   f"{atom['coords'][0]:8.3f}{atom['coords'][1]:8.3f}{atom['coords'][2]:8.3f}\n")
        
        f.write(f"{box_line}\n")

def main():
    parser = argparse.ArgumentParser(description="Translate pose-oriented ligand to tunnel entrance position")
    parser.add_argument("--pose_gro", required=True, help="Pose-oriented ligand GRO file")
    parser.add_argument("--tunnel_position", required=True, help="File containing tunnel entrance position (x y z in nm)")
    parser.add_argument("--output", required=True, help="Output GRO file with ligand at tunnel position")
    
    args = parser.parse_args()
    
    # Read pose-oriented ligand
    title, atoms, box_line = read_gro_coordinates(args.pose_gro)
    
    if not atoms:
        raise ValueError(f"No atoms found in {args.pose_gro}")
    
    # Read tunnel position
    tunnel_pos = read_tunnel_position(args.tunnel_position)
    
    # Compute current center of mass
    current_com = compute_center_of_mass(atoms)
    
    # Compute translation vector
    translation = tunnel_pos - current_com
    
    print(f"Current ligand COM: [{current_com[0]:.3f}, {current_com[1]:.3f}, {current_com[2]:.3f}] nm")
    print(f"Target tunnel position: [{tunnel_pos[0]:.3f}, {tunnel_pos[1]:.3f}, {tunnel_pos[2]:.3f}] nm")
    print(f"Translation vector: [{translation[0]:.3f}, {translation[1]:.3f}, {translation[2]:.3f}] nm")
    
    # Apply translation to all atoms
    for atom in atoms:
        atom['coords'] += translation
    
    # Verify final position
    final_com = compute_center_of_mass(atoms)
    distance_to_target = np.linalg.norm(final_com - tunnel_pos)
    
    print(f"Final ligand COM: [{final_com[0]:.3f}, {final_com[1]:.3f}, {final_com[2]:.3f}] nm")
    print(f"Distance to target: {distance_to_target:.6f} nm")
    
    if distance_to_target > 0.001:
        print("Warning: Final position does not match target position!")
    
    # Write output
    output_dir = os.path.dirname(args.output)
    if output_dir:  # Only create directory if path contains a directory
        os.makedirs(output_dir, exist_ok=True)
    write_gro_file(title, atoms, box_line, args.output)
    
    print(f"Translated ligand written to: {args.output}")

if __name__ == "__main__":
    main() 