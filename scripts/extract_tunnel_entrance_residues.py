#!/usr/bin/env python3
"""
Extract tunnel entrance residues from MOLE2 XML output.
This script identifies residues at or near the tunnel entrance coordinates
from the best tunnel selected by MOLE2.
"""

import argparse
import json
import xml.etree.ElementTree as ET
import numpy as np
from typing import Dict, List, Tuple


def parse_mole2_tunnels(xml_path: str) -> List[Dict]:
    """
    Parse MOLE2 tunnels.xml into a list of tunnels with entrance coordinates.
    """
    tree = ET.parse(xml_path)
    root = tree.getroot()
    
    tunnels = []
    
    for tunnel in root.iter():
        if tunnel.tag.lower() != "tunnel":
            continue
            
        tinfo = {}
        
        # Get tunnel length
        for key in ["Length", "length", "Len", "len", "Distance", "distance"]:
            if key in tunnel.attrib:
                try:
                    tinfo["length_A"] = float(tunnel.attrib[key])
                    break
                except ValueError:
                    pass
        
        # Get bottleneck radius
        for key in ["BottleneckRadius", "bottleneckRadius", "bottleneck_radius", 
                   "MinRadius", "minRadius", "min_radius"]:
            if key in tunnel.attrib:
                try:
                    tinfo["bottleneck_radius_A"] = float(tunnel.attrib[key])
                    break
                except ValueError:
                    pass
        
        # Collect node coordinates
        node_coords = []
        for child in tunnel.iter():
            if child.tag.lower() == "node":
                try:
                    x = float(child.attrib.get("x", child.attrib.get("X")))
                    y = float(child.attrib.get("y", child.attrib.get("Y")))
                    z = float(child.attrib.get("z", child.attrib.get("Z")))
                    node_coords.append((x, y, z))
                except Exception:
                    pass
        
        # Entrance is the first node
        if node_coords:
            tinfo["entrance_A"] = list(node_coords[0])
            tinfo["nodes"] = node_coords
            tunnels.append(tinfo)
    
    return tunnels


def select_best_tunnel(tunnels: List[Dict], ligand_radius_A: float = 5.0) -> Dict:
    """
    Select the best tunnel based on bottleneck radius and length.
    Prefer tunnels that can accommodate the ligand and are shortest.
    """
    if not tunnels:
        raise ValueError("No tunnels found in MOLE2 output")
    
    # Filter tunnels that can accommodate the ligand
    suitable_tunnels = []
    for tunnel in tunnels:
        bottleneck = tunnel.get("bottleneck_radius_A", 0)
        if bottleneck >= ligand_radius_A:
            suitable_tunnels.append(tunnel)
    
    # If no suitable tunnels, use the one with largest bottleneck
    if not suitable_tunnels:
        suitable_tunnels = sorted(tunnels, 
                                key=lambda t: t.get("bottleneck_radius_A", 0), 
                                reverse=True)[:1]
    
    # Among suitable tunnels, select the shortest
    best_tunnel = min(suitable_tunnels, 
                     key=lambda t: t.get("length_A", float('inf')))
    
    return best_tunnel


def read_pdb_residues(pdb_path: str) -> List[Tuple[int, str, np.ndarray]]:
    """
    Read PDB file and extract residue information.
    Returns list of (residue_number, chain, ca_coordinates) tuples.
    """
    residues = []
    
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                try:
                    residue_num = int(line[22:26].strip())
                    chain = line[21]
                    x = float(line[30:38])
                    y = float(line[38:46]) 
                    z = float(line[46:54])
                    coords = np.array([x, y, z])
                    residues.append((residue_num, chain, coords))
                except (ValueError, IndexError):
                    continue
    
    return residues


def find_entrance_residues(entrance_coords_A: np.ndarray, 
                          pdb_residues: List[Tuple[int, str, np.ndarray]], 
                          cutoff_A: float = 15.0) -> List[int]:
    """
    Find residues within cutoff distance of tunnel entrance.
    """
    entrance_residues = []
    
    for residue_num, chain, ca_coords in pdb_residues:
        # Calculate distance from CA to tunnel entrance
        distance = np.linalg.norm(ca_coords - entrance_coords_A)
        
        if distance <= cutoff_A:
            entrance_residues.append(residue_num)
    
    return sorted(entrance_residues)


def main():
    parser = argparse.ArgumentParser(
        description="Extract tunnel entrance residues from MOLE2 output"
    )
    parser.add_argument("--tunnels_xml", required=True, 
                       help="Path to MOLE2 tunnels.xml")
    parser.add_argument("--pdb", required=True, 
                       help="PDB file for residue coordinate mapping")
    parser.add_argument("--output", required=True, 
                       help="Output file for tunnel entrance residues")
    parser.add_argument("--cutoff", type=float, default=15.0,
                       help="Distance cutoff in Angstroms for entrance residues (default: 15.0)")
    parser.add_argument("--ligand_radius", type=float, default=5.0,
                       help="Ligand radius in Angstroms for tunnel selection (default: 5.0)")
    parser.add_argument("--json_out", help="Optional JSON output with tunnel info")
    
    args = parser.parse_args()
    
    # Parse MOLE2 tunnels
    tunnels = parse_mole2_tunnels(args.tunnels_xml)
    print(f"Found {len(tunnels)} tunnels from MOLE2")
    
    # Select best tunnel
    best_tunnel = select_best_tunnel(tunnels, args.ligand_radius)
    length_val = best_tunnel.get('length_A', 'N/A')
    bottleneck_val = best_tunnel.get('bottleneck_radius_A', 'N/A')
    
    if length_val != 'N/A':
        length_str = f"{length_val:.1f}Å"
    else:
        length_str = "N/A"
        
    if bottleneck_val != 'N/A':
        bottleneck_str = f"{bottleneck_val:.1f}Å"
    else:
        bottleneck_str = "N/A"
    
    print(f"Selected tunnel: length={length_str}, bottleneck={bottleneck_str}")
    
    # Get entrance coordinates
    if "entrance_A" not in best_tunnel:
        raise ValueError("Selected tunnel has no entrance coordinates")
    
    entrance_coords = np.array(best_tunnel["entrance_A"])
    print(f"Tunnel entrance: [{entrance_coords[0]:.1f}, {entrance_coords[1]:.1f}, {entrance_coords[2]:.1f}]Å")
    
    # Read PDB residues
    pdb_residues = read_pdb_residues(args.pdb)
    print(f"Read {len(pdb_residues)} residues from PDB")
    
    # Find entrance residues
    entrance_residues = find_entrance_residues(entrance_coords, pdb_residues, args.cutoff)
    print(f"Found {len(entrance_residues)} residues within {args.cutoff}Å of tunnel entrance")
    
    # Write output
    with open(args.output, 'w') as f:
        f.write(','.join(map(str, entrance_residues)))
    
    print(f"Tunnel entrance residues: {entrance_residues}")
    print(f"Written to: {args.output}")
    
    # Optional JSON output
    if args.json_out:
        tunnel_info = {
            "best_tunnel": best_tunnel,
            "entrance_coordinates_A": entrance_coords.tolist(),
            "entrance_residues": entrance_residues,
            "cutoff_A": args.cutoff,
            "total_tunnels": len(tunnels)
        }
        
        with open(args.json_out, 'w') as f:
            json.dump(tunnel_info, f, indent=2)
        
        print(f"Tunnel info written to: {args.json_out}")


if __name__ == "__main__":
    main()
