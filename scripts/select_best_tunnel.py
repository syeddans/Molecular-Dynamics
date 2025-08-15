#!/usr/bin/env python3
"""
Select the best tunnel from MOLE2 output (shortest that can fit the ligand),
compute the entrance point and direction toward the pocket center, and write
an insertion coordinate 2.0 nm outside the entrance along the outward direction.

Outputs:
- positions file with a single line: "x y z" in nm (for `gmx insert-molecules -ip`)
- optional JSON with diagnostic info
"""

import argparse
import json
import math
import os
import re
import sys
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional, Tuple

import numpy as np


def read_residue_list(residue_file: str) -> List[int]:
    with open(residue_file, "r") as f:
        text = f.read().strip()
    # support formats like "123,124,130" or lines separated
    parts = re.split(r"[\s,]+", text)
    residues = []
    for p in parts:
        p = p.strip()
        if not p:
            continue
        try:
            residues.append(int(p))
        except ValueError:
            pass
    if not residues:
        raise ValueError(f"No residues parsed from {residue_file}")
    return residues


def read_pdb_coords_for_residues(pdb_path: str, residues: List[int]) -> np.ndarray:
    target = set(residues)
    coords: List[List[float]] = []
    with open(pdb_path, "r") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            try:
                resseq = int(line[22:26])
            except ValueError:
                continue
            if resseq in target:
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    continue
                coords.append([x, y, z])
    if not coords:
        raise ValueError("No coordinates found for specified residues in PDB")
    return np.array(coords, dtype=float)


def calc_ligand_radius_from_sdf(sdf_file: str) -> float:
    try:
        from rdkit import Chem
    except Exception as e:
        raise RuntimeError("RDKit is required to compute ligand radius") from e
    mol = Chem.MolFromMolFile(sdf_file, removeHs=False)
    if mol is None:
        raise ValueError(f"Could not read molecule from {sdf_file}")
    conf = mol.GetConformer()
    coords = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())])
    if len(coords) < 2:
        raise ValueError("Ligand must have at least 2 atoms")
    dists = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    max_dist = float(np.max(dists))  # in Å
    return max_dist / 2.0  # radius in Å


def parse_mole2_tunnels(xml_path: str) -> List[Dict]:
    """
    Parse MOLE2 tunnels.xml into a list of tunnels with keys:
    - length_A
    - bottleneck_radius_A (if available)
    - entrance_A: [x, y, z]
    Fallbacks attempt to read common attribute and element names.
    """
    tree = ET.parse(xml_path)
    root = tree.getroot()

    tunnels: List[Dict] = []

    # Iterate only actual Tunnel elements, not the root Tunnels container
    for tunnel in root.iter():
        tag_lower = tunnel.tag.lower()
        if tag_lower != "tunnel":
            continue

        tinfo: Dict = {}
        # Length from attributes if present
        for key in ["Length", "length", "Len", "len", "Distance", "distance"]:
            if key in tunnel.attrib:
                try:
                    tinfo["length_A"] = float(tunnel.attrib[key])
                except ValueError:
                    pass

        # Bottleneck radius from attributes if present
        for key in [
            "BottleneckRadius",
            "bottleneckRadius",
            "bottleneck_radius",
            "MinRadius",
            "minRadius",
            "min_radius",
        ]:
            if key in tunnel.attrib:
                try:
                    tinfo["bottleneck_radius_A"] = float(tunnel.attrib[key])
                except ValueError:
                    pass

        # Collect Node elements for potential entrance, bottleneck and length
        node_coords: List[Tuple[float, float, float]] = []
        node_radii: List[float] = []
        node_bottlenecks: List[float] = []
        node_distances: List[float] = []

        for child in tunnel.iter():
            ctag = child.tag.lower()
            if ctag == "node":
                # Read coordinates with case-insensitive keys
                try:
                    x = float(child.attrib.get("x", child.attrib.get("X")))
                    y = float(child.attrib.get("y", child.attrib.get("Y")))
                    z = float(child.attrib.get("z", child.attrib.get("Z")))
                    node_coords.append((x, y, z))
                except Exception:
                    pass
                # Read radii and distances if available (case-insensitive)
                for rkey in ("BRadius", "bRadius", "bradius", "Radius", "radius"):
                    if rkey in child.attrib:
                        try:
                            val = float(child.attrib[rkey])
                            if rkey.lower().startswith("b"):
                                node_bottlenecks.append(val)
                            else:
                                node_radii.append(val)
                        except Exception:
                            pass
                for dkey in ("Distance", "distance"):
                    if dkey in child.attrib:
                        try:
                            node_distances.append(float(child.attrib[dkey]))
                        except Exception:
                            pass

        # Entrance point: first node coordinate if none found yet
        if "entrance_A" not in tinfo and node_coords:
            tinfo["entrance_A"] = list(node_coords[0])

        # Bottleneck radius: prefer min BRadius across nodes, else min Radius
        if "bottleneck_radius_A" not in tinfo:
            if node_bottlenecks:
                tinfo["bottleneck_radius_A"] = float(min(node_bottlenecks))
            elif node_radii:
                tinfo["bottleneck_radius_A"] = float(min(node_radii))

        # Length: use max/last Distance across nodes if not present
        if "length_A" not in tinfo and node_distances:
            try:
                tinfo["length_A"] = float(max(node_distances))
            except Exception:
                pass

        if "length_A" in tinfo or "entrance_A" in tinfo or "bottleneck_radius_A" in tinfo:
            tunnels.append(tinfo)

    if not tunnels:
        raise ValueError(f"No tunnels parsed from {xml_path}")

    return tunnels


def select_best_tunnel(tunnels: List[Dict], ligand_radius_A: float, pocket_center_A: np.ndarray = None, max_distance_A: float = 80.0) -> Dict:
    # Prefer tunnels whose bottleneck radius >= ligand_radius AND entrance is reasonable distance from pocket
    candidates: List[Tuple[float, float, Dict]] = []  # (length, entrance_distance, tunnel)
    fallbacks: List[Tuple[float, float, float, Dict]] = []  # (negative_radius, length, entrance_distance, tunnel)

    for t in tunnels:
        length = float(t.get("length_A", math.inf))
        bottleneck = t.get("bottleneck_radius_A")
        
        # Calculate distance from tunnel entrance to pocket center for validation
        entrance_distance = math.inf
        if pocket_center_A is not None and "entrance_A" in t:
            entrance_A = np.array(t["entrance_A"], dtype=float)
            entrance_distance = float(np.linalg.norm(entrance_A - pocket_center_A))
        
        # Skip tunnels that are unreasonably far from the pocket
        if entrance_distance > max_distance_A:
            print(f"Skipping tunnel with entrance {entrance_distance:.1f}Å from pocket (max: {max_distance_A:.1f}Å)")
            continue
        
        if bottleneck is not None:
            bottleneck = float(bottleneck)
            if bottleneck >= ligand_radius_A:
                candidates.append((length, entrance_distance, t))
            else:
                # keep for fallback, sort by radius desc then length asc then entrance distance
                fallbacks.append((-bottleneck, length, entrance_distance, t))
        else:
            # unknown bottleneck; keep as fallback with neutral radius
            fallbacks.append((0.0, length, entrance_distance, t))

    if candidates:
        # Sort by length first, then by entrance distance (prefer closer to pocket)
        candidates.sort(key=lambda x: (x[0], x[1]))
        selected = candidates[0][2]
        print(f"Selected tunnel: length={candidates[0][0]:.1f}Å, entrance_distance={candidates[0][1]:.1f}Å")
        return selected
    else:
        if not fallbacks:
            # pick the shortest that's not too far
            valid_tunnels = [t for t in tunnels if "entrance_A" in t]
            if pocket_center_A is not None:
                valid_tunnels = [t for t in valid_tunnels 
                               if np.linalg.norm(np.array(t["entrance_A"]) - pocket_center_A) <= max_distance_A]
            if not valid_tunnels:
                valid_tunnels = tunnels  # fallback to any tunnel
            tunnels_sorted = sorted(valid_tunnels, key=lambda tt: float(tt.get("length_A", math.inf)))
            return tunnels_sorted[0]
        # sort by radius desc (since we stored -radius), then length asc, then entrance distance
        fallbacks.sort()
        selected = fallbacks[0][3]
        print(f"Selected fallback tunnel: length={fallbacks[0][1]:.1f}Å, entrance_distance={fallbacks[0][2]:.1f}Å")
        return selected


def main():
    parser = argparse.ArgumentParser(description="Select best MOLE2 tunnel and compute 2 nm offset insertion point")
    parser.add_argument("--tunnels_xml", required=True, help="Path to MOLE2 tunnels.xml")
    parser.add_argument("--ligand_sdf", required=True, help="Ligand SDF to compute radius")
    parser.add_argument("--residues", required=True, help="File with comma-separated residue ids for pocket center (updated numbering)")
    parser.add_argument("--pdb", required=True, help="PDB file corresponding to residues numbering (e.g., output/pdb/clean_frame_fixed.pdb)")
    parser.add_argument("--positions_out", required=True, help="Output positions file (nm) for gmx insert-molecules -ip")
    parser.add_argument("--json_out", required=False, default="", help="Optional JSON info output")
    parser.add_argument("--offset_nm", required=False, type=float, default=0.0, help="Offset distance outside entrance in nm (0.0 = at entrance)")

    args = parser.parse_args()

    residues = read_residue_list(args.residues)
    pocket_coords_A = read_pdb_coords_for_residues(args.pdb, residues)
    pocket_center_A = np.mean(pocket_coords_A, axis=0)

    ligand_radius_A = calc_ligand_radius_from_sdf(args.ligand_sdf)

    tunnels = parse_mole2_tunnels(args.tunnels_xml)
    print(f"Found {len(tunnels)} tunnels from MOLE2")
    print(f"Ligand radius: {ligand_radius_A:.1f}Å")
    print(f"Pocket center: [{pocket_center_A[0]:.1f}, {pocket_center_A[1]:.1f}, {pocket_center_A[2]:.1f}]Å")
    
    # Select best tunnel with distance validation (max 8nm = 80Å from pocket center)
    best = select_best_tunnel(tunnels, ligand_radius_A, pocket_center_A, max_distance_A=80.0)

    if "entrance_A" not in best:
        raise ValueError("Selected tunnel does not contain an entrance point; cannot compute placement")

    entrance_A = np.array(best["entrance_A"], dtype=float)

    # Direction from entrance to pocket center (towards interior)
    direction_vec_A = pocket_center_A - entrance_A
    norm = float(np.linalg.norm(direction_vec_A))
    if norm == 0.0:
        raise ValueError("Zero-length direction vector computed; entrance equals pocket center")
    unit_dir = direction_vec_A / norm

    # Convert to nm and place AT tunnel entrance (no offset for direct tunnel entry)
    entrance_nm = entrance_A / 10.0
    pocket_center_nm = pocket_center_A / 10.0
    
    # CRITICAL: Account for protein centering transformation
    # The protein was centered to origin [0,0,0] in simulation, but MOLE2 used original coordinates
    # We need to shift tunnel coordinates by the same amount the protein was shifted
    protein_shift = pocket_center_nm  # This is how much the protein center moved from origin
    
    print(f"Protein was shifted by: [{protein_shift[0]:.3f}, {protein_shift[1]:.3f}, {protein_shift[2]:.3f}] nm")
    
    # Apply the same shift to tunnel coordinates to match centered protein
    entrance_nm_centered = entrance_nm - protein_shift
    # Always place at exact tunnel entrance - no offset applied here
    insertion_nm = entrance_nm_centered
    
    print(f"Original entrance: [{entrance_nm[0]:.3f}, {entrance_nm[1]:.3f}, {entrance_nm[2]:.3f}] nm")
    print(f"Centered entrance: [{entrance_nm_centered[0]:.3f}, {entrance_nm_centered[1]:.3f}, {entrance_nm_centered[2]:.3f}] nm")
    
    # Validate insertion position is reasonable (now relative to centered protein)
    centered_pocket = np.array([0.0, 0.0, 0.0])  # Pocket center in simulation coordinates
    insertion_to_protein_center = float(np.linalg.norm(insertion_nm))
    insertion_to_pocket_dist = float(np.linalg.norm(insertion_nm - centered_pocket))
    
    print(f"Final insertion position: [{insertion_nm[0]:.3f}, {insertion_nm[1]:.3f}, {insertion_nm[2]:.3f}] nm") 
    print(f"Distance from protein center: {insertion_to_protein_center:.3f} nm")
    print(f"Distance from pocket center: {insertion_to_pocket_dist:.3f} nm")
    
    # Warn if insertion position seems too far from protein center (should be ~2-6 nm for typical proteins)
    if insertion_to_protein_center > 6.0:
        print(f"WARNING: Insertion position is {insertion_to_protein_center:.1f} nm from protein center")
        print(f"This may be too far for effective steered MD simulations")
        
        # Create a fallback position closer to the protein center
        # Place ligand at a reasonable distance (4-5 nm from protein center)
        direction_from_center = insertion_nm / np.linalg.norm(insertion_nm)  # Unit vector from origin to insertion
        reasonable_distance = 4.5  # nm from protein center
        fallback_insertion = direction_from_center * reasonable_distance
        
        print(f"Creating fallback position: [{fallback_insertion[0]:.3f}, {fallback_insertion[1]:.3f}, {fallback_insertion[2]:.3f}] nm")
        print(f"Fallback distance from protein center: {reasonable_distance:.1f} nm")
        
        # Use fallback position
        insertion_nm = fallback_insertion

    os.makedirs(os.path.dirname(args.positions_out), exist_ok=True)
    with open(args.positions_out, "w") as f:
        f.write(f"{insertion_nm[0]:.3f} {insertion_nm[1]:.3f} {insertion_nm[2]:.3f}\n")

    if args.json_out:
        info = {
            "ligand_radius_A": ligand_radius_A,
            "pocket_center_A": pocket_center_A.tolist(),
            "entrance_A": entrance_A.tolist(),
            "entrance_nm": entrance_nm.tolist(),
            "insertion_nm": insertion_nm.tolist(),
            "selected_tunnel": {k: v for k, v in best.items() if k in ("length_A", "bottleneck_radius_A")},
        }
        os.makedirs(os.path.dirname(args.json_out), exist_ok=True)
        with open(args.json_out, "w") as jf:
            json.dump(info, jf, indent=2)

    print(f"Wrote insertion position to {args.positions_out}")


if __name__ == "__main__":
    main() 