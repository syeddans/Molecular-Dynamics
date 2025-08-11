#!/usr/bin/env python3

import argparse
import os
import re
import sys
from typing import List, Tuple, Dict, Optional

import numpy as np


AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL", "SEC", "PYL"
}

WATER_AND_IONS = {"SOL", "WAT", "HOH", "NA", "CL", "K", "MG", "CA"}


def parse_residue_list(path: str) -> List[int]:
    if not os.path.exists(path):
        raise FileNotFoundError(f"Residue list not found: {path}")
    text = open(path, "r").read().strip()
    if not text:
        return []
    parts = re.split(r"[\s,]+", text)
    residues: List[int] = []
    for p in parts:
        if not p:
            continue
        try:
            residues.append(int(p))
        except ValueError:
            continue
    return residues


def _parse_box_matrix_from_line(line: str) -> np.ndarray:
    """Parse the last .gro line into a 3x3 triclinic box matrix (nm).
    Accepts 3 or 9 floats. Mapping for 9 floats follows GROMACS order:
    vxx vyy vzz vxy vxz vyx vyz vzx vzy, and box matrix B is
    [[vxx, vyx, vzx], [vxy, vyy, vzy], [vxz, vyz, vzz]].
    """
    vals = [float(x) for x in line.strip().split()]
    if len(vals) >= 9:
        vxx, vyy, vzz, vxy, vxz, vyx, vyz, vzx, vzy = vals[:9]
        B = np.array(
            [
                [vxx, vyx, vzx],
                [vxy, vyy, vzy],
                [vxz, vyz, vzz],
            ],
            dtype=float,
        )
    elif len(vals) >= 3:
        vxx, vyy, vzz = vals[:3]
        B = np.diag([vxx, vyy, vzz]).astype(float)
    else:
        # Fallback to 10 nm cube
        B = np.diag([10.0, 10.0, 10.0]).astype(float)
    return B


def read_gro_atoms_and_box(gro_path: str) -> Tuple[List[Dict], np.ndarray]:
    if not os.path.exists(gro_path):
        raise FileNotFoundError(f"GRO not found: {gro_path}")

    with open(gro_path, "r") as f:
        lines = f.readlines()

    if len(lines) < 3:
        raise ValueError(f"Not enough lines in GRO file: {gro_path}")

    try:
        num_atoms = int(lines[1].strip())
    except Exception as e:
        raise ValueError(f"Failed to read atom count from {gro_path}") from e

    atoms: List[Dict] = []
    start = 2
    end = min(2 + num_atoms, len(lines) - 1)  # last line is box

    for i in range(start, end):
        line = lines[i]
        if len(line) < 44:
            continue
        try:
            resid_str = line[0:5].strip()
            resname = line[5:10].strip()
            atomname = line[10:15].strip()
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
            resid = int(resid_str)
        except Exception:
            continue
        atoms.append({
            "resid": resid,
            "resname": resname,
            "atom": atomname,
            "coords": np.array([x, y, z], dtype=float)
        })

    # Box matrix (nm)
    box_line = lines[-1]
    B = _parse_box_matrix_from_line(box_line)

    return atoms, B


def split_system_atoms(atoms: List[Dict]) -> Tuple[List[Dict], List[Dict], List[Dict]]:
    ligand_atoms: List[Dict] = []
    protein_atoms: List[Dict] = []
    other_atoms: List[Dict] = []

    for a in atoms:
        resname = a["resname"].upper()
        if resname == "MOL":
            ligand_atoms.append(a)
        elif resname in AMINO_ACIDS:
            protein_atoms.append(a)
        elif resname in WATER_AND_IONS:
            other_atoms.append(a)
        else:
            if len(resname) == 3 and resname.isalpha() and resname.isupper():
                protein_atoms.append(a)
            else:
                other_atoms.append(a)

    return ligand_atoms, protein_atoms, other_atoms


def select_active_site_atoms(atoms: List[Dict], residue_ids: List[int]) -> List[Dict]:
    target = set(residue_ids)
    return [a for a in atoms if a["resid"] in target]


def _wrap_min_image(diff_frac: np.ndarray) -> np.ndarray:
    """Wrap fractional coordinate differences to [-0.5, 0.5)."""
    return diff_frac - np.round(diff_frac)


def compute_min_distance_pbc(set_a: np.ndarray, set_b: np.ndarray, box: np.ndarray) -> float:
    if set_a.size == 0 or set_b.size == 0:
        return float("nan")
    try:
        invB = np.linalg.inv(box)
    except np.linalg.LinAlgError:
        # Degenerate box; fall back to non-PBC
        diff = set_a[:, None, :] - set_b[None, :, :]
        d = np.linalg.norm(diff, axis=-1)
        return float(np.min(d))

    # fractional coordinates
    a_frac = set_a @ invB
    b_frac = set_b @ invB
    diff_frac = a_frac[:, None, :] - b_frac[None, :, :]
    diff_frac = _wrap_min_image(diff_frac)
    diff_cart = diff_frac @ box
    d = np.linalg.norm(diff_cart, axis=-1)
    return float(np.min(d))


def geometric_center(coords: np.ndarray) -> np.ndarray:
    if coords.size == 0:
        return np.array([np.nan, np.nan, np.nan], dtype=float)
    return np.mean(coords, axis=0)


def pbc_distance_point_to_point(p: np.ndarray, q: np.ndarray, box: np.ndarray) -> float:
    if np.any(np.isnan(p)) or np.any(np.isnan(q)):
        return float("nan")
    try:
        invB = np.linalg.inv(box)
    except np.linalg.LinAlgError:
        return float(np.linalg.norm(p - q))
    dp_frac = (p - q) @ invB
    dp_frac = _wrap_min_image(dp_frac)
    dp = dp_frac @ box
    return float(np.linalg.norm(dp))


def compute_metrics(gro: str, residues_file: str) -> Dict[str, float]:
    atoms, box = read_gro_atoms_and_box(gro)
    ligand, protein, _ = split_system_atoms(atoms)

    ligand_coords = np.array([a["coords"] for a in ligand], dtype=float) if ligand else np.empty((0, 3))
    protein_coords = np.array([a["coords"] for a in protein], dtype=float) if protein else np.empty((0, 3))

    active_resids = parse_residue_list(residues_file) if residues_file else []
    active_site_atoms = select_active_site_atoms(protein, active_resids) if active_resids else []
    active_coords = np.array([a["coords"] for a in active_site_atoms], dtype=float) if active_site_atoms else np.empty((0, 3))

    min_lig_prot = compute_min_distance_pbc(ligand_coords, protein_coords, box)
    min_lig_active = compute_min_distance_pbc(ligand_coords, active_coords, box) if active_coords.size else float("nan")

    lig_com = geometric_center(ligand_coords)
    active_com = geometric_center(active_coords)
    com_dist = pbc_distance_point_to_point(lig_com, active_com, box) if active_coords.size else float("nan")

    return {
        "num_ligand_atoms": float(len(ligand)),
        "num_protein_atoms": float(len(protein)),
        "num_active_site_atoms": float(len(active_site_atoms)),
        "min_ligand_protein_nm": min_lig_prot,
        "min_ligand_active_site_nm": min_lig_active,
        "ligand_COM_to_active_site_COM_nm": com_dist,
    }


def format_metrics(step_label: str, gro: str, metrics: Dict[str, float]) -> str:
    lines = [
        f"[pose-metrics] Step: {step_label}",
        f"[pose-metrics] File: {gro}",
        f"[pose-metrics] Ligand atoms          : {int(metrics['num_ligand_atoms'])}",
        f"[pose-metrics] Protein atoms          : {int(metrics['num_protein_atoms'])}",
        f"[pose-metrics] Active-site atoms      : {int(metrics['num_active_site_atoms'])}",
        f"[pose-metrics] Min ligand-protein     : {metrics['min_ligand_protein_nm']:.3f} nm",
    ]
    if not np.isnan(metrics["min_ligand_active_site_nm"]):
        lines.append(f"[pose-metrics] Min ligand-active     : {metrics['min_ligand_active_site_nm']:.3f} nm")
    if not np.isnan(metrics["ligand_COM_to_active_site_COM_nm"]):
        lines.append(f"[pose-metrics] COM(lig)-COM(active) : {metrics['ligand_COM_to_active_site_COM_nm']:.3f} nm")
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description="Compute per-pose ligand distances to protein and active site from a GRO file.")
    parser.add_argument("--gro", required=True, help="Path to the .gro structure file")
    parser.add_argument("--residues", required=False, default="output/autodock_vina/active_site_residues_updated.txt", help="File containing active site residue IDs (one list, comma/space separated)")
    parser.add_argument("--step", required=False, default="", help="Optional step label to display")
    parser.add_argument("--out", required=False, default="", help="Optional output file to write metrics text")

    args = parser.parse_args()

    metrics = compute_metrics(args.gro, args.residues)
    label = args.step if args.step else os.path.basename(args.gro)
    report = format_metrics(label, args.gro, metrics)

    print(report)
    if args.out:
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        with open(args.out, "w") as f:
            f.write(report + "\n")


if __name__ == "__main__":
    main() 