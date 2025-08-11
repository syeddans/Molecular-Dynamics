#!/usr/bin/env python3
import argparse
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger

RDLogger.DisableLog('rdApp.*')

def read_gro_coords(gro_path):
    with open(gro_path, 'r') as f:
        lines = f.readlines()
    natoms = int(lines[1].strip())
    atom_lines = lines[2:2+natoms]
    coords = []
    for line in atom_lines:
        x = float(line[20:28])
        y = float(line[28:36])
        z = float(line[36:44])
        coords.append([x, y, z])
    return np.array(coords), lines


def write_gro_coords(gro_lines, new_coords, out_path):
    lines = gro_lines.copy()
    natoms = int(lines[1].strip())
    assert new_coords.shape[0] == natoms
    for i in range(natoms):
        x, y, z = new_coords[i]
        line = lines[2 + i]
        lines[2 + i] = f"{line[:20]}{x:8.3f}{y:8.3f}{z:8.3f}{line[44:]}"
    with open(out_path, 'w') as f:
        f.writelines(lines)


def get_coords_from_mol(mol):
    conf = mol.GetConformer()
    return np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())])


def kabsch(P, Q):
    # Both Nx3, returns R (3x3), t (3,)
    Pc = P.mean(axis=0)
    Qc = Q.mean(axis=0)
    P0 = P - Pc
    Q0 = Q - Qc
    H = P0.T @ Q0
    U, S, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[2, :] *= -1
        R = Vt.T @ U.T
    t = Qc - R @ Pc
    return R, t


def main():
    ap = argparse.ArgumentParser(description='Apply pose orientation to base ligand .gro')
    ap.add_argument('--base_sdf', required=True, help='Base ligand SDF used for parametrization')
    ap.add_argument('--pose_sdf', required=True, help='Pose SDF (converted from Vina pose)')
    ap.add_argument('--base_gro', required=True, help='Base ligand .gro to transform')
    ap.add_argument('--out_gro', required=True, help='Output pose .gro with base atom order')
    args = ap.parse_args()

    base = Chem.MolFromMolFile(args.base_sdf, removeHs=False)
    pose = Chem.MolFromMolFile(args.pose_sdf, removeHs=False)
    if base is None or pose is None:
        raise SystemExit('Failed to read SDF files')

    # Ensure 3D and hydrogens
    if not base.GetConformer().Is3D():
        AllChem.EmbedMolecule(base)
    if not pose.GetConformer().Is3D():
        AllChem.EmbedMolecule(pose)

    # Make working copies
    base_orig = Chem.Mol(base)
    base_to_pose = Chem.Mol(base)

    # Use O3A to align base onto pose (handles mapping)
    try:
        from rdkit.Chem import rdMolAlign
        o3a = rdMolAlign.GetO3A(base_to_pose, pose)
        o3a.Align()
    except Exception:
        # Fallback: MMFF-based alignment
        AllChem.AlignMol(base_to_pose, pose)

    # Compute rigid transform between base_orig and base_to_pose
    P = get_coords_from_mol(base_orig)
    Q = get_coords_from_mol(base_to_pose)
    # Use heavy atoms only for stability
    heavy = [i for i in range(base.GetNumAtoms()) if base.GetAtomWithIdx(i).GetAtomicNum() > 1]
    if len(heavy) >= 3:
        P = P[heavy]
        Q = Q[heavy]
    R, t = kabsch(P, Q)

    # Apply transform (Å) to base_gro coords (nm -> convert to Å, apply, back to nm)
    gro_coords, gro_lines = read_gro_coords(args.base_gro)
    gro_coords_A = gro_coords * 10.0
    transformed_A = (gro_coords_A @ R.T) + t
    transformed_nm = transformed_A / 10.0

    write_gro_coords(gro_lines, transformed_nm, args.out_gro)


if __name__ == '__main__':
    main() 