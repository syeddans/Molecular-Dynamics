#!/usr/bin/env python3
import argparse
import json


def main():
    parser = argparse.ArgumentParser(description="Update pull MDP from tunnel info JSON")
    parser.add_argument("--json", required=True, help="Path to output/ligand/tunnel_info.json")
    parser.add_argument("--template", required=True, help="Template MDP file (e.g., mdp/pull.mdp)")
    parser.add_argument("--output", required=True, help="Output updated MDP file")
    parser.add_argument("--mode", choices=["inward", "outward"], default="inward", help="Direction to pull relative to tunnel entrance")
    args = parser.parse_args()

    with open(args.json, "r") as f:
        info = json.load(f)

    entrance = info["entrance_nm"]
    pocket = info["pocket_center_A"]

    # Convert pocket center from Å to nm
    pocket_nm = [c / 10.0 for c in pocket]

    # Direction unit vector from entrance to pocket (inward)
    # This defines the pulling direction for GROMACS SMD
    # With negative pull rate, ligand will be pulled along this vector (toward binding site)
    dir_vec = [p - e for p, e in zip(pocket_nm, entrance)]
    norm = sum(d * d for d in dir_vec) ** 0.5
    if norm == 0:
        raise SystemExit("Zero direction vector in tunnel info JSON")
    unit = [d / norm for d in dir_vec]

    if args.mode == "outward":
        unit = [-u for u in unit]

    with open(args.template, "r") as f:
        lines = f.readlines()

    updated = False
    vec_str = f"{unit[0]:.3f} {unit[1]:.3f} {unit[2]:.3f}"
    for i, line in enumerate(lines):
        if line.strip().startswith("pull_coord1_vec"):
            lines[i] = f"pull_coord1_vec         = {vec_str}    ; Set from tunnel entrance to pocket center\n"
            updated = True
            break

    if not updated:
        # Append if not present
        lines.append(f"pull_coord1_vec         = {vec_str}    ; Set from tunnel entrance to pocket center\n")

    with open(args.output, "w") as f:
        f.writelines(lines)

    print(f"Updated {args.output} with pull vector {vec_str}")


if __name__ == "__main__":
    main() 