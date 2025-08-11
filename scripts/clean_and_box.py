import sys
from pymol import cmd
import argparse
import math

def rgyrate(selection):
    # Get the atoms for the selection
    model = cmd.get_model(selection).atom
    # Extract the coordinates
    x = [i.coord for i in model]
    # Get the masses
    mass = [i.get_mass() for i in model]
    # Mass-weighted coordinates
    xm = [(m*i, m*j, m*k) for (i, j, k), m in zip(x, mass)]
    # Sum of masses
    tmass = sum(mass)
    # First part of the sum under the sqrt
    rr = sum(mi*i + mj*j + mk*k for (i, j, k), (mi, mj, mk) in zip(x, xm))
    # Second part of the sum under the sqrt
    mm = sum((sum(i)/tmass)**2 for i in zip(*xm))
    # Radius of gyration
    rg = math.sqrt(rr/tmass - mm)
    # Print it...
    print(f"Radius of gyration: {rg}")
    return rg

parser = argparse.ArgumentParser()
parser.add_argument("--input", required=True, help="Input protein PDB file")
parser.add_argument("--clean_output", required=True, help="Output cleaned protein PDB file")
parser.add_argument("--box_output", required=True, help="Output box info text file")
parser.add_argument("--ligand_resn", default="LIG", help="Ligand residue name (default: LIG)")
parser.add_argument("--active_site_residues", help="File with active site residue numbers")
args = parser.parse_args()

cmd.load(args.input)

# Determine box center based on active site residues or ligand
if args.active_site_residues:
    # Read active site residues from file
    with open(args.active_site_residues, 'r') as f:
        content = f.read().strip()
        if content:
            active_site_residues = [int(r.strip()) for r in content.split(',') if r.strip()]
            print(f"Using active site residues: {active_site_residues}")
            
            # Create selection for active site residues
            residue_selection = " or ".join([f"resi {r}" for r in active_site_residues])
            cmd.select("active_site", f"({residue_selection}) and chain A")
            
            # Calculate center of active site
            center = cmd.centerofmass("active_site")
            print(f"Active site center: {center}")
            
            # Calculate radius of gyration for active site
            try:
                cmd.extend("rgyrate", rgyrate)
                radius_of_gyration = cmd.rgyrate("active_site")
            except:
                radius_of_gyration = rgyrate("active_site")
        else:
            print("Warning: Active site residues file is empty, falling back to ligand-based box")
            cmd.select("ligand", f"resn {args.ligand_resn}")
            center = cmd.centerofmass("ligand")
            try:
                cmd.extend("rgyrate", rgyrate)
                radius_of_gyration = cmd.rgyrate("ligand")
            except:
                radius_of_gyration = rgyrate("ligand")
else:
    print("Warning: Active site residues file not found, falling back to ligand-based box")
    cmd.select("ligand", f"resn {args.ligand_resn}")
    center = cmd.centerofmass("ligand")
    try:
        cmd.extend("rgyrate", rgyrate)
        radius_of_gyration = cmd.rgyrate("ligand")
    except:
        radius_of_gyration = rgyrate("ligand")



# Calculate box size based on radius of gyration (typically 2-3 times the radius)
box_size = max(20.0, radius_of_gyration * 3.0)  # Minimum 20 Å, or 3x radius of gyration

# Get coordinates to calculate bounding box
if args.active_site_residues:
    # Use active site coordinates
    coords = cmd.get_model('active_site').get_coord_list()

if coords:
    min_coords = [min(coord[i] for coord in coords) for i in range(3)]
    max_coords = [max(coord[i] for coord in coords) for i in range(3)]
else:
    # Fallback if no coordinates
    min_coords = [c - box_size/2 for c in center]
    max_coords = [c + box_size/2 for c in center]

# Write box information in the format expected by the Snakemake workflow
with open(args.box_output, "w") as f:
    f.write(f"min: {min_coords}\n")
    f.write(f"max: {max_coords}\n")
    f.write(f"center: {center}\n")
    f.write(f"size: [{box_size}, {box_size}, {box_size}]\n")
    f.write(f"radius_of_gyration: {radius_of_gyration}\n")

#cmd.remove("not chain A")
#cmd.remove("alt !A")
cmd.remove("solvent")
# Only remove ligand if we're not using active site residues
if not args.active_site_residues:
    cmd.remove(f"resn {args.ligand_resn}")
#cmd.remove("resn NA+CL+MG+CA")
cmd.remove("hydro")
# Remove problematic elements that AutoDock Vina cannot handle
cmd.remove("elem Mo")  # Molybdenum
cmd.remove("elem W")   # Tungsten
cmd.remove("elem Re")  # Rhenium
cmd.remove("elem Os")  # Osmium
cmd.remove("elem Ir")  # Iridium
cmd.remove("elem Pt")  # Platinum
cmd.remove("elem Au")  # Gold
cmd.remove("elem Hg")  # Mercury
cmd.remove("elem Pb")  # Lead
cmd.remove("elem Bi")  # Bismuth
cmd.remove("elem Po")  # Polonium
cmd.remove("elem At")  # Astatine
cmd.remove("elem Rn")  # Radon
cmd.remove("elem Fr")  # Francium
cmd.remove("elem Ra")  # Radium
cmd.remove("elem Ac")  # Actinium
cmd.remove("elem Th")  # Thorium
cmd.remove("elem Pa")  # Protactinium
cmd.remove("elem U")   # Uranium
# Keep common cofactor ions like Zn, Fe, Mn, Cu, Mg, Ca, Na, K, etc.
cmd.h_add()
cmd.save(args.clean_output)
cmd.quit() 
