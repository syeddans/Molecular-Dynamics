configfile: "config.yaml"

# Use conda environments automatically


rule all:
    input:
        # Apo simulation outputs
        "output/md/production_apo.xtc",
        "output/md/production_apo.tpr", 
        "output/analysis/rmsd_apo.xvg",
        "output/ligand/ligand.gro",

        # Best pose selection
        "output/best_pose/selection.success",
        
        # Sequential workflow completion
        "output/workflow_complete.txt",

        # Steered MD outputs
        "output/smd/smd.xtc",
        "output/smd/smd.tpr",

        # MM/PBSA final binding affinity results
        "output/mmpbsa/binding_analysis.txt",
        "output/mmpbsa/binding_plot.png",

        # Per-pose orientation screening report
        "output/pose_runs/pose_forces_summary.txt"

# Clean protein by removing ligand and using PDBFixer for robust cleaning
rule create_apo_protein:
    input:
        holo_pdb = config["holo_protein_pdb"]
    output:
        apo_pdb = "output/protein/apo_protein.pdb",
        apo_success = "output/protein/apo_protein.success",
        active_site_residues = "output/protein/active_site_residues.txt"
    conda: "envs/environment.yml"
    shell:
        """
        set -e  # Exit on any error
        
        mkdir -p output/protein
        
        python scripts/create_apo_protein_clean.py \
            --input {input.holo_pdb} \
            --output {output.apo_pdb} \
            --active_site_output {output.active_site_residues}
            
        # Create success marker
        echo "Apo protein creation completed successfully" > {output.apo_success}
        
        """



# Process cleaned apo protein with pdb2gmx (handles both single and multi-chain)
rule process_apo_protein:
    input:
        apo_pdb = "output/protein/apo_protein.pdb",
        apo_success = "output/protein/apo_protein.success"
    output:
        protein_gro = "output/protein/protein.gro",
        protein_top = "output/protein/topol.top",
        protein_success = "output/protein/protein.success"
    params:
        ff = config.get("force_field", "amber99sb-ildn"),
        water = config.get("water_model", "tip3p")
    conda: "envs/environment.yml"
    shell:
        """
        set -e  # Exit on any error
        
        # Run GROMACS pdb2gmx
        cd output/protein
        gmx pdb2gmx -f apo_protein.pdb -o protein.gro -p topol.top -ff {params.ff} -water {params.water} -ignh
        
        # Validate success
        if [ ! -f protein.gro ] || [ ! -f topol.top ]; then
            echo "Error: Protein processing failed - output files not created"
            exit 1
        fi
        
        # Create success marker
        echo "Protein processing completed successfully" > protein.success
        

        """

#  Create simulation box and solvate apo protein
rule solvate_apo:
    input:
        protein_gro = "output/protein/protein.gro",
        protein_top = "output/protein/topol.top",
        protein_success = "output/protein/protein.success"
    output:
        boxed_gro = "output/solvation/apo_boxed.gro",
        solvated_gro = "output/solvation/apo_solvated.gro",
        solvated_top = "output/solvation/topol.top",
        solvate_success = "output/solvation/solvate.success"
    params:
        box_distance = config.get("box_distance", "1.2")
    conda: "envs/environment.yml"
    shell:
        """
        set -e  # Exit on any error
        
        mkdir -p output/solvation
        
        # Start with a clean topology (remove any existing water)
        sed '/^SOL/d' {input.protein_top} > output/solvation/clean_topol.top
        
        # Create box using rhombic dodecahedron for better performance (15% fewer water molecules)
        BOX_TYPE=$(echo "{config[box_type]}" 2>/dev/null || echo "dodecahedron")
        gmx editconf -f {input.protein_gro} -o {output.boxed_gro} -c -d {params.box_distance} -bt $BOX_TYPE
        
        # Solvate using clean topology
        gmx solvate -cp {output.boxed_gro} -cs spc216.gro -o {output.solvated_gro} -p output/solvation/clean_topol.top
        
        # Copy final topology and all required include files
        cp output/solvation/clean_topol.top {output.solvated_top}
        cp output/protein/*.itp output/solvation/ 2>/dev/null || true
        
        # Validate success
        if [ ! -f {output.solvated_gro} ] || [ ! -f {output.solvated_top} ]; then
            echo "Error: Solvation failed - output files not created"
            exit 1
        fi
        
        # Create success marker
        echo "Solvation completed successfully" > {output.solvate_success}
        

        """

# Add ions for neutralization 
rule add_ions_apo:
    input:
        solvated_gro = "output/solvation/apo_solvated.gro",
        solvated_top = "output/solvation/topol.top",
        solvate_success = "output/solvation/solvate.success",
        ions_mdp = "mdp/ions.mdp"
    output:
        ions_tpr = "output/solvation/ions.tpr",
        system_gro = "output/solvation/apo_system.gro",
        system_top = "output/solvation/system_topol.top",
        ions_success = "output/solvation/ions.success"
    conda: "envs/environment.yml"
    shell:
        """
        set -e  # Exit on any error
        
        # Create tpr for genion
        gmx grompp -f {input.ions_mdp} -c {input.solvated_gro} -p {input.solvated_top} -o {output.ions_tpr} 
        
        # Add ions (select SOL when prompted)
        echo "SOL" | gmx genion -s {output.ions_tpr} -o {output.system_gro} -p {input.solvated_top} -pname NA -nname CL -neutral -conc 0.15
        
        # Copy final topology
        cp {input.solvated_top} {output.system_top}
        
        # Validate success
        if [ ! -f {output.system_gro} ] || [ ! -f {output.system_top} ]; then
            echo "Error: Ion addition failed - output files not created"
            exit 1
        fi
        
        # Create success marker
        echo "Ion addition completed successfully" > {output.ions_success}
        

        """

# Energy minimization
rule energy_minimization:
    input:
        structure = "output/solvation/apo_system.gro",
        topology = "output/solvation/system_topol.top",
        ions_success = "output/solvation/ions.success",
        mdp = "mdp/em.mdp"
    output:
        tpr = "output/em/em.tpr",
        gro = "output/em/em.gro",
        edr = "output/em/em.edr",
        em_success = "output/em/em.success"
    conda: "envs/environment.yml"
    shell:
        """
        set -e  # Exit on any error
        
        mkdir -p output/em

        # Prepare EM
        gmx grompp -f {input.mdp} -c {input.structure} -p {input.topology} -o {output.tpr}

        # Run EM - use consistent threading
        export OMP_NUM_THREADS=4
        gmx mdrun -v -deffnm output/em/em -ntmpi 1 -ntomp 4
        
        # Validate success
        if [ ! -f {output.gro} ]; then
            echo "Error: Energy minimization failed - output file not created"
            exit 1
        fi
        
        # Create success marker
        echo "Energy minimization completed successfully" > {output.em_success}

        """

# NVT equilibration (temperature coupling)
rule nvt_equilibration:
    input:
        em_gro = "output/em/em.gro",
        em_success = "output/em/em.success",
        system_top = "output/solvation/system_topol.top",
        nvt_mdp = "mdp/nvt.mdp"
    output:
        nvt_tpr = "output/nvt/nvt.tpr",
        nvt_gro = "output/nvt/nvt.gro",
        nvt_edr = "output/nvt/nvt.edr",
        nvt_cpt = "output/nvt/nvt.cpt",
        nvt_success = "output/nvt/nvt.success"
    conda: "envs/environment.yml"
    shell:
        """
        set -e  # Exit on any error
        
        mkdir -p output/nvt
        
        # Prepare NVT
        gmx grompp -f {input.nvt_mdp} -c {input.em_gro} -r {input.em_gro} -p {input.system_top} -o {output.nvt_tpr} -maxwarn 1 
        
        # Run NVT - GPU-optimized threading
        export OMP_NUM_THREADS=4
        gmx mdrun -v -deffnm output/nvt/nvt -ntmpi 1 -ntomp 4
        
        # Validate success
        if [ ! -f {output.nvt_gro} ]; then
            echo "Error: NVT equilibration failed - output file not created"
            exit 1
        fi
        
        # Create success marker
        echo "NVT equilibration completed successfully" > {output.nvt_success}
        

        """

# NPT equilibration (pressure coupling)
rule npt_equilibration:
    input:
        nvt_gro = "output/nvt/nvt.gro",
        nvt_cpt = "output/nvt/nvt.cpt",
        nvt_success = "output/nvt/nvt.success",
        system_top = "output/solvation/system_topol.top",
        npt_mdp = "mdp/npt.mdp"
    output:
        npt_tpr = "output/npt/npt.tpr",
        npt_gro = "output/npt/npt.gro",
        npt_edr = "output/npt/npt.edr",
        npt_cpt = "output/npt/npt.cpt",
        npt_success = "output/npt/npt.success"
    conda: "envs/environment.yml"
    shell:
        """
        set -e  # Exit on any error
        
        mkdir -p output/npt
        
        # Prepare NPT  
        gmx grompp -f {input.npt_mdp} -c {input.nvt_gro} -t {input.nvt_cpt} -r {input.nvt_gro} -p {input.system_top} -o {output.npt_tpr} -maxwarn 1 
        
        # Run NPT - GPU-optimized threading
        export OMP_NUM_THREADS=4
        gmx mdrun -v -deffnm output/npt/npt -ntmpi 1 -ntomp 4
        
        # Validate success
        if [ ! -f {output.npt_gro} ]; then
            echo "Error: NPT equilibration failed - output file not created"
            exit 1
        fi
        
        # Create success marker
        echo "NPT equilibration completed successfully" > {output.npt_success}
        

        """

# Production MD simulation (GPU-optimized long run)
rule production_md_apo:
    input:
        npt_gro = "output/npt/npt.gro",
        npt_cpt = "output/npt/npt.cpt",
        npt_success = "output/npt/npt.success",
        system_top = "output/solvation/system_topol.top",
        prod_mdp = "mdp/production.mdp"
    output:
        prod_tpr = "output/md/production_apo.tpr",
        prod_xtc = "output/md/production_apo.xtc",
        prod_gro = "output/md/production_apo.gro",
        prod_edr = "output/md/production_apo.edr",
        prod_log = "output/md/production_apo.log",
        prod_cpt = "output/md/production_apo.cpt",
        prod_success = "output/md/production_apo.success"
    conda: "envs/environment.yml"
    shell:
        """
        set -e  # Exit on any error
        
        mkdir -p output/md
        
        # Prepare production MD
        gmx grompp -f {input.prod_mdp} -c {input.npt_gro} -t {input.npt_cpt} -p {input.system_top} -o {output.prod_tpr} 
        
        # Run production MD - GPU-optimized threading
        export OMP_NUM_THREADS=4
        gmx mdrun -v -deffnm output/md/production_apo -ntmpi 1 -ntomp 4
        
        # Validate success
        if [ ! -f {output.prod_xtc} ]; then
            echo "Error: Production MD failed - output file not created"
            exit 1
        fi
        
        # Create success marker
        echo "Production MD completed successfully" > {output.prod_success}
        

        """


# Basic analysis - RMSD
rule analyze_rmsd:
    input:
        prod_xtc = "output/md/production_apo.xtc",
        prod_tpr = "output/md/production_apo.tpr",
        prod_success = "output/md/production_apo.success"
    output:
        rmsd_xvg = "output/analysis/rmsd_apo.xvg",
        rmsd_plot = "output/analysis/rmsd_apo_plot.png",
        analysis_success = "output/analysis/analysis.success"
    conda: "envs/environment.yml"
    shell:
        """
        set -e  # Exit on any error
        
        mkdir -p output/analysis
        
        # Calculate RMSD (Backbone selection: 4)
        echo "4 4" | gmx rms -s {input.prod_tpr} -f {input.prod_xtc} -o {output.rmsd_xvg} -tu ns
        python scripts/plot_rmsd.py 
        
        # Validate success
        if [ ! -f {output.rmsd_xvg} ]; then
            echo "Error: Analysis failed - output file not created"
            exit 1
        fi
        
        # Create success marker
        echo "Analysis completed successfully" > {output.analysis_success}
        """



rule prepare_ligand:
    input:
        ligand_sdf = config["new_ligand_file"],
        protein_gro = "output/md/production_apo.gro",
        original_pdb = config["holo_protein_pdb"],
        analysis_success = "output/analysis/analysis.success"
    output:
        ligand_gro = "output/ligand/ligand.gro",
        ligand_itp = "output/ligand/ligand.itp",
        ligand_success = "output/ligand/ligand.success"
    conda: "envs/environment.yml"
    shell:
        """
        set -e
        mkdir -p output/ligand
        python scripts/setup_amber_ligand_complex.py --ligand {input.ligand_sdf} --output-dir output/ligand
        echo "Ligand setup completed successfully" > {output.ligand_success}
        """

# Create clean apo protein structure (remove water, ions, and non-protein components)
rule create_clean_apo:
    input:
        protein_gro = "output/md/production_apo.gro",
        production_top = "output/solvation/system_topol.top"
    output:
        clean_apo_gro = "output/protein/clean_apo.gro",
        clean_apo_top = "output/protein/clean_apo.top",
        clean_apo_success = "output/protein/clean_apo.success"
    conda: "envs/environment.yml"
    shell:
        """
        set -e
        
        # Extract only protein atoms (group 1) and exclude water, ions, and other molecules
        # Group 1 is the Protein group as shown in the make_ndx output
        echo "1" | gmx trjconv -s {input.protein_gro} -f {input.protein_gro} -o {output.clean_apo_gro}
        
        # Create clean topology by removing water and ions from the production topology
        sed '/^SOL/d' {input.production_top} | sed '/^NA/d' | sed '/^CL/d' > {output.clean_apo_top}
        
        echo "Clean apo protein creation completed successfully" > {output.clean_apo_success}
        """

rule autodock_vina: 
    input: 
        ligand_sdf = config["new_ligand_file"],
        production_gro = "output/md/production_apo.gro",
        production_tpr = "output/md/production_apo.tpr",
        production_success = "output/md/production_apo.success",
        original_pdb = config["holo_protein_pdb"],
        active_site_residues = "output/protein/active_site_residues.txt"
    output:
        vina_success = "output/autodock_vina/vina.success",
        original_residues = "output/autodock_vina/active_site_residues_original.txt",
        updated_residues = "output/autodock_vina/active_site_residues_updated.txt"
    conda: "envs/environment.yml"
    shell:
        """
        set -e  # Exit on any error
        
        echo "=== DEBUG: Starting AutoDock Vina workflow ==="
        
        mkdir -p output/autodock_vina
        mkdir -p output/pdb
        
        echo "=== DEBUG: Converting GRO to PDB ==="
        # Convert GRO to PDB format for AutoDock Vina
        echo "1\n1" | gmx trjconv -f {input.production_gro} -s {input.production_tpr} -o output/pdb/clean_frame.pdb -pbc mol -center
        
        echo "=== DEBUG: Fixing PDB structure ==="
        # Fix the PDB structure using PDBFixer
        python scripts/fix_pdb_for_caver.py --input output/pdb/clean_frame.pdb --output output/pdb/clean_frame_fixed.pdb
        
        echo "=== DEBUG: Converting residue numbers from original to GROMACS PDB ==="
        # Keep original residue numbers for MOLE2 (which uses original PDB coordinate system)
        cp {input.active_site_residues} output/autodock_vina/active_site_residues_original.txt
        
        # Convert residue numbers from original PDB to GROMACS-processed PDB
        python scripts/convert_residue_numbers.py \
            --original_pdb {input.original_pdb} \
            --gromacs_pdb output/pdb/clean_frame_fixed.pdb \
            --active_site_residues {input.active_site_residues} \
            --output output/autodock_vina/active_site_residues_updated.txt
        
        echo "=== DEBUG: Reading converted residues ==="
        residues=$(cat output/autodock_vina/active_site_residues_updated.txt)
        echo "=== DEBUG: Converted active site residues: $residues ==="
        
        # Convert ligand to PDBQT format
        obabel {input.ligand_sdf} -O output/autodock_vina/ligand.pdbqt --gen3d
        
        
        # Create binding box based on active site residues
        python scripts/clean_and_box.py \
            --input output/pdb/clean_frame_fixed.pdb \
            --clean_output output/autodock_vina/clean_protein.pdb \
            --box_output output/autodock_vina/box.txt \
            --active_site_residues output/autodock_vina/active_site_residues_updated.txt
        
        # Skip the problematic prepare_pdb_split_alt_confs.py and use prepare_receptor4.py directly
        # The prepare_receptor4.py script can handle the protein preparation without the split step
        mgltools/mgltools_x86_64Linux2_1.5.6/bin/pythonsh mgltools/mgltools_x86_64Linux2_1.5.6/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r output/autodock_vina/clean_protein.pdb -o output/autodock_vina/clean_protein.pdbqt
        
        center_line=$(grep "center:" output/autodock_vina/box.txt)
        center_coords=$(echo $center_line | sed 's/center: \[//; s/\]//; s/,/ /g')
        center_x=$(echo $center_coords | cut -d' ' -f1)
        center_y=$(echo $center_coords | cut -d' ' -f2)
        center_z=$(echo $center_coords | cut -d' ' -f3)
        
        size_line=$(grep "size:" output/autodock_vina/box.txt)
        size_coords=$(echo $size_line | sed 's/size: \[//; s/\]//; s/,/ /g')
        size_x=$(echo $size_coords | cut -d' ' -f1)
        size_y=$(echo $size_coords | cut -d' ' -f2)
        size_z=$(echo $size_coords | cut -d' ' -f3)
        
        # Run vina docking
        autodock_vina_1_1_2_linux_x86/bin/vina --receptor output/autodock_vina/clean_protein.pdbqt \
             --ligand output/autodock_vina/ligand.pdbqt \
             --center_x $center_x --center_y $center_y --center_z $center_z \
             --size_x $size_x --size_y $size_y --size_z $size_z \
             --out output/autodock_vina/docked.pdbqt \
             --log output/autodock_vina/vina.log

        echo "Autodock VINA completed successfully" > {output.vina_success}
        """

# Split Vina multi-model PDBQT into per-pose SDFs
checkpoint split_vina_poses:
    input:
        vina_success = "output/autodock_vina/vina.success"
    output:
        directory("output/autodock_vina/poses")
    conda: "envs/environment.yml"
    shell:
        """
        set -e
        rm -rf output/autodock_vina/poses
        mkdir -p output/autodock_vina/poses
        # Split models: pose_1.sdf, pose_2.sdf, ...
        obabel output/autodock_vina/docked.pdbqt -O output/autodock_vina/poses/pose_.sdf -m
        # Normalize consistent naming pose_1.sdf, pose_2.sdf ... already done by obabel
        ls -1 output/autodock_vina/poses/*.sdf | nl -v 1 > output/autodock_vina/poses/index.txt
        """

# Helper to get dynamic list of poses after checkpoint
def pose_pullfs_input(wildcards):
    ck = checkpoints.split_vina_poses.get(**wildcards)
    pose_dir = ck.output[0]
    import glob, os
    sdf_files = sorted(glob.glob(os.path.join(pose_dir, "pose_*.sdf")))
    return [f"output/pose_runs/{os.path.splitext(os.path.basename(p))[0]}/smd/pullf.xvg" for p in sdf_files]

# Prepare per-pose ligand parameters (GROMACS) in pose directory
# Per-pose ligand parameters are not needed; reuse global output/ligand/clean_ligand.itp
# and only orient geometry per pose.

# Prepare oriented ligand .gro per pose by applying pose transform to base ligand
rule prepare_pose_gro:
    input:
        base_sdf = config["new_ligand_file"],
        pose_sdf = "output/autodock_vina/poses/{pose}.sdf",
        base_gro = "output/ligand/ligand.gro"
    output:
        pose_gro = "output/poses/{pose}.gro"
    conda: "envs/environment.yml"
    shell:
        """
        set -e
        mkdir -p output/poses
        python scripts/apply_pose_to_ligand.py --base_sdf {input.base_sdf} --pose_sdf {input.pose_sdf} --base_gro {input.base_gro} --out_gro {output.pose_gro}
        """

# Build per-pose complex at tunnel entrance using the oriented pose geometry
rule create_pose_complex:
    input:
        clean_apo_gro = "output/protein/clean_apo.gro",
        clean_apo_top = "output/protein/clean_apo.top",
        placement = "output/ligand/ligand_offset_vector.dat",
        pose_gro = "output/poses/{pose}.gro",
        ligand_itp = "output/ligand/ligand.itp"
    output:
        complex_gro = "output/pose_runs/{pose}/complex/complex.gro",
        complex_top = "output/pose_runs/{pose}/complex/complex.top",
        complex_success = "output/pose_runs/{pose}/complex/complex.success",
        pose_clean_ligand_itp = "output/pose_runs/{pose}/complex/clean_ligand.itp",
        pose_posre_itp = "output/pose_runs/{pose}/complex/posre.itp"
    conda: "envs/environment.yml"
    shell:
        """
        set -e
        mkdir -p output/pose_runs/{wildcards.pose}/complex
        
        # Center protein and create pose-specific ligand positioned at tunnel entrance
        gmx editconf -f {input.clean_apo_gro} -o output/pose_runs/{wildcards.pose}/complex/clean_apo_centered.gro -c -center 0 0 0
        
        # CRITICAL FIX: Move the pose-oriented ligand to the tunnel entrance position
        # The pose_gro has the correct orientation but wrong position (Vina binding site)
        # We need to translate it to the tunnel entrance coordinates
        python scripts/place_pose_at_tunnel.py \
            --pose_gro {input.pose_gro} \
            --tunnel_position {input.placement} \
            --output output/pose_runs/{wildcards.pose}/complex/pose_at_tunnel.gro
        
        # Directly combine protein and ligand coordinates to preserve exact positioning
        python scripts/combine_protein_ligand.py \
            --protein output/pose_runs/{wildcards.pose}/complex/clean_apo_centered.gro \
            --ligand output/pose_runs/{wildcards.pose}/complex/pose_at_tunnel.gro \
            --output {output.complex_gro}

        # Build per-pose complex topology using cleaned apo top and RAW ligand itp to inject atomtypes
        python scripts/create_complex_topology.py --protein-top {input.clean_apo_top} --ligand-itp {input.ligand_itp} --output {output.complex_top}

        # Generate a cleaned ligand ITP (without atomtypes/defaults) in the pose directory as referenced by the topology
        python scripts/fix_ligand_topology.py --input {input.ligand_itp} --output {output.pose_clean_ligand_itp} --add-posres

        # Ensure posre include resolves from the pose directory
        cp output/solvation/posre.itp {output.pose_posre_itp}

        # Make includes absolute to avoid CWD issues in parallel
        sed -i "s|#include \"clean_ligand.itp\"|#include \"$(pwd)/output/pose_runs/{wildcards.pose}/complex/clean_ligand.itp\"|" {output.complex_top}
        # Note: posre.itp path will be handled by the topology creation script

        echo "Protein-ligand complex creation completed." > {output.complex_success}

        # Pose metrics for complex structure
        python scripts/pose_distance_metrics.py \
            --gro {output.complex_gro} \
            --residues output/autodock_vina/active_site_residues_updated.txt \
            --step create_pose_complex \
            --out output/pose_runs/{wildcards.pose}/complex/metrics_complex.txt
        """

rule solvate_pose:
    input:
        complex_gro = "output/pose_runs/{pose}/complex/complex.gro",
        complex_top = "output/pose_runs/{pose}/complex/complex.top",
        pose_clean_ligand_itp = "output/pose_runs/{pose}/complex/clean_ligand.itp"
    output:
        solv_gro = "output/pose_runs/{pose}/solvation/solv.gro",
        solv_top = "output/pose_runs/{pose}/solvation/solv.top"
    params:
        box_distance = config.get("box_distance", "1.2")
    conda: "envs/environment.yml"
    shell:
        """
        set -e
        mkdir -p output/pose_runs/{wildcards.pose}/solvation
        
        BOX_TYPE=$(echo "{config[box_type]}" 2>/dev/null || echo "dodecahedron")
        
        # CRITICAL FIX: Use -noc to prevent automatic centering during box creation
        # This preserves the carefully calculated ligand-protein relative positions
        gmx editconf -f {input.complex_gro} -o output/pose_runs/{wildcards.pose}/solvation/boxed.gro -d {params.box_distance} -bt $BOX_TYPE -noc
        
        gmx solvate -cp output/pose_runs/{wildcards.pose}/solvation/boxed.gro -cs spc216.gro -o {output.solv_gro} -p {input.complex_top}
        
        # Fix include paths in the topology for the new location
        sed 's|#include "posre.itp"|#include "../complex/posre.itp"|g; s|#include "clean_ligand.itp"|#include "../complex/clean_ligand.itp"|g' {input.complex_top} > {output.solv_top}

        # Pose metrics after solvation
        python scripts/pose_distance_metrics.py \
            --gro {output.solv_gro} \
            --residues output/autodock_vina/active_site_residues_updated.txt \
            --step solvate_pose \
            --out output/pose_runs/{wildcards.pose}/solvation/metrics_solvate.txt
        """

rule ions_pose:
    input:
        solv_gro = "output/pose_runs/{pose}/solvation/solv.gro",
        solv_top = "output/pose_runs/{pose}/solvation/solv.top",
        ions_mdp = "mdp/ions_complex.mdp"
    output:
        ion_gro = "output/pose_runs/{pose}/ions/ion.gro",
        ion_top = "output/pose_runs/{pose}/ions/ion.top"
    conda: "envs/environment.yml"
    shell:
        """
        set -e
        mkdir -p output/pose_runs/{wildcards.pose}/ions
        
        gmx grompp -f {input.ions_mdp} -c {input.solv_gro} -r {input.solv_gro} -p {input.solv_top} -o output/pose_runs/{wildcards.pose}/ions/temp.tpr -maxwarn 1000
        echo "SOL" | gmx genion -s output/pose_runs/{wildcards.pose}/ions/temp.tpr -o {output.ion_gro} -p {input.solv_top} -pname NA -nname CL -neutral
        cp {input.solv_top} {output.ion_top}
        rm -f output/pose_runs/{wildcards.pose}/ions/temp.tpr

        # Pose metrics after ion addition
        python scripts/pose_distance_metrics.py \
            --gro {output.ion_gro} \
            --residues output/autodock_vina/active_site_residues_updated.txt \
            --step ions_pose \
            --out output/pose_runs/{wildcards.pose}/ions/metrics_ions.txt
        """

rule em_pose:
    input:
        ion_gro = "output/pose_runs/{pose}/ions/ion.gro",
        ion_top = "output/pose_runs/{pose}/ions/ion.top",
        em_mdp = "mdp/em_complex.mdp"
    output:
        em_gro = "output/pose_runs/{pose}/em/em.gro",
        em_tpr = "output/pose_runs/{pose}/em/em.tpr",
        em_log = "output/pose_runs/{pose}/em/em.log",
        em_vis_pdb = "output/pose_runs/{pose}/em/em_vis.pdb"
    conda: "envs/environment.yml"
    shell:
        """
        set -e
        mkdir -p output/pose_runs/{wildcards.pose}/em
        
        gmx grompp -f {input.em_mdp} -c {input.ion_gro} -r {input.ion_gro} -p {input.ion_top} -o {output.em_tpr} -maxwarn 1000
        export OMP_NUM_THREADS=4
        gmx mdrun -v -deffnm output/pose_runs/{wildcards.pose}/em/em -ntmpi 1 -ntomp 4

        # Create visualization PDB with proper PBC correction
        printf "1\n0\n" | gmx trjconv -s {output.em_tpr} -f {output.em_gro} -o {output.em_vis_pdb} -pbc mol -center -ur compact
        
        # Pose metrics after EM
        python scripts/pose_distance_metrics.py \
            --gro {output.em_gro} \
            --residues output/autodock_vina/active_site_residues_updated.txt \
            --step em_pose \
            --out output/pose_runs/{wildcards.pose}/em/metrics_em.txt
        """

rule nvt_pose:
    input:
        em_gro = "output/pose_runs/{pose}/em/em.gro",
        ion_top = "output/pose_runs/{pose}/ions/ion.top",
        nvt_mdp = "mdp/nvt_complex.mdp"
    output:
        nvt_gro = "output/pose_runs/{pose}/nvt/nvt.gro",
        nvt_tpr = "output/pose_runs/{pose}/nvt/nvt.tpr",
        nvt_log = "output/pose_runs/{pose}/nvt/nvt.log",
        nvt_cpt = "output/pose_runs/{pose}/nvt/nvt.cpt",
        nvt_vis_pdb = "output/pose_runs/{pose}/nvt/nvt_vis.pdb"
    conda: "envs/environment.yml"
    shell:
        """
        set -e
        mkdir -p output/pose_runs/{wildcards.pose}/nvt
        
        gmx grompp -f {input.nvt_mdp} -c {input.em_gro} -r {input.em_gro} -p {input.ion_top} -o {output.nvt_tpr} -maxwarn 1000
        export OMP_NUM_THREADS=4
        gmx mdrun -v -deffnm output/pose_runs/{wildcards.pose}/nvt/nvt -ntmpi 1 -ntomp 4

        # Create visualization PDB with proper PBC correction
        printf "1\n0\n" | gmx trjconv -s {output.nvt_tpr} -f {output.nvt_gro} -o {output.nvt_vis_pdb} -pbc mol -center -ur compact
        
        # Pose metrics after NVT
        python scripts/pose_distance_metrics.py \
            --gro {output.nvt_gro} \
            --residues output/autodock_vina/active_site_residues_updated.txt \
            --step nvt_pose \
            --out output/pose_runs/{wildcards.pose}/nvt/metrics_nvt.txt
        """

rule npt_pose:
    input:
        nvt_gro = "output/pose_runs/{pose}/nvt/nvt.gro",
        nvt_cpt = "output/pose_runs/{pose}/nvt/nvt.cpt",
        ion_top = "output/pose_runs/{pose}/ions/ion.top",
        npt_mdp = "mdp/npt_complex.mdp"
    output:
        npt_gro = "output/pose_runs/{pose}/npt/npt.gro",
        npt_tpr = "output/pose_runs/{pose}/npt/npt.tpr",
        npt_log = "output/pose_runs/{pose}/npt/npt.log",
        npt_cpt = "output/pose_runs/{pose}/npt/npt.cpt",
        npt_vis_pdb = "output/pose_runs/{pose}/npt/npt_vis.pdb"
    conda: "envs/environment.yml"
    shell:
        """
        set -e
        mkdir -p output/pose_runs/{wildcards.pose}/npt
        
        gmx grompp -f {input.npt_mdp} -c {input.nvt_gro} -r {input.nvt_gro} -t {input.nvt_cpt} -p {input.ion_top} -o {output.npt_tpr} -maxwarn 1000
        export OMP_NUM_THREADS=4
        gmx mdrun -v -deffnm output/pose_runs/{wildcards.pose}/npt/npt -ntmpi 1 -ntomp 4

        # Create visualization PDB with proper PBC correction
        printf "1\n0\n" | gmx trjconv -s {output.npt_tpr} -f {output.npt_gro} -o {output.npt_vis_pdb} -pbc mol -center -ur compact
        
        # Pose metrics after NPT
        python scripts/pose_distance_metrics.py \
            --gro {output.npt_gro} \
            --residues output/autodock_vina/active_site_residues_updated.txt \
            --step npt_pose \
            --out output/pose_runs/{wildcards.pose}/npt/metrics_npt.txt
        """

rule short_smd_pose:
    input:
        npt_gro = "output/pose_runs/{pose}/npt/npt.gro",
        npt_cpt = "output/pose_runs/{pose}/npt/npt.cpt", 
        ion_top = "output/pose_runs/{pose}/ions/ion.top",
        template_mdp = "mdp/pull_short.mdp",
        active_site_residues = "output/autodock_vina/active_site_residues_updated.txt"
    output:
        index = "output/pose_runs/{pose}/smd/index.ndx",
        tpr = "output/pose_runs/{pose}/smd/smd.tpr",
        smd_gro = "output/pose_runs/{pose}/smd/smd.gro",
        smd_vis_pdb = "output/pose_runs/{pose}/smd/smd_vis.pdb",
        pullf = "output/pose_runs/{pose}/smd/pullf.xvg",
        pullx = "output/pose_runs/{pose}/smd/pullx.xvg",
        success = "output/pose_runs/{pose}/smd/smd.success"
    conda: "envs/environment.yml"
    shell:
        """
        set -e
        mkdir -p output/pose_runs/{wildcards.pose}/smd
        
        # Create index file with ActiveSite group for SMD
        python scripts/add_active_site_group.py \
            --gro {input.npt_gro} \
            --topology {input.ion_top} \
            --residues {input.active_site_residues} \
            --output output/pose_runs/{wildcards.pose}/smd
        
        # Use template MDP directly - no vector calculation needed for distance geometry
        cp {input.template_mdp} output/pose_runs/{wildcards.pose}/smd/pull.mdp
        
        # Run SMD with crash recovery (let it crash when ligand reaches active site)
        python scripts/run_smd_with_recovery.py \
            --mdp output/pose_runs/{wildcards.pose}/smd/pull.mdp \
            --gro {input.npt_gro} \
            --cpt {input.npt_cpt} \
            --top {input.ion_top} \
            --ndx {output.index} \
            --output output/pose_runs/{wildcards.pose}/smd/smd
        
        # Create visualization PDB with proper PBC correction  
        printf "1\n0\n" | gmx trjconv -s output/pose_runs/{wildcards.pose}/smd/smd.tpr -f output/pose_runs/{wildcards.pose}/smd/smd.gro -o output/pose_runs/{wildcards.pose}/smd/smd_vis.pdb -pbc mol -center -ur compact
        
        # Rename pull output files to match expected names
        mv output/pose_runs/{wildcards.pose}/smd/smd_pullf.xvg {output.pullf}
        mv output/pose_runs/{wildcards.pose}/smd/smd_pullx.xvg {output.pullx}
        
        echo OK > {output.success}
        """

rule analyze_pose_forces:
    input:
        pullfs = pose_pullfs_input
    output:
        report = "output/pose_runs/pose_forces_summary.txt"
    conda: "envs/environment.yml"
    shell:
        """
        set -e
        python scripts/analyze_force_spikes.py --inputs {input.pullfs} --output {output.report}
        """

rule run_mole2:
    input:
        production_tpr = "output/md/production_apo.tpr",
        production_gro = "output/md/production_apo.gro",
        ligand_sdf = config["new_ligand_file"],
        mole2_config = "mole2_config.xml",
        active_site_residues = "output/protein/active_site_residues.txt",
        original_residues = "output/autodock_vina/active_site_residues_original.txt",
        vina_success = "output/autodock_vina/vina.success"
    output:
        mole2_xml = "output/mole2_output/mole2_config.xml",
        mole2_log = "output/mole2_output/mole2.log",
        mole2_success = "output/mole2_output/mole2.success",
        tunnels_xml = "output/mole2_output/xml/tunnels.xml"
    conda: "envs/environment.yml"
    shell:
        """
        set -e  # Exit on any error
        
        echo "=== DEBUG: Starting Mole2 workflow ==="
        
        mkdir -p output/mole2_output
        
        echo "=== DEBUG: Creating Mole2 XML configuration ==="
        # Copy the template XML file
        cp {input.mole2_config} {output.mole2_xml}
        
        # Keep using original apo_protein.pdb path (no change needed)
        
        echo "=== DEBUG: Original XML content ==="
        cat {output.mole2_xml}
        
        # Read original residues for MOLE2 (using original PDB coordinate system)
        residues=$(cat {input.original_residues})
        echo "=== DEBUG: Original active site residues for MOLE2: $residues ==="
        
        # Add origins for each pocket residue to the XML
        IFS=',' read -ra RESIDUES <<< "$residues"
        for residue in "${{RESIDUES[@]}}"; do
            echo "=== DEBUG: Adding origin for residue $residue ==="
            # Insert origin before the closing </Origins> tag
            sed -i "s|  </Origins>|    <Origin><Residue SequenceNumber=\"$residue\" Chain=\"A\" /></Origin>\\n  </Origins>|" {output.mole2_xml}
        done
        
        echo "=== DEBUG: Final XML content ==="
        cat {output.mole2_xml}
        
        echo "=== DEBUG: Fixing XML attributes ==="
        # Fix the XML attributes by properly quoting them
        python scripts/fix_mole2_xml.py {output.mole2_xml} {output.mole2_xml}
        
        echo "=== DEBUG: Fixed XML content ==="
        cat {output.mole2_xml}
        
        echo "=== DEBUG: Running Mole2 ==="
        # Run Mole2 with the configuration
        cd Mole2_cmd
        echo "=== DEBUG: Current directory: $(pwd) ==="
        echo "=== DEBUG: Mole2 executable exists: $(ls -la mole2.exe) ==="
        echo "=== DEBUG: XML file exists: $(ls -la ../{output.mole2_xml}) ==="
        
        mono mole2.exe ../{output.mole2_xml} > ../{output.mole2_log} 2>&1
        mole2_exit_code=$?
        echo "=== DEBUG: Mole2 exit code: $mole2_exit_code ==="
        
        echo "=== DEBUG: Mole2 log content ==="
        cat ../{output.mole2_log}
        
        echo "=== DEBUG: Checking Mole2 output files ==="
        ls -la mole2_output/ || echo "mole2_output directory not found"
        ls -la mole2_output/xml/ || echo "xml directory not found"
        
        # Ensure project-level output path contains the XML files
        mkdir -p ../output/mole2_output/xml
        cp -r mole2_output/xml/* ../output/mole2_output/xml/
        
        # Validate success (now check via project-level path)
        if [ ! -f ../output/mole2_output/xml/tunnels.xml ]; then
            echo "=== DEBUG: tunnels.xml not found at project-level ==="
            ls -la ../output/mole2_output/ || true
            ls -la ../output/mole2_output/xml/ || true
            echo "Error: Mole2 failed - output files not created at expected location"
            exit 1
        fi
        
        echo "=== DEBUG: Mole2 completed successfully ==="
        # Create success marker
        echo "Mole2 tunnel detection completed successfully" > ../{output.mole2_success}
        """

# Create centered protein PDB for MOLE2 analysis (same coordinate system as simulation)
rule create_centered_protein_pdb:
    input:
        protein_gro = "output/protein/clean_apo.gro"
    output:
        centered_pdb = "output/protein/clean_apo_centered.pdb"
    conda: "envs/environment.yml"
    shell:
        """
        set -e
        # Center the protein at origin (same as done in complex creation)
        gmx editconf -f {input.protein_gro} -o output/protein/clean_apo_centered.gro -c -center 0 0 0
        # Convert centered GRO to PDB for MOLE2 analysis
        gmx editconf -f output/protein/clean_apo_centered.gro -o {output.centered_pdb}
        echo "Created centered protein PDB for MOLE2 analysis"
        """

# Compute tunnel-based ligand placement using centered coordinates (no transformation needed)
rule compute_tunnel_placement:
    input:
        tunnels_xml = "output/mole2_output/xml/tunnels.xml",
        residues = "output/autodock_vina/active_site_residues_updated.txt",
        original_pdb = "output/protein/apo_protein.pdb",
        ligand_sdf = config["new_ligand_file"],
        mole2_success = "output/mole2_output/mole2.success"
    output:
        positions = "output/ligand/ligand_offset_vector.dat",
        info_txt = "output/ligand/tunnel_info.txt"
    conda: "envs/environment.yml"
    shell:
        """
        set -e
        python scripts/select_best_tunnel.py \
            --tunnels_xml {input.tunnels_xml} \
            --ligand_sdf {input.ligand_sdf} \
            --residues {input.residues} \
            --pdb {input.original_pdb} \
            --protein_pdb {input.original_pdb} \
            --positions_out {output.positions} \
            --tunnel_info_out {output.info_txt} \
            --offset_nm 1.5
        """

rule select_best_pose:
    input:
        pose_analysis = "output/pose_runs/pose_forces_summary.txt"
    output:
        best_pose = directory("output/best_pose/"),
        selection_success = "output/best_pose/selection.success"
    conda: "envs/environment.yml"
    shell:
        """
        set -e
        
        # Determine the best pose from the analysis results
        best_pose=$(grep "Best window:" {input.pose_analysis} | cut -d'/' -f3)
        echo "Selected best pose: $best_pose"
        
        # Create best_pose directory
        mkdir -p output/best_pose
        
        # Copy the equilibrated files from the best pose (ready for SMD)
        cp -r "output/pose_runs/$best_pose/"* {output.best_pose}
        
        echo "Best pose files copied successfully from $best_pose" > {output.selection_success}
        echo "Best pose: $best_pose" >> {output.selection_success}
        """

rule resume_simulation:
    input:
        best_pose_dir = "output/best_pose/",
        selection_success = "output/best_pose/selection.success",
        prod_mdp = "mdp/production.mdp"
    output:
        smd_tpr = "output/smd/smd.tpr",
        smd_xtc = "output/smd/smd.xtc",
        smd_log = "output/smd/smd.log",
        smd_success = "output/smd/smd.success"
    conda: "envs/environment.yml"
    shell:
        """
        set -e
        mkdir -p output/smd
        
        # Create a 100 ns production MDP for free binding simulation
        sed 's/nsteps.*/nsteps = 50000000  ; 100 ns free MD simulation/' {input.prod_mdp} > output/smd/production_100ns.mdp
        
        # Prepare 100 ns free MD simulation (no pulling forces)
        gmx grompp -f output/smd/production_100ns.mdp -c output/best_pose/npt/npt.gro -t output/best_pose/npt/npt.cpt -p output/best_pose/ions/ion.top -o {output.smd_tpr} -maxwarn 1000
        
        # Run 100 ns free MD to observe natural ligand binding
        export OMP_NUM_THREADS=4
        gmx mdrun -v -deffnm output/smd/smd -ntmpi 1 -ntomp 4
        
        # Validate success
        if [ ! -f {output.smd_tpr} ] || [ ! -f {output.smd_xtc} ]; then
            echo "Error: Free MD simulation failed - output files not created"
            exit 1
        fi
        
        echo "100 ns free MD simulation completed successfully" > {output.smd_success}
        """

# Sequential workflow rule to force execution order
rule merged_workflow:
    input:
        selection_success = "output/best_pose/selection.success",
        smd_success = "output/smd/smd.success",
        mmpbsa_success = "output/mmpbsa/mmpbsa.success"
    output:
        workflow_complete = "output/workflow_complete.txt"
    shell:
        """
        echo "Complete workflow (NPT + SMD + MM-PBSA) completed successfully" > {output.workflow_complete}
        """


rule extract_frames:
    input:
        tpr = "output/smd/smd.tpr",
        xtc = "output/smd/smd.xtc",
        smd_success = "output/smd/smd.success"
    output:
        pdb_dir = directory("output/mmpbsa/frames")
    conda: "envs/environment.yml"
    shell:
        """
        set -e
        mkdir -p {output.pdb_dir}
        
        # Extract frames from steered MD trajectory
        echo "0" | gmx trjconv -s {input.tpr} -f {input.xtc} -o {output.pdb_dir}/frame.pdb -sep
        
        echo "Frame extraction completed successfully"
        """


rule mmpbsa_analysis:
    input:
        trajectory = "output/smd/smd.xtc",
        topology = "output/smd/smd.tpr",
        best_pose_dir = "output/best_pose/",
        smd_success = "output/smd/smd.success"
    output:
        results = "output/mmpbsa/binding_analysis.txt",
        plot = "output/mmpbsa/binding_plot.png",
        mmpbsa_success = "output/mmpbsa/mmpbsa.success"
    conda: "envs/environment.yml"
    shell:
        """
        set -e
        mkdir -p output/mmpbsa
        
        # Run MM-PBSA analysis
        python scripts/setup_mmpbsa.py \
            --trajectory {input.trajectory} \
            --topology {input.topology} \
            --system_top output/best_pose/ions/ion.top \
            --output_results {output.results} \
            --output_plot {output.plot}
        
        # Validate success
        if [ ! -f {output.results} ] || [ ! -f {output.plot} ]; then
            echo "Error: MM-PBSA analysis failed - output files not created"
            exit 1
        fi
        
        echo "MM-PBSA analysis completed successfully" > {output.mmpbsa_success}
        """






