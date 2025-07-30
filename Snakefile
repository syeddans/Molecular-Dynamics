configfile: "config.yaml"

# Use conda environments automatically


rule all:
    input:
        "md/production_apo.xtc",
        "md/production_apo.tpr",
        "analysis/rmsd_apo.xvg",
        "analysis/rmsf_apo.xvg",
        # Include ligand preparation
        "ligand/ligand.gro",
        # Include analysis plots
        "analysis/rmsd_apo_plot.png",
        "analysis/rmsf_apo_plot.png"

# Extract original ligand from holo structure (for reference only)
rule extract_original_ligand:
    input:
        holo_pdb = config["holo_protein_pdb"]
    output:
        original_ligand_pdb = "ligand_original/original_ligand.pdb"
    params:
        ligand_resn = config.get("original_ligand_resn", "EST")
    conda: "envs/environment.yml"
    shell:
        """
        mkdir -p ligand_original
        python scripts/extract_ligand.py --input {input.holo_pdb} --output {output.original_ligand_pdb} --ligand_resn {params.ligand_resn}
        """

# Clean protein by removing ligand and using PDBFixer for robust cleaning
rule create_apo_protein:
    input:
        holo_pdb = config["holo_protein_pdb"]
    output:
        apo_pdb = "protein/apo_protein.pdb"
    params:
        ligand_resn = config.get("original_ligand_resn", "EST")
        # target_chain parameter removed - now keeping all protein chains
    conda: "envs/environment.yml"
    shell:
        """
        python scripts/create_apo_protein_clean.py --input {input.holo_pdb} --output {output.apo_pdb} --ligand_resn {params.ligand_resn}
        """

# Process ligand for docking
rule process_ligand:
    input:
        ligand = "inputs/ICI.sdf"
    output:
        ligand_mol2 = "ligand/ligand.mol2",
        gaff_mol2 = "ligand/ligand_gaff.mol2",
        frcmod = "ligand/ligand.frcmod"
    params:
        charge = config.get("ligand_charge", "auto")
    conda: "envs/environment.yml"
    shell:
        """
        mkdir -p ligand
        cd ligand
        
        # Determine input format from file extension
        LIGAND_EXT=$(echo "{input.ligand}" | sed 's/.*\./\./' | tr '[:upper:]' '[:lower:]')
        
        if [ "$LIGAND_EXT" = ".mol2" ]; then
            INPUT_FORMAT="mol2"
        elif [ "$LIGAND_EXT" = ".sdf" ]; then
            INPUT_FORMAT="mdl"
        else
            INPUT_FORMAT="pdb"  # Default for .pdb or others
        fi
        
        # Copy input file to working directory
        cp ../{input.ligand} input_ligand$LIGAND_EXT
        
        # Determine net charge
        if [ "{params.charge}" = "auto" ]; then
            NET_CHARGE="0"  # Default to 0, antechamber will adjust if needed
            echo "No charge specified - using default charge: $NET_CHARGE"
        else
            NET_CHARGE="{params.charge}"
            echo "Using specified charge: $NET_CHARGE"
        fi
        
        # Convert to MOL2 using Gasteiger charges (faster, more reliable than AM1-BCC)
        if [ "$INPUT_FORMAT" != "mol2" ]; then
            echo "Converting to MOL2 and calculating Gasteiger charges..."
            antechamber -i input_ligand$LIGAND_EXT -fi $INPUT_FORMAT -o ligand.mol2 -fo mol2 -c gas -s 2 -nc $NET_CHARGE
        else
            echo "Input is already MOL2, recalculating charges with Gasteiger method..."
            antechamber -i input_ligand$LIGAND_EXT -fi mol2 -o ligand.mol2 -fo mol2 -c gas -s 2 -nc $NET_CHARGE
        fi
        
        # Generate GAFF parameters (use simple method)
        echo "Generating GAFF parameters..."
        antechamber -i ligand.mol2 -fi mol2 -o ligand_gaff.mol2 -fo mol2 -c gas -s 2
        
        # Generate additional parameters
        echo "Generating additional force field parameters..."
        parmchk2 -i ligand_gaff.mol2 -f mol2 -o ligand.frcmod
        
        echo "Ligand parameterization complete!"
        echo "Net molecular charge used: $NET_CHARGE"
        echo "Atomic charges calculated using Gasteiger method"
        """

# Convert ligand to GROMACS format
rule ligand_to_gromacs:
    input:
        gaff_mol2 = "ligand/ligand_gaff.mol2",
        frcmod = "ligand/ligand.frcmod"
    output:
        ligand_gro = "ligand/ligand.gro",
        ligand_itp = "ligand/ligand.itp",
        ligand_top = "ligand/ligand.top"
    conda: "envs/environment.yml"
    shell:
        """
        cd ligand
        
        # Create AMBER topology files first
        antechamber -i ligand_gaff.mol2 -fi mol2 -o ligand.prep -fo prepi
        antechamber -i ligand_gaff.mol2 -fi mol2 -o ligand.ac -fo ac
        
        # Use tleap to create AMBER topology
        cat > leap.in << 'EOF'
source leaprc.gaff
loadAmberParams ligand.frcmod
loadAmberPrep ligand.prep
mol = loadmol2 ligand_gaff.mol2
saveAmberParm mol ligand.prmtop ligand.rst7
quit
EOF
        
        tleap -f leap.in
        
        # Convert to GROMACS using parmed
        cat > convert_to_gromacs.py << 'EOF'
import parmed as pmd
parm = pmd.load_file('ligand.prmtop', 'ligand.rst7')
parm.save('ligand.gro', format='gro')
parm.save('ligand.top', format='gromacs')
parm.save('ligand.itp', format='gromacs')
EOF
        
        python convert_to_gromacs.py
        
        echo "Ligand conversion to GROMACS format complete!"
        """

# Process cleaned apo protein with pdb2gmx (handles both single and multi-chain)
rule process_apo_protein:
    input:
        apo_pdb = "protein/apo_protein.pdb"
    output:
        protein_gro = "protein/protein.gro",
        protein_top = "protein/topol.top",
        posre_itp = "protein/posre.itp"
    params:
        ff = config.get("force_field", "amber99sb-ildn"),
        water = config.get("water_model", "tip3p")
    conda: "envs/environment.yml"
    shell:
        """
        # Run GROMACS pdb2gmx
        cd protein
        gmx pdb2gmx -f apo_protein.pdb -o protein.gro -p topol.top -ff {params.ff} -water {params.water} -ignh
        
        # Handle position restraint files (single vs multi-chain)
        python ../scripts/handle_posre_files.py --protein_dir .
        """

# Create simulation box and solvate apo protein
rule solvate_apo:
    input:
        protein_gro = "protein/protein.gro",
        protein_top = "protein/topol.top"
    output:
        boxed_gro = "solvation/apo_boxed.gro",
        solvated_gro = "solvation/apo_solvated.gro",
        solvated_top = "solvation/topol.top"
    params:
        box_distance = config.get("box_distance", "1.2")
    conda: "envs/environment.yml"
    shell:
        """
        mkdir -p solvation
        
        # Start with a clean topology (remove any existing water)
        sed '/^SOL/d' {input.protein_top} > solvation/clean_topol.top
        
        # Create box using rhombic dodecahedron for better performance (15% fewer water molecules)
        gmx editconf -f {input.protein_gro} -o {output.boxed_gro} -c -d {params.box_distance} -bt dodecahedron
        
        # Solvate using clean topology
        gmx solvate -cp {output.boxed_gro} -cs spc216.gro -o {output.solvated_gro} -p solvation/clean_topol.top
        
        # Copy final topology and all required include files
        cp solvation/clean_topol.top {output.solvated_top}
        cp protein/*.itp solvation/ 2>/dev/null || true
        """

# Add ions for neutralization 
rule add_ions_apo:
    input:
        solvated_gro = "solvation/apo_solvated.gro",
        solvated_top = "solvation/topol.top",
        ions_mdp = "mdp/ions.mdp"
    output:
        ions_tpr = "solvation/ions.tpr",
        system_gro = "solvation/apo_system.gro",
        system_top = "solvation/system_topol.top"
    conda: "envs/environment.yml"
    shell:
        """
        # Create tpr for genion
        gmx grompp -f {input.ions_mdp} -c {input.solvated_gro} -p {input.solvated_top} -o {output.ions_tpr} -maxwarn 2
        
        # Add ions (select SOL when prompted)
        echo "SOL" | gmx genion -s {output.ions_tpr} -o {output.system_gro} -p {input.solvated_top} -pname NA -nname CL -neutral -conc 0.15
        
        # Copy final topology
        cp {input.solvated_top} {output.system_top}
        """

# Energy minimization
rule energy_minimization:
    input:
        system_gro = "solvation/apo_system.gro",
        system_top = "solvation/system_topol.top",
        em_mdp = "mdp/em.mdp"
    output:
        em_tpr = "em/em.tpr",
        em_gro = "em/em.gro",
        em_edr = "em/em.edr"
    threads: config.get("eq_threads", 8)
    params:
        ntomp = config.get("openmp_threads", 4)
    conda: "envs/environment.yml"
    shell:
        """
        mkdir -p em
        
        # Prepare EM
        gmx grompp -f {input.em_mdp} -c {input.system_gro} -p {input.system_top} -o {output.em_tpr} -maxwarn 1
        
        # Run EM with threading optimization
        # Calculate MPI ranks for energy minimization
        NTMPI=$((({threads} + {params.ntomp} - 1) / {params.ntomp}))
        export OMP_NUM_THREADS={params.ntomp}
        gmx mdrun -v -deffnm em/em -nt {threads} -ntomp {params.ntomp} -ntmpi $NTMPI -pin on
        """

# NVT equilibration (temperature coupling)
rule nvt_equilibration:
    input:
        em_gro = "em/em.gro",
        system_top = "solvation/system_topol.top",
        nvt_mdp = "mdp/nvt.mdp"
    output:
        nvt_tpr = "nvt/nvt.tpr",
        nvt_gro = "nvt/nvt.gro",
        nvt_edr = "nvt/nvt.edr",
        nvt_cpt = "nvt/nvt.cpt"
    threads: config.get("eq_threads", 8)
    params:
        ntomp = config.get("openmp_threads", 4)
    conda: "envs/environment.yml"
    shell:
        """
        mkdir -p nvt
        
        # Prepare NVT
        gmx grompp -f {input.nvt_mdp} -c {input.em_gro} -r {input.em_gro} -p {input.system_top} -o {output.nvt_tpr} -maxwarn 1
        
        # Run NVT (100 ps) with CPU performance optimization
        # Calculate MPI ranks for equilibration
        NTMPI=$((({threads} + {params.ntomp} - 1) / {params.ntomp}))
        export OMP_NUM_THREADS={params.ntomp}
        gmx mdrun -v -deffnm nvt/nvt -nt {threads} -ntomp {params.ntomp} -ntmpi $NTMPI -pin on
        """

# NPT equilibration (pressure coupling)
rule npt_equilibration:
    input:
        nvt_gro = "nvt/nvt.gro",
        nvt_cpt = "nvt/nvt.cpt",
        system_top = "solvation/system_topol.top",
        npt_mdp = "mdp/npt.mdp"
    output:
        npt_tpr = "npt/npt.tpr",
        npt_gro = "npt/npt.gro",
        npt_edr = "npt/npt.edr",
        npt_cpt = "npt/npt.cpt"
    threads: config.get("eq_threads", 8)
    params:
        ntomp = config.get("openmp_threads", 4)
    conda: "envs/environment.yml"
    shell:
        """
        mkdir -p npt
        
        # Prepare NPT  
        gmx grompp -f {input.npt_mdp} -c {input.nvt_gro} -t {input.nvt_cpt} -r {input.nvt_gro} -p {input.system_top} -o {output.npt_tpr} -maxwarn 1
        
        # Run NPT (100 ps) with CPU performance optimization
        # Calculate MPI ranks for equilibration  
        NTMPI=$((({threads} + {params.ntomp} - 1) / {params.ntomp}))
        export OMP_NUM_THREADS={params.ntomp}
        gmx mdrun -v -deffnm npt/npt -nt {threads} -ntomp {params.ntomp} -ntmpi $NTMPI -pin on
        """

# Production MD simulation (long run for apo relaxation)
rule production_md_apo:
    input:
        npt_gro = "npt/npt.gro",
        npt_cpt = "npt/npt.cpt", 
        system_top = "solvation/system_topol.top",
        prod_mdp = "mdp/production.mdp"
    output:
        prod_tpr = "md/production_apo.tpr",
        prod_xtc = "md/production_apo.xtc",
        prod_gro = "md/production_apo.gro",
        prod_edr = "md/production_apo.edr",
        prod_log = "md/production_apo.log",
        prod_cpt = "md/production_apo.cpt"
    threads: config.get("md_threads", 16)
    params:
        pme_ranks = config.get("pme_ranks", 4),
        ntomp = config.get("openmp_threads", 4)
    conda: "envs/environment.yml"
    shell:
        """
        mkdir -p md
        
        # Prepare production MD
        gmx grompp -f {input.prod_mdp} -c {input.npt_gro} -t {input.npt_cpt} -p {input.system_top} -o {output.prod_tpr} -maxwarn 1
        
        # Run CPU-optimized production MD with PME load balancing
        # Calculate MPI ranks: total_threads / openmp_threads_per_rank
        NTMPI=$((({threads} + {params.ntomp} - 1) / {params.ntomp}))
        
        if [ "{params.pme_ranks}" != "auto" ]; then
            # Ensure PME ranks < total MPI ranks
            if [ {params.pme_ranks} -ge $NTMPI ]; then
                PME_RANKS=$((NTMPI - 1))
            else
                PME_RANKS={params.pme_ranks}
            fi
            PME_ARGS="-ntmpi $NTMPI -npme $PME_RANKS"
        else
            PME_ARGS="-ntmpi $NTMPI"
        fi
        
        export OMP_NUM_THREADS={params.ntomp}
        gmx mdrun -v -deffnm md/production_apo -nt {threads} -ntomp {params.ntomp} -pin on $PME_ARGS
        """

# Basic analysis - RMSD
rule analyze_rmsd:
    input:
        prod_xtc = "md/production_apo.xtc",
        prod_tpr = "md/production_apo.tpr"
    output:
        rmsd_xvg = "analysis/rmsd_apo.xvg"
    conda: "envs/environment.yml"
    shell:
        """
        mkdir -p analysis
        
        # Calculate RMSD (Backbone selection: 4)
        echo "4 4" | gmx rms -s {input.prod_tpr} -f {input.prod_xtc} -o {output.rmsd_xvg} -tu ns
        """

# Basic analysis - RMSF  
rule analyze_rmsf:
    input:
        prod_xtc = "md/production_apo.xtc",
        prod_tpr = "md/production_apo.tpr"
    output:
        rmsf_xvg = "analysis/rmsf_apo.xvg"
    conda: "envs/environment.yml"
    shell:
        """
        mkdir -p analysis
        
        # Calculate RMSF (Backbone selection: 4)
        echo "4" | gmx rmsf -s {input.prod_tpr} -f {input.prod_xtc} -o {output.rmsf_xvg} -res
        """

# Generate RMSD plot
rule plot_rmsd:
    input:
        rmsd_xvg = "analysis/rmsd_apo.xvg"
    output:
        rmsd_plot = "analysis/rmsd_apo_plot.png"
    conda: "envs/environment.yml"
    shell:
        """
        python plot_rmsd.py
        """

# Generate RMSF plot
rule plot_rmsf:
    input:
        rmsf_xvg = "analysis/rmsf_apo.xvg"
    output:
        rmsf_plot = "analysis/rmsf_apo_plot.png"
    conda: "envs/environment.yml"
    shell:
        """
        python scripts/plot_rmsf.py --input {input.rmsf_xvg} --output {output.rmsf_plot}
        """

# Extract representative apo conformations for docking
rule extract_apo_conformations:
    input:
        prod_xtc = "md/production_apo.xtc",
        prod_tpr = "md/production_apo.tpr"
    output:
        apo_conf1 = "docking_prep/apo_conf1.pdb",
        apo_conf2 = "docking_prep/apo_conf2.pdb",
        apo_conf3 = "docking_prep/apo_conf3.pdb"
    conda: "envs/environment.yml"
    shell:
        """
        mkdir -p docking_prep
        
        # Extract conformations at different time points for ensemble docking
        # Frame at 250 ns
        echo "1" | gmx trjconv -s {input.prod_tpr} -f {input.prod_xtc} -o {output.apo_conf1} -dump 250000
        
        # Frame at 500 ns  
        echo "1" | gmx trjconv -s {input.prod_tpr} -f {input.prod_xtc} -o {output.apo_conf2} -dump 500000
        
        # Frame at 750 ns
        echo "1" | gmx trjconv -s {input.prod_tpr} -f {input.prod_xtc} -o {output.apo_conf3} -dump 750000
        """
