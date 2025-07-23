# Testing Guide for gRINN

## Overview
This guide helps you systematically test gRINN and report useful feedback to improve the tool for everyone.

gRINN (protein Residue Interaction Networks from Nanosecond molecular dynamics) is a comprehensive workflow that combines GROMACS molecular dynamics simulations with sophisticated energy decomposition analysis. It can both run complete MD simulations from PDB structures AND analyze existing trajectory data to construct protein energy networks.

## gRINN Capabilities

### 🧬 Molecular Dynamics Preparation
- **PDB Processing**: Automatic fixing of PDB files using PDBFixer
- **Topology Generation**: GROMACS topology creation with multiple force fields
- **System Preparation**: Solvation, ionization, and energy minimization
- **Equilibration**: NVT and NPT equilibration protocols

### ⚗️ Simulation Options
- **Vacuum Simulations**: Fast analysis for small proteins
- **Solvated Systems**: Full explicit water MD simulations
- **GPU Acceleration**: CUDA support for faster calculations
- **Custom Parameters**: Flexible MD parameters and force field selection

### 📊 Energy Analysis
- **Interaction Decomposition**: Per-residue pairwise energy analysis
- **Energy Types**: Van der Waals and electrostatic contributions
- **Frame-by-frame Analysis**: Temporal evolution of interactions
- **Network Construction**: Protein Energy Networks (PENs) with centrality metrics

### 🎯 Advanced Features
- **Residue Selection**: Focus on specific protein regions or binding sites
- **Multiple Cutoffs**: Network construction with various energy thresholds
- **Betweenness Centrality**: Identification of allosteric communication pathways
- **Interactive Dashboard**: 3D visualization with trajectory playback

## Quick Testing Checklist

### ✅ Basic Functionality Test
- [ ] Docker image builds successfully (15-20 minutes)
- [ ] Built-in test data runs without errors
- [ ] Simple PDB simulation completes (vacuum MD)
- [ ] Dashboard launches and displays results
- [ ] Energy matrices look reasonable
- [ ] Network plots are generated

### ✅ Performance Test
- [ ] Note runtime for small protein (< 50 residues)
- [ ] Note runtime for medium protein (50-150 residues)
- [ ] Test vacuum vs solvated simulation times
- [ ] Monitor memory usage during calculation
- [ ] Test with/without GPU acceleration

### ✅ Advanced Features
- [ ] Solvated MD simulations with NPT equilibration
- [ ] Custom trajectory analysis with existing data
- [ ] Energy network construction (--create_pen)
- [ ] Residue selection analysis (--source_sel, --target_sel)
- [ ] Multiple energy cutoffs and centrality metrics

### ✅ Edge Cases
- [ ] Multi-chain proteins
- [ ] Proteins with missing residues
- [ ] Very small proteins (< 20 residues)
- [ ] Proteins with unusual amino acids
- [ ] Empty or corrupted input files
- [ ] Large proteins (>300 residues)

## Systematic Testing Protocol

### Step 1: Environment Setup
```bash
# Test Docker installation
docker --version
docker info

# Build gRINN image (this takes ~15-20 minutes due to GROMACS compilation)
git clone https://github.com/osercinoglu/grinn.git
cd grinn
docker build -t grinn .

# Verify build and check included test data
docker images | grep grinn
docker run grinn help

# Test that built-in data is available
docker run grinn bash -c "ls -la prot_lig_1/"
# Should show: npt_dry_wCh.pdb, topol_dry.top, md_0_1000ns_merged_dry_nojump_skip1000.xtc, etc.
```

### Step 2: Functional Testing

#### A. Built-in Trajectory Analysis (Fastest test)
```bash
# Test with built-in trajectory data (Docker includes test data from Google Cloud)
# This analyzes a pre-computed protein-ligand MD trajectory
mkdir test_results

# Run trajectory analysis on built-in test data
# This will decompose interaction energies from existing trajectory
docker run -v $(pwd)/test_results:/app/test_results grinn \
  workflow prot_lig_1/npt_dry_wCh.pdb test_prot_lig_1 \
  --top prot_lig_1/topol_dry.top \
  --traj prot_lig_1/md_0_1000ns_merged_dry_nojump_skip1000.xtc \
  --nt 4 \
  --ff_folder prot_lig_1/charmm36.ff

# Verify trajectory analysis outputs
ls test_results/
# Should contain: energies_intEnTotal.csv, energies_intEnVdW.csv, energies_intEnElec.csv, 
#                system_dry.pdb, traj_superposed.xtc, calc.log, gromacs.log

# Test dashboard with trajectory analysis results
docker run -p 8051:8051 -v $(pwd)/test_results:/results grinn dashboard /results
# Check: http://localhost:8051
```

#### B. Complete MD Simulation from PDB (Full workflow test)
```bash
# Test complete gRINN workflow starting from PDB file
# This tests GROMACS preparation, simulation, and analysis
wget https://files.rcsb.org/download/1L2Y.pdb
mkdir results_full_workflow

# Option 1: Quick vacuum simulation (~5-10 minutes)
docker run -v $(pwd):/data -v $(pwd)/results_full_workflow:/results grinn \
  workflow /data/1L2Y.pdb /results --nt 4

# Option 2: Full solvated simulation with equilibration (~30-60 minutes)
docker run -v $(pwd):/data -v $(pwd)/results_full_workflow:/results grinn \
  workflow /data/1L2Y.pdb /results --solvate --npt --nt 4 --gpu

# Verify complete workflow outputs
ls results_full_workflow/
# Should contain: All trajectory analysis files PLUS
#                gromacs_files/, system_prepared.pdb, topol.top, etc.
```

### Step 3: Performance Benchmarking

#### A. Simulation Performance Testing
```bash
# Small protein complete workflow (Trp-cage)
wget https://files.rcsb.org/download/1L2Y.pdb
mkdir results_small

# Test 1: Vacuum MD (fast baseline, ~5-10 minutes)
time docker run -v $(pwd):/data -v $(pwd)/results_small:/results grinn \
  workflow /data/1L2Y.pdb /results --nt 4

# Test 2: Solvated MD with equilibration (~30-60 minutes)
time docker run -v $(pwd):/data -v $(pwd)/results_small:/results grinn \
  workflow /data/1L2Y.pdb /results --solvate --npt --nt 4 --gpu

# Medium protein workflow (Lysozyme)
wget https://files.rcsb.org/download/1LYS.pdb
mkdir results_medium

# Full workflow with energy networks
time docker run -v $(pwd):/data -v $(pwd)/results_medium:/results grinn \
  workflow /data/1LYS.pdb /results --solvate --npt --nt 4 --gpu \
  --create_pen --pen_cutoffs 0.5 1.0 2.0 --pen_include_covalents True False
```

#### B. Trajectory Analysis Performance
```bash
# Test trajectory analysis with different sampling rates
mkdir results_sampling

# Test 1: Analyze every frame (slower, more complete)
time docker run -v $(pwd)/test_results:/input -v $(pwd)/results_sampling:/results grinn \
  workflow /input/system_dry.pdb /results \
  --top prot_lig_1/topol_dry.top \
  --traj /input/traj_superposed.xtc \
  --skip 1 --nt 4

# Test 2: Analyze every 10th frame (faster)
time docker run -v $(pwd)/test_results:/input -v $(pwd)/results_sampling:/results grinn \
  workflow /input/system_dry.pdb /results \
  --top prot_lig_1/topol_dry.top \
  --traj /input/traj_superposed.xtc \
  --skip 10 --nt 4

# Test 3: Your own trajectory data
# Prepare your GROMACS trajectory files in a folder
mkdir my_trajectory_test
# Place your files: protein.pdb, topol.top, trajectory.xtc

time docker run -v $(pwd)/my_trajectory_test:/data -v $(pwd)/results_custom:/results grinn \
  workflow /data/protein.pdb /results \
  --top /data/topol.top \
  --traj /data/trajectory.xtc \
  --ff_folder /data/forcefield.ff \
  --nt 8

# Note: Monitor memory usage with: docker stats
```

### Step 4: Validation Testing

#### A. Energy Analysis Validation
```bash
# Check energy decomposition results (should be chemically reasonable)
head -5 results_small/energies_intEnTotal.csv

# Example output format:
# ResidueA,ResidueB,Energy_kcal_mol
# ALA_1,VAL_2,-0.45
# ALA_1,PRO_3,0.12
# VAL_2,PRO_3,-1.23

# Strong interactions should be < -1.0 kcal/mol
# Most interactions should be in range -5.0 to +2.0 kcal/mol

# Compare frame-by-frame vs averaged energies
wc -l results_small/energies_*.csv
# intEnTotal should have all pairwise interactions
# intEnVdW and intEnElec should sum to intEnTotal
```

#### B. GROMACS Workflow Validation
```bash
# Check that GROMACS preparation was successful
ls results_small/gromacs_files/
# Should contain: system_prepared.pdb, topol.top, minimization logs, etc.

# Verify topology file integrity
docker run -v $(pwd)/results_small:/data grinn \
  gmx check -s /data/gromacs_files/topol.tpr

# Check simulation completion
grep "Finished mdrun" results_small/gromacs.log
# Should show successful completion without fatal errors
```

#### C. Scientific Validation
```bash
# Verify network makes sense biologically
# Strong interactions should be chemically reasonable:
# - Salt bridges (charged residues: ARG, LYS, ASP, GLU)
# - Hydrogen bonds (polar residues: SER, THR, TYR, ASN, GLN)  
# - Hydrophobic contacts (nonpolar residues: ALA, VAL, LEU, ILE, PHE)

# If you used --create_pen, check network metrics
ls results_medium/pen_*
# Should contain: pen_cutoff_0.5_cov_True.csv, pen_centralities_*.csv

# Verify trajectory processing was successful
ls -la results_small/
# Should contain processed trajectory: traj_superposed.xtc
# Check frame count matches expected: skip=1 → all frames, skip=10 → every 10th frame
```

### Step 5: Advanced Testing Options

#### A. Complete Workflow Features
```bash
# Test PDB fixing and preparation
docker run -v $(pwd):/data -v $(pwd)/results:/results grinn \
  workflow /data/protein.pdb /results \
  --nofixpdb \
  --solvate --npt --gpu

# Test custom force fields
docker run -v $(pwd):/data -v $(pwd)/results:/results grinn \
  workflow /data/protein.pdb /results \
  --ff_folder /data/amber99sb.ff \
  --solvate --nt 8

# Test with custom MD parameters
# Place your .mdp files in mdp_custom/
docker run -v $(pwd):/data -v $(pwd)/results:/results grinn \
  workflow /data/protein.pdb /results \
  --solvate --npt \
  --include_files /data/mdp_custom/npt_custom.mdp
```

#### B. Advanced Trajectory Analysis
```bash
# Test energy networks from trajectory analysis
docker run -v $(pwd)/test_results:/input -v $(pwd)/results_networks:/results grinn \
  workflow /input/system_dry.pdb /results \
  --top prot_lig_1/topol_dry.top \
  --traj /input/traj_superposed.xtc \
  --create_pen --pen_cutoffs 0.5 1.0 2.0 \
  --pen_include_covalents True False \
  --nt 8

# Test specific residue pair analysis
docker run -v $(pwd)/my_trajectory_test:/data -v $(pwd)/results_binding:/results grinn \
  workflow /data/protein.pdb /results \
  --top /data/topol.top \
  --traj /data/trajectory.xtc \
  --source_sel "resid 50 to 60" \
  --target_sel "resname LIG" \
  --nt 4

# Test large trajectory with frame skipping
docker run -v $(pwd)/my_trajectory_test:/data -v $(pwd)/results_large:/results grinn \
  workflow /data/protein.pdb /results \
  --top /data/topol.top \
  --traj /data/long_trajectory.xtc \
  --skip 50 \
  --nt 16
```

#### C. GPU and Performance Testing
```bash
# Test GPU acceleration (if available)
docker run --gpus all -v $(pwd):/data -v $(pwd)/results:/results grinn \
  workflow /data/protein.pdb /results --gpu --solvate --npt

# Test only validation (check system compatibility)
docker run -v $(pwd):/data grinn \
  workflow /data/protein.pdb /tmp/test --test-only

# Test multi-threading performance
for nt in 1 2 4 8; do
  echo "Testing with $nt threads..."
  time docker run -v $(pwd):/data -v $(pwd)/results_nt$nt:/results grinn \
    workflow /data/1L2Y.pdb /results --nt $nt
done
```

## What to Look For

### ✅ Success Indicators
- **Build time**: 15-20 minutes (includes GROMACS compilation)
- **MD simulation**: Completes without GROMACS fatal errors
- **PDB preparation**: Successful topology generation and system setup
- **Trajectory processing**: Successful frame loading and superposition
- **Energy decomposition**: Reasonable interaction energies (-5.0 to +2.0 kcal/mol)
- **Strong interactions**: < -1.0 kcal/mol, chemically reasonable pairs
- **GROMACS integration**: Proper .top, .tpr, .mdp file handling
- **Network topology**: If using --create_pen, betweenness centrality files created
- **Dashboard functionality**: Smooth 3D trajectory visualization, interactive energy matrices
- **File outputs**: CSV energy files, processed trajectories, GROMACS files, log files

### 🚨 Warning Signs
- **GROMACS build failures**: Fatal errors during topology generation or minimization
- **PDB processing errors**: Cannot fix or process input structure files
- **Trajectory loading fails**: Cannot read .xtc/.trr files, topology mismatch errors
- **Extreme energies**: > +10 or < -20 kcal/mol (possible force field issues)
- **Empty outputs**: Missing CSV files, zero-size trajectory files
- **Frame count mismatch**: Expected vs actual frames processed
- **Simulation crashes**: GROMACS mdrun fatal errors during MD steps
- **Dashboard errors**: Cannot load trajectory, blank visualizations
- **Memory issues**: Docker container killed (exit code 137) = out of memory
- **Permission errors**: Cannot write to mounted volumes

### 🐛 Common Issues & Solutions
- **Docker build fails**: 
  - Check Docker version ≥ 20.10
  - Ensure sufficient disk space (≥10GB free)
  - Try `docker system prune` to free space
- **Permission errors with volumes**: 
  - Use `chmod 777 results_folder` before running
  - On SELinux systems: add `:Z` flag `docker run -v $(pwd)/results:/results:Z`
- **GROMACS topology errors**:
  - Check PDB file format and completeness
  - Try with `--nofixpdb` to disable automatic PDB fixing
  - Verify force field compatibility with your system
- **Trajectory/topology mismatch**:
  - Ensure PDB structure matches trajectory atoms
  - Check that topology (.top) file corresponds to trajectory system
  - Verify force field files are included with --ff_folder
- **MD simulation failures**:
  - Check system preparation: solvation, minimization logs
  - Verify .mdp parameters are appropriate for your system
  - Monitor memory usage: large solvated systems need more RAM
- **Dashboard won't load trajectory**: 
  - Check that traj_superposed.xtc was generated successfully
  - Try different port: `docker run -p 8052:8051`
  - Check browser console for JavaScript errors
- **Out of memory during analysis**:
  - Use `--skip N` to analyze fewer frames (start with --skip 10)
  - Reduce trajectory length or protein size
  - Monitor with `docker stats <container_id>`
- **Force field issues**:
  - Ensure force field directory structure is correct
  - Check for missing .itp files or parameters
  - Try different force field: amber99sb, charmm36, gromos54a7

## Reporting Results

### Success Report Template
```markdown
## gRINN Test Report: [Protein Name]

**System**: [OS, Docker version, CPU cores, RAM]
**Protein**: [PDB ID, size (residues), description]
**Workflow**: [complete MD workflow/trajectory analysis only, energy networks, frame sampling]
**Runtime**: [Build: X min, MD Simulation: X min, Analysis: X min, Total: X min]
**Memory**: [Peak memory usage if monitored]

### Results Summary
- **System preparation**: [vacuum/solvated, force field used]
- **Simulation details**: [ns simulated, frames generated] OR [frames analyzed from existing trajectory]
- **Energy range**: [min to max kcal/mol]
- **Strong interactions**: [count of interactions < -1.0 kcal/mol]
- **Network metrics**: [if --create_pen used: centrality stats]
- **Notable interactions**: [interesting chemical findings from analysis]

### Performance Metrics
- **Docker build**: [X minutes, any issues?]
- **GROMACS preparation**: [X minutes, topology generation success]
- **MD simulation**: [X minutes for Y ns] OR [Trajectory loading: X minutes]
- **Energy decomposition**: [X minutes, memory usage]
- **Dashboard loading**: [trajectory visualization responsive/slow/errors]

### File Outputs ✅
- [x] energies_intEnTotal.csv ([X interactions])
- [x] energies_intEnVdW.csv 
- [x] energies_intEnElec.csv
- [x] system_dry.pdb
- [x] traj_superposed.xtc ([X frames])
- [x] calc.log (no fatal errors)
- [x] gromacs.log (successful completion)

### Issues Found
- [None / List specific problems with details]

### Recommendation
- [Ready for publication / Needs fixes / Requires more testing]
```

### Bug Report Template
```markdown
## gRINN Bug Report: [Brief Description]

**Environment**:
- OS: [Linux/Windows/macOS + version]
- Docker: [version from `docker --version`]
- Hardware: [CPU cores, RAM, GPU if used]
- gRINN version: [git commit hash if known]

**Command Used**:
```bash
# Exact docker run command that failed
docker run -v $(pwd):/data -v $(pwd)/results:/results grinn \
  workflow /data/protein.pdb /results --your --flags --here
```

**Input Files**:
- PDB structure: [name, size, source - attach if < 1MB]
- Topology: [.top file details, force field used] OR [auto-generated from PDB]
- Trajectory: [.xtc/.trr file details, frames, time range] OR [generated by gRINN MD]
- Force field: [charmm36, amber99sb, gromos54a7, etc.]
- MD parameters: [custom .mdp files if used]

**Expected Behavior**:
[What should have happened - complete MD workflow, trajectory loaded, energy files generated, dashboard shows frames, etc.]

**Actual Behavior**:
[What actually happened - GROMACS crashed, trajectory loading failed, missing energy decomposition, etc.]

**Error Messages**:
```
[Full error output from terminal]
[Include relevant lines from calc.log and gromacs.log]
```

**Reproduction Steps**:
1. [Exact sequence of commands to reproduce]
2. [Include file preparations if needed]
3. [When exactly the error occurs]

**Workarounds Tried**:
- [ ] Used --test-only flag for system validation
- [ ] Tried different skip rate (--skip 10, --skip 50)
- [ ] Used --nofixpdb flag to disable PDB fixing
- [ ] Tested with vacuum simulation instead of solvated
- [ ] Checked trajectory/topology compatibility
- [ ] Verified file permissions and disk space
- [ ] Tested with different force field folder
- [ ] Tried without GPU acceleration (removed --gpu)
- [ ] Checked available memory and CPU resources

**Additional Context**:
- [MD simulation details: system size, simulation time, equilibration protocol]
- [Trajectory source: simulation software used, system composition]
- [Screenshots for dashboard trajectory issues]
- [Memory/CPU usage when error occurred]
```

## Advanced Testing

### Complete MD Workflow Testing
```bash
# Test full protein preparation pipeline
docker run -v $(pwd):/data -v $(pwd)/results:/results grinn \
  workflow /data/protein.pdb /results \
  --solvate --npt \
  --ff_folder /data/amber99sb.ff \
  --create_pen --pen_cutoffs 0.5 1.0 2.0

# Test with custom MD parameters
# Create custom_mdp/ folder with your .mdp files
docker run -v $(pwd):/data -v $(pwd)/results:/results grinn \
  workflow /data/protein.pdb /results \
  --solvate --npt \
  --include_files /data/custom_mdp/npt_custom.mdp /data/custom_mdp/nvt_custom.mdp

# Test multi-chain proteins and complexes
docker run -v $(pwd):/data -v $(pwd)/results:/results grinn \
  workflow /data/protein_complex.pdb /results \
  --source_sel "chain A" \
  --target_sel "chain B" \
  --solvate --nt 8
```

### Custom MD Trajectory Analysis
```bash
# If you have existing GROMACS files from other simulations
docker run -v $(pwd):/data -v $(pwd)/results:/results grinn \
  workflow /data/your_protein.pdb /results \
  --top /data/topol.top \
  --traj /data/production.xtc \
  --ff_folder /data/charmm36.ff \
  --skip 5 \
  --nt 8

# Include additional parameter files if needed
docker run -v $(pwd):/data -v $(pwd)/results:/results grinn \
  workflow /data/complex.pdb /results \
  --top /data/topol.top \
  --traj /data/trajectory.xtc \
  --include_files /data/ligand.itp /data/posre.itp /data/ligand.prm

# Test with AMBER or other trajectory formats
docker run -v $(pwd):/data -v $(pwd)/results:/results grinn \
  workflow /data/amber_system.pdb /results \
  --top /data/system.prmtop \
  --traj /data/production.dcd
```

### Stress Testing
```bash
# Large protein test (>500 residues) - monitor memory carefully
docker run -v $(pwd):/data -v $(pwd)/results:/results grinn \
  workflow /data/large_protein.pdb /results \
  --solvate --nt 8 --skip 10

# Very long trajectory analysis
docker run -v $(pwd):/data -v $(pwd)/results:/results grinn \
  workflow /data/protein.pdb /results \
  --traj /data/very_long_trajectory.xtc \
  --skip 100 \
  --nt 16

# Multiple simultaneous containers (test resource limits)
for i in {1..3}; do
  docker run -d --name grinn_test_$i \
    -v $(pwd)/test$i:/data -v $(pwd)/results$i:/results \
    grinn workflow /data/protein$i.pdb /results --nt 2
done

# Monitor all containers
docker stats $(docker ps --format "table {{.Names}}" | grep grinn_test)

# Force field compatibility testing
for ff in charmm36 amber99sb gromos54a7; do
  echo "Testing with $ff force field..."
  docker run -v $(pwd):/data -v $(pwd)/results_$ff:/results grinn \
    workflow /data/protein.pdb /results \
    --ff_folder /data/${ff}.ff --solvate
done
```

### Scientific Validation
```bash
# Compare gRINN results with experimental data
# 1. Known allosteric networks (e.g., PDZ domains, kinases)
docker run -v $(pwd):/data -v $(pwd)/results:/results grinn \
  workflow /data/pdz_domain.pdb /results \
  --create_pen --pen_cutoffs 0.5 1.0 1.5 2.0 \
  --source_sel "resid 10 to 20" \
  --target_sel "resid 80 to 90" \
  --solvate --npt

# 2. Drug binding site analysis
docker run -v $(pwd):/data -v $(pwd)/results:/results grinn \
  workflow /data/protein_ligand.pdb /results \
  --source_sel "resname LIG" \
  --target_sel "protein and within 5 of resname LIG" \
  --create_pen --solvate

# 3. Cross-validation with literature
# Run gRINN on published protein systems and compare energy networks
# with reported experimental interaction networks

# 4. Comparison with other MD analysis tools
# Run same system with VMD, PyMOL, or other analysis packages
# Compare energy decomposition and network metrics
```

## Getting Help

### Quick Questions
- [GitHub Discussions](https://github.com/osercinoglu/grinn-ismb-2025/discussions)
- Check existing [issues](https://github.com/osercinoglu/grinn/issues)

### Detailed Reports
- File [bug reports](https://github.com/osercinoglu/grinn/issues/new)
- Share results in [CONTRIBUTING.md](../CONTRIBUTING.md#reporting)

### During CollaborationFest
- Ask project mentors directly
- Join debugging sessions
- Coordinate with other testers

---

**Remember**: Every bug found and reported makes gRINN better for the entire community! 🚀
