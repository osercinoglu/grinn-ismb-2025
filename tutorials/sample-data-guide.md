# Suggested Test Systems - gRINN CollaborationFest

This guide provides a list of recommended protein systems for testing gRINN. Contributors can use these PDB structures to generate their own MD trajectories and test the tool with fresh input data.

## Why Use Different Test Systems?

The goal is to test gRINN with newly-generated input from different researchers, which helps:
- Identify issues with diverse protein types and sizes
- Test compatibility with different GROMACS setups
- Validate results across different simulation protocols
- Discover edge cases and unexpected behaviors

## Recommended Test Systems

### Small Proteins (< 100 residues)
**Purpose**: Quick testing, validation, learning the tool

| Protein | PDB ID | Residues | Description | Why Test? |
|---------|--------|----------|-------------|-----------|
| Villin Headpiece | `2F4K` | 35 | Fast-folding mini-protein | Very fast to simulate, well-studied |
| Trp-cage | `1L2Y` | 20 | Ultra-fast folding protein | Minimal system, good for debugging |
| WW Domain | `1PIN` | 37 | Small β-sheet protein | Tests β-sheet interactions |
| Chignolin | `1UAO` | 10 | Minimal β-hairpin | Extremely small test case |

### Medium Proteins (100-300 residues)
**Purpose**: Standard testing scenarios

| Protein | PDB ID | Residues | Description | Why Test? |
|---------|--------|----------|-------------|-----------|
| Lysozyme | `1LYS` | 129 | Well-studied enzyme | Classic test system, stable |
| Ubiquitin | `1UBQ` | 76 | Small regulatory protein | Compact, all-β structure |
| Protein G | `1PGB` | 56 | Immunoglobulin-binding | Mixed α/β structure |
| Crambin | `1CRN` | 46 | Plant seed protein | Contains disulfide bonds |

### Protein-Ligand Complexes
**Purpose**: Test protein-ligand interaction analysis

| Protein | PDB ID | Residues | Ligand | Why Test? |
|---------|--------|----------|--------|-----------|
| Lysozyme + inhibitor | `1LYZ` | 129 | Tri-NAG | Clear binding site |
| Trypsin + inhibitor | `1PPH` | 223 | Benzamidine | Strong binding |
| HIV Protease | `1HTM` | 198 | Inhibitor | Drug target |
| Carbonic Anhydrase | `1CA2` | 259 | Inhibitor | Well-defined active site |

### Large Systems (> 300 residues)
**Purpose**: Test scalability and performance

| Protein | PDB ID | Residues | Description | Why Test? |
|---------|--------|----------|-------------|-----------|
| Antibody Fab | `1IGY` | ~450 | Heavy + light chains | Multi-chain system |
| Barnase-Barstar | `1BRS` | ~250 | Protein-protein complex | Interface interactions |
| T4 Lysozyme | `1L63` | 164 | Larger enzyme | More complex topology |
| Immunoglobulin | `1IGT` | ~220 | Antibody fragment | β-sandwich structures |

### Challenging Systems
**Purpose**: Test edge cases and robustness

| Protein | PDB ID | Residues | Description | Challenge |
|---------|--------|----------|-------------|-----------|
| Membrane protein | `1OKC` | 348 | Rhodopsin | Membrane environment |
| Intrinsically disordered | `1F0R` | 140 | p53 TAD | Flexible regions |
| Multi-domain | `1TEN` | 89+76 | Fibronectin | Domain interfaces |
| Metal-binding | `1ZNC` | 65 | Zinc finger | Metal coordination |

## How to Use These Systems

### 1. Choose Your Test System
- **New users**: Start with small proteins (Villin, Trp-cage)
- **Experienced users**: Try medium or large systems
- **Method developers**: Use challenging systems

### 2. Download Structure
```bash
# Download PDB file
wget https://files.rcsb.org/download/1L2Y.pdb

# Or use your preferred method to get the structure
```

### 3. Prepare MD Simulation
- Clean the structure (remove waters, heterogens as needed)
- Generate topology with GROMACS
- Run your simulation protocol
- Generate trajectory for gRINN analysis

### 4. Test with gRINN
```bash
python grinn_workflow.py your_protein.pdb results/ \
  --top your_protein.top \
  --traj your_trajectory.xtc \
  --nt 4
```

## Testing Scenarios

### Basic Functionality Tests
- Does gRINN run without errors?
- Are output files generated correctly?
- Does the dashboard load and display results?

### Scientific Validation Tests
- Do energy values make physical sense?
- Are strong interactions captured correctly?
- Do network representations look reasonable?

### Performance Tests
- How long does analysis take?
- Memory usage with different system sizes?
- Effect of frame skipping on results?

### Robustness Tests
- Different force fields (AMBER, CHARMM, OPLS)
- Various simulation lengths
- Different trajectory formats (.xtc, .trr)
- Unusual residue names or modifications

## What to Report

### Successful Tests
- System details (protein, size, simulation parameters)
- Runtime and memory usage
- Interesting scientific findings
- Workflow that worked well

### Issues Found
- Error messages and when they occur
- Unexpected results or behaviors
- Performance problems
- Suggestions for improvement

### Comparison Studies
- Results from different force fields
- Effect of simulation length
- Comparison with experimental data
- Network analysis insights

## Creating Your Test Dataset

### Recommended Simulation Protocol
1. **Structure preparation**: Clean PDB, add hydrogens
2. **Solvation**: Add water box (if desired)
3. **Minimization**: Energy minimize the system
4. **Equilibration**: Short NVT/NPT equilibration
5. **Production**: 10-100 ns simulation (depending on system)
6. **Post-processing**: Extract protein-only trajectory

### GROMACS Example Workflow
```bash
# 1. Process structure
gmx pdb2gmx -f protein.pdb -o processed.pdb -p topol.top -ff amber99sb-ildn -water tip3p

# 2. Create box and solvate
gmx editconf -f processed.pdb -o boxed.pdb -c -d 1.0 -bt cubic
gmx solvate -cp boxed.pdb -cs spc216.gro -p topol.top -o solvated.pdb

# 3. Add ions if needed
gmx grompp -f ions.mdp -c solvated.pdb -p topol.top -o ions.tpr
echo "SOL" | gmx genion -s ions.tpr -o ions.pdb -p topol.top -pname NA -nname CL -neutral

# 4. Energy minimization
gmx grompp -f minim.mdp -c ions.pdb -p topol.top -o em.tpr
gmx mdrun -deffnm em

# 5. NVT equilibration
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt

# 6. NPT equilibration
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr
gmx mdrun -deffnm npt

# 7. Production MD
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -deffnm md

# 8. Extract protein trajectory for gRINN
echo "Protein" | gmx trjconv -s md.tpr -f md.xtc -o protein_traj.xtc -pbc mol -center
```

## Tips for Effective Testing

### Start Small
- Begin with small, well-behaved systems
- Use short trajectories initially
- Verify basic functionality before scaling up

### Document Everything
- Keep track of what you test
- Note any modifications to standard protocols
- Record both successes and failures

### Share Results
- Post successful workflows in discussions
- Report bugs with clear reproduction steps
- Share interesting scientific findings

### Collaborate
- Coordinate with others to avoid duplicate testing
- Share computational resources if possible
- Help troubleshoot others' issues

## Resources

### Structure Databases
- [RCSB PDB](https://www.rcsb.org/)
- [PDB Europe](https://www.ebi.ac.uk/pdbe/)
- [Membrane Proteins](https://blanco.biomol.uci.edu/mpstruc/)

### Simulation Tutorials
- [GROMACS Tutorials](http://www.mdtutorials.com/gmx/)
- [Best Practices](https://www.mdtutorials.com/gmx/lysozyme/index.html)
- [Force Field Comparison](https://www.mdtutorials.com/gmx/complex/index.html)

### Analysis Tools
- [VMD](https://www.ks.uiuc.edu/Research/vmd/)
- [PyMOL](https://pymol.org/)
- [ChimeraX](https://www.cgl.ucsf.edu/chimerax/)

## Getting Help

- **Technical issues**: GitHub issues on main repository
- **Simulation questions**: GROMACS user forum
- **gRINN usage**: GitHub discussions
- **CollaborationFest**: Ask coordinators during event

Happy testing with real systems! 🧬

## Example Workflows

### Basic Energy Analysis
```bash
# Simple energy calculation
python grinn_workflow.py small_protein.pdb basic_output/ \
  --top small_protein.top \
  --traj small_protein.xtc \
  --nt 4

# Check outputs
ls basic_output/
head basic_output/energies_intEnTotal.csv
```

### Protein Energy Network (PEN) Analysis
```bash
# Generate PENs with multiple cutoffs
python grinn_workflow.py small_protein.pdb pen_output/ \
  --top small_protein.top \
  --traj small_protein.xtc \
  --create_pen \
  --pen_cutoffs 1.0 2.0 3.0 \
  --nt 4

# Outputs include:
# - Network files (.gml)
# - Betweenness centrality analysis
# - Dashboard-ready files
```

### Advanced Selection Examples
```bash
# Specific residue types
python grinn_workflow.py protein_ligand.pdb aromatic_output/ \
  --top protein_ligand.top \
  --traj protein_ligand.xtc \
  --source_sel "resname PHE TRP TYR" \
  --target_sel "protein" \
  --nt 4

# Binding site analysis
python grinn_workflow.py protein_ligand.pdb binding_site/ \
  --top protein_ligand.top \
  --traj protein_ligand.xtc \
  --source_sel "protein and within 8 of resname LIG" \
  --target_sel "resname LIG" \
  --nt 4

# Loop-domain interactions
python grinn_workflow.py large_complex.pdb loop_analysis/ \
  --top large_complex.top \
  --traj large_complex.xtc \
  --source_sel "resid 50:75" \
  --target_sel "resid 100:150" \
  --skip 5 \
  --nt 6
```

### Dashboard Visualization
```bash
# After running analysis, launch dashboard
python gRINN_Dashboard/grinn_dashboard.py results_folder/

# Open browser to http://localhost:8051
# Upload and explore your results interactively
```

## Performance Optimization Examples

### Memory-Conscious Analysis
```bash
# For systems that might use too much memory
python grinn_workflow.py large_complex.pdb memory_safe/ \
  --top large_complex.top \
  --traj large_complex.xtc \
  --skip 20 \
  --initpairfiltercutoff 8.0 \
  --nt 4
```

### High-Performance Computing
```bash
# For clusters or high-end workstations
python grinn_workflow.py large_complex.pdb hpc_output/ \
  --top large_complex.top \
  --traj large_complex.xtc \
  --skip 5 \
  --nt 16 \
  --create_pen \
  --pen_cutoffs 0.5 1.0 1.5 2.0 2.5
```

## Testing Different Scenarios

### Error Handling Tests
```bash
# Test with missing files (should fail gracefully)
python grinn_workflow.py nonexistent.pdb output/

# Test with mismatched topology
python grinn_workflow.py small_protein.pdb output/ \
  --top protein_ligand.top \
  --traj small_protein.xtc

# Test input validation
python grinn_workflow.py small_protein.pdb output/ \
  --top small_protein.top \
  --traj small_protein.xtc \
  --test_only
```

### Edge Cases
```bash
# Very restrictive cutoff (should find few interactions)
python grinn_workflow.py small_protein.pdb restrictive/ \
  --top small_protein.top \
  --traj small_protein.xtc \
  --initpairfiltercutoff 3.0

# Single residue selection
python grinn_workflow.py small_protein.pdb single_res/ \
  --top small_protein.top \
  --traj small_protein.xtc \
  --source_sel "resid 10" \
  --target_sel "protein"
```

## Expected Output Files

### Basic Analysis Output
```
results/
├── calc.log                           # Detailed log file
├── energies_intEnTotal.csv             # Total interaction energies
├── energies_intEnElec.csv              # Electrostatic energies
├── energies_intEnVdW.csv               # van der Waals energies
├── energies_intEnTotal.pickle          # Binary format
├── system_dry.pdb                      # Processed structure
├── traj_superposed.xtc                 # Aligned trajectory
├── topol_dry.top                       # Processed topology
└── grinn_workflow_summary.json         # Analysis summary
```

### PEN Analysis Output
```
results/
├── [basic files above]
├── pen_cov_True_cutoff_1.0_frame_0.gml    # Network files
├── pen_cov_False_cutoff_2.0_frame_0.gml   # (various combinations)
└── pen_betweenness_centralities.csv       # Centrality analysis
```

## Data Interpretation Guide

### Energy CSV Format
```csv
pair,frame_0,frame_1,frame_2,...,res1_index,res2_index,res1_chain,res2_chain
0-1,-5.23,-4.87,-6.12,...,0,1,A,A
0-2,-2.11,-1.95,-2.34,...,0,2,A,A
```

### Understanding Energy Values
- **Negative values**: Attractive interactions (favorable)
- **Positive values**: Repulsive interactions (unfavorable)
- **Magnitude**: Strength of interaction (|E| > 1 kcal/mol often significant)

### Network Analysis
- **Nodes**: Residues in the protein
- **Edges**: Significant interactions (above cutoff)
- **Betweenness centrality**: Importance in communication pathways

## Troubleshooting Sample Data

### Common Issues

#### File Not Found
```bash
# Ensure you're in the correct directory
ls sample_data/
# Should show .pdb, .top, and .xtc files
```

#### GROMACS Errors
```bash
# Check GROMACS installation
gmx --version

# Verify topology format
head -20 sample_data/small_protein.top
```

#### Memory Issues
```bash
# Use smaller test case first
python grinn_workflow.py small_protein.pdb test/ \
  --top small_protein.top \
  --traj small_protein.xtc \
  --nt 2
```

## Creating Your Own Test Cases

### Preparing MD Data
1. **Clean PDB structure**: Remove water, ions (for dry analysis)
2. **Generate topology**: Use GROMACS tools
3. **Run short simulation**: 1-10 ns sufficient for testing
4. **Save trajectory**: Use .xtc format for efficiency

### Example GROMACS Workflow
```bash
# Generate topology
gmx pdb2gmx -f protein.pdb -o protein_processed.pdb -p topol.top

# Create simulation box
gmx editconf -f protein_processed.pdb -o boxed.pdb -c -d 1.0

# Run energy minimization
gmx grompp -f minim.mdp -c boxed.pdb -p topol.top -o em.tpr
gmx mdrun -deffnm em

# Short MD simulation
gmx grompp -f md.mdp -c em.gro -p topol.top -o md.tpr
gmx mdrun -deffnm md

# Extract trajectory for gRINN
gmx trjconv -f md.xtc -o trajectory_for_grinn.xtc -s md.tpr
```

## Performance Benchmarks

### Expected Runtimes (4 CPU cores)
- **Small protein (35 res, 100 frames)**: 2-5 minutes
- **Protein-ligand (129 res, 200 frames)**: 5-10 minutes
- **Large complex (450 res, 500 frames, skip=10)**: 20-40 minutes

### Memory Usage
- **Small systems**: < 1 GB RAM
- **Medium systems**: 1-4 GB RAM
- **Large systems**: 4-16 GB RAM (with frame skipping)

## Next Steps

After testing with sample data:
1. **Try your own systems**: Apply gRINN to your research data
2. **Explore parameters**: Test different cutoffs and selections
3. **Contribute feedback**: Report any issues or suggestions
4. **Share results**: Discuss interesting findings with the community

## Resources

- [Main gRINN Repository](https://github.com/osercinoglu/grinn)
- [GROMACS Tutorials](http://www.mdtutorials.com/gmx/)
- [ProDy Selection Syntax](http://prody.csb.pitt.edu/manual/reference/atomic/select.html)
- [Network Analysis with NetworkX](https://networkx.org/documentation/stable/tutorial.html)

Happy analyzing! 🧬
