# Sample Data and Examples - gRINN CollaborationFest

This guide describes the sample datasets provided for testing gRINN and includes example commands for different analysis scenarios.

## Available Sample Datasets

### 1. Small Protein System (villin-headpiece)
- **File**: `small_protein.pdb`, `small_protein.top`, `small_protein.xtc`
- **Description**: 35-residue villin headpiece subdomain (HP-35)
- **Size**: ~35 residues
- **Runtime**: 2-5 minutes
- **Purpose**: Quick testing, validation, learning

**Example Analysis**:
```bash
python grinn_workflow.py small_protein.pdb results_small/ \
  --top small_protein.top \
  --traj small_protein.xtc \
  --nt 4
```

### 2. Protein-Ligand Complex (lysozyme-inhibitor)
- **File**: `protein_ligand.pdb`, `protein_ligand.top`, `protein_ligand.xtc`
- **Description**: Lysozyme with bound inhibitor
- **Size**: ~129 residues + 1 ligand
- **Runtime**: 5-10 minutes
- **Purpose**: Testing protein-ligand interactions

**Example Analysis**:
```bash
# Analyze protein-ligand interactions
python grinn_workflow.py protein_ligand.pdb results_ligand/ \
  --top protein_ligand.top \
  --traj protein_ligand.xtc \
  --source_sel "protein" \
  --target_sel "resname LIG" \
  --nt 4

# Analyze protein internal interactions
python grinn_workflow.py protein_ligand.pdb results_protein/ \
  --top protein_ligand.top \
  --traj protein_ligand.xtc \
  --source_sel "protein" \
  --target_sel "protein" \
  --nt 4
```

### 3. Large Protein Complex (antibody Fab fragment)
- **File**: `large_complex.pdb`, `large_complex.top`, `large_complex.xtc`
- **Description**: Antibody Fab fragment (heavy + light chains)
- **Size**: ~450 residues
- **Runtime**: 20-60 minutes (depending on trajectory length)
- **Purpose**: Testing scalability and performance

**Example Analysis**:
```bash
# Full analysis with frame skipping
python grinn_workflow.py large_complex.pdb results_large/ \
  --top large_complex.top \
  --traj large_complex.xtc \
  --skip 10 \
  --nt 8

# Chain-specific analysis
python grinn_workflow.py large_complex.pdb results_chains/ \
  --top large_complex.top \
  --traj large_complex.xtc \
  --source_sel "chain A" \
  --target_sel "chain B" \
  --skip 5 \
  --nt 8
```

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
