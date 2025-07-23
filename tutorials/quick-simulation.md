# Quick Simulation Guide for gRINN Testing

## Overview
If you don't have existing MD trajectories, this guide helps you generate a quick simulation for testing gRINN. **Note**: For CollaborationFest testing, the built-in Docker test data is usually sufficient.

## When You Need This Guide
- Testing gRINN with your own protein structure
- Comparing different simulation parameters
- Generating custom trajectories for specific research questions
- Learning the complete MD → gRINN workflow

## Prerequisites
- Docker with gRINN image built
- A protein structure (PDB file)
- Basic familiarity with molecular dynamics concepts

## Quick Simulation Workflow

### Option 1: Using gRINN's Built-in Simulation (Recommended)

```bash
# gRINN includes GROMACS and can run simulations automatically
# This is the easiest approach for testing

# Download a protein
wget https://files.rcsb.org/download/1L2Y.pdb

# Run complete workflow: simulation + analysis
mkdir results
docker run -v $(pwd):/data -v $(pwd)/results:/results grinn \
  workflow /data/1L2Y.pdb /results --gpu --solvate --npt

# This will:
# 1. Fix PDB file
# 2. Add solvent (water)
# 3. Run energy minimization
# 4. Run NPT equilibration
# 5. Calculate interaction energies
# 6. Generate network analysis
```

### Option 2: External Simulation → gRINN Analysis

If you prefer to run simulations separately:

```bash
# 1. Generate topology and coordinates with GROMACS
docker run -v $(pwd):/data gromacs/gromacs:2024 \
  gmx pdb2gmx -f /data/protein.pdb -o /data/protein.gro -p /data/topol.top

# 2. Add solvent and run simulation (simplified example)
# [Full GROMACS workflow - see GROMACS tutorials for details]

# 3. Use existing trajectory with gRINN
docker run -v $(pwd):/data -v $(pwd)/results:/results grinn \
  workflow /data/protein.gro /results --top /data/topol.top --traj /data/trajectory.xtc
```

## Simulation Parameters for Testing

### Quick Test (Development/Debugging)
```bash
# Minimal simulation for fast testing
--nsteps 1000      # 2 ps simulation
--skip 1           # Analyze every frame
```

### Standard Test (Validation)
```bash
# Short but meaningful simulation
--nsteps 50000     # 100 ps simulation  
--skip 5           # Analyze every 5th frame
```

### Production Test (Research)
```bash
# Longer simulation for real results
--nsteps 500000    # 1 ns simulation
--skip 10          # Analyze every 10th frame
```

## Input File Requirements

### PDB File Preparation
```bash
# Download from PDB
wget https://files.rcsb.org/download/PDBID.pdb

# Or clean existing PDB (remove water, ligands if needed)
grep "^ATOM" protein_original.pdb > protein_clean.pdb
```

### File Format Support
- **PDB files**: Standard protein structure format
- **GRO files**: GROMACS coordinate format (preferred for analysis)
- **Topology files**: Required for custom trajectories

## Common Issues & Solutions

### 🚨 Simulation Fails to Start
```bash
# Check PDB file format
head -10 protein.pdb
# Should start with HEADER, ATOM, or MODEL lines

# Try with PDB fixer enabled (default)
docker run -v $(pwd):/data -v $(pwd)/results:/results grinn \
  workflow /data/protein.pdb /results --fixpdb
```

### 🚨 Simulation Too Slow
```bash
# Use smaller system
grep "^ATOM" large_protein.pdb | head -500 > small_protein.pdb

# Skip solvation for testing
docker run -v $(pwd):/data -v $(pwd)/results:/results grinn \
  workflow /data/protein.pdb /results --no-solvate

# Use CPU only if GPU causes issues
docker run -v $(pwd):/data -v $(pwd)/results:/results grinn \
  workflow /data/protein.pdb /results --no-gpu
```

### 🚨 Memory Issues
```bash
# Reduce trajectory analysis frequency
docker run -v $(pwd):/data -v $(pwd)/results:/results grinn \
  workflow /data/protein.pdb /results --skip 20

# Use smaller cutoff for initial filtering
docker run -v $(pwd):/data -v $(pwd)/results:/results grinn \
  workflow /data/protein.pdb /results --initpairfiltercutoff 5.0
```

## Validation Checklist

### ✅ Simulation Quality
- [ ] System equilibrated (stable temperature/pressure)
- [ ] No clashes or extreme forces
- [ ] Protein structure remains stable
- [ ] Trajectory looks reasonable in viewer

### ✅ gRINN Analysis
- [ ] Energy matrices populated
- [ ] Interaction energies in reasonable range (-10 to +5 kcal/mol)
- [ ] Strong interactions make chemical sense
- [ ] Networks show expected connectivity

## Example Test Systems

### Tiny (< 1 minute)
- **Trp-cage** (1L2Y): 20 residues, very fast
- **Mini-protein** (1WRP): 23 residues

### Small (< 5 minutes)
- **Villin headpiece** (1YRF): 35 residues
- **Protein G** (1PGB): 56 residues

### Medium (< 30 minutes)
- **Ubiquitin** (1UBQ): 76 residues
- **Lysozyme** (1LYS): 129 residues

## Tips for Efficient Testing

### 🚀 Speed Up Testing
1. **Start small**: Use proteins < 50 residues for initial tests
2. **Short simulations**: 100 ps is enough for functionality testing
3. **Skip frames**: Analyze every 10th frame to reduce computation
4. **No solvation**: Skip water for pure protein analysis

### 🎯 Focus on Validation
1. **Known systems**: Use well-studied proteins with known interactions
2. **Compare tools**: Cross-check results with other analysis methods
3. **Literature validation**: Compare networks with experimental data
4. **Chemical intuition**: Strong interactions should make sense

### 📊 Systematic Testing
1. **Document everything**: Note parameters and results
2. **Test systematically**: Vary one parameter at a time
3. **Share results**: Report both successes and failures
4. **Help others**: Add working examples to the results folder

## Next Steps

**Simulation successful?**
- ✅ Proceed to [quickstart guide](quickstart.md) for analysis
- ✅ Try the interactive dashboard
- ✅ Test with larger systems

**Encountered issues?**
- 🔧 Check [testing guide](testing-guide.md) for troubleshooting
- 🐛 Report bugs via [CONTRIBUTING.md](../CONTRIBUTING.md#reporting)
- 💬 Ask for help in GitHub discussions

**Want to go deeper?**
- 📖 Read [development guide](development-guide.md) for customization
- 🧪 Try advanced analysis parameters
- 🎯 Validate with your research systems

---

**Remember**: The goal is testing gRINN, not perfect simulations! Short, quick runs are perfectly adequate for finding bugs and validating functionality. 🚀
