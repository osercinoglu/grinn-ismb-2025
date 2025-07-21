# Testing Guide for gRINN

## Overview
This guide helps you systematically test gRINN with different protein systems and report your findings effectively.

## Prerequisites
- Docker installed and running
- gRINN Docker image built (`docker build -f Dockerfile.dev -t grinn-dev .`)
- Basic understanding of molecular dynamics simulations

## Quick Testing Checklist

### ✅ Basic Functionality Test
1. **Small protein test** (< 50 residues)
   - Use 1L2Y (Trp-cage) - fast and reliable
   - Generate 5 ns trajectory
   - Run gRINN analysis
   - Check dashboard loads

2. **Medium protein test** (50-150 residues)
   - Use 1UBQ (Ubiquitin) or 1LYS (Lysozyme)
   - Generate 10 ns trajectory
   - Run gRINN analysis
   - Compare performance

3. **Large protein test** (> 150 residues)
   - Use 1IGY (Antibody) or similar
   - Test memory usage and runtime
   - Check for crashes or errors

## Systematic Testing Protocol

### Step 1: Environment Test
```bash
# Verify Docker setup
docker run grinn-dev help

# Check GROMACS functionality
docker run grinn-dev gmx --version
```

### Step 2: Simulation Testing
Follow [quick-simulation.md](quick-simulation.md) for each test protein:

```bash
# Create test directory
mkdir grinn_test_PROTEIN && cd grinn_test_PROTEIN

# Download protein (replace PDBID)
wget https://files.rcsb.org/download/PDBID.pdb

# Run simulation pipeline
docker run -v $(pwd):/data -w /data grinn-dev gmx pdb2gmx -f PDBID.pdb -o protein.gro -p topol.top -ff amber99sb-ildn -water tip3p
# ... complete simulation steps
```

### Step 3: gRINN Analysis Testing
```bash
docker run -v $(pwd):/data grinn-dev workflow /data/em.gro /data/results 
  --top /data/topol.top --traj /data/protein_traj.xtc

# Strict cutoffs for dense networks
docker run -v $(pwd):/data grinn-dev workflow /data/em.gro /data/results_strict 
  --top /data/topol.top --traj /data/protein_traj.xtc --cutoff -9

# Fast processing with simplified algorithm
docker run -v $(pwd):/data grinn-dev workflow /data/em.gro /data/results_fast 
  --top /data/topol.top --traj /data/protein_traj.xtc --method bfs
```

### Step 4: Dashboard Testing
```bash
# Launch dashboard
docker run -p 8051:8051 -v $(pwd):/data grinn-dev dashboard /data/results

# Test in browser at http://localhost:8051
# Check: plots load, interactions are reasonable, navigation works
```

## What to Test and Report

### ✅ Success Criteria
- [ ] Analysis completes without errors
- [ ] Output files are generated (`interaction_energies.csv`, etc.)
- [ ] Dashboard loads and displays data
- [ ] Interaction energies are in reasonable range (-50 to +10 kJ/mol typically)
- [ ] Network visualization makes biological sense

### 🐛 Common Issues to Watch For
- Memory errors with large systems
- File format incompatibilities
- Dashboard loading failures
- Unreasonable energy values
- Missing or corrupted output files
- Docker permission issues

### 📊 Performance Metrics to Track
- **Runtime**: Record analysis time vs protein size
- **Memory usage**: Monitor peak memory consumption
- **File sizes**: Note output file sizes
- **Success rate**: Track which proteins work vs fail

## Testing Different Scenarios

### File Format Testing
```bash
# Test with different input formats
docker run -v $(pwd):/data grinn-dev workflow /data/protein.gro /data/results_gro
docker run -v $(pwd):/data grinn-dev workflow /data/protein.pdb /data/results_pdb
# Compare results for consistency
```

### Parameter Testing
```bash
# Different energy cutoffs
for cutoff in -5 -10 -15 -20; do
  docker run -v $(pwd):/data grinn-dev workflow /data/em.gro /data/results_cut$cutoff \
    --top /data/topol.top --traj /data/protein_traj.xtc --cutoff $cutoff
done
```

### Trajectory Length Testing
```bash
# Test with different trajectory lengths
# Short (1-2 ns), Medium (5-10 ns), Long (20+ ns)
# Monitor how analysis time scales
```

## Reporting Your Results

### For Successful Tests
Create a report including:
- **Protein tested**: PDB ID, size, type
- **Simulation details**: Length, force field, box size
- **gRINN performance**: Runtime, memory usage
- **Key findings**: Interesting interactions, network features
- **Files produced**: List of output files and sizes

### For Failed Tests
Report should include:
- **Full error message**: Copy-paste complete error
- **System details**: OS, Docker version, available memory
- **Input files**: Protein used, simulation parameters
- **Steps to reproduce**: Exact commands that led to failure
- **Partial results**: Any files that were created before failure

### For Performance Issues
Document:
- **Resource usage**: CPU, memory, disk space
- **Scaling behavior**: How performance varies with protein size
- **Bottlenecks identified**: Which steps are slowest
- **Suggested improvements**: Ideas for optimization

## Advanced Testing Scenarios

### Edge Cases
- Very small proteins (< 20 residues)
- Very large proteins (> 500 residues)
- Membrane proteins
- Multi-chain complexes
- Proteins with non-standard residues

### Stress Testing
- Very long trajectories (> 100 ns)
- High-resolution trajectories (many frames)
- Multiple parallel analyses
- Low-memory environments

### Integration Testing
- Test complete workflow from PDB to dashboard
- Test with different GROMACS versions
- Test on different operating systems
- Test with different Docker configurations

## Automated Testing (For Developers)

### Test Suite Creation
Create scripts that automatically test multiple scenarios:

```bash
#!/bin/bash
# auto_test.sh - Automated gRINN testing

proteins=("1L2Y" "1UBQ" "1LYS")
for protein in "${proteins[@]}"; do
  echo "Testing $protein..."
  mkdir test_$protein && cd test_$protein
  
  # Download and simulate
  wget https://files.rcsb.org/download/$protein.pdb
  # ... run simulation pipeline
  
  # Test gRINN
  docker run -v $(pwd):/data grinn-dev workflow /data/em.gro /data/results \
    --top /data/topol.top --traj /data/protein_traj.xtc || echo "FAILED: $protein"
    
  cd ..
done
```

## Quality Assurance

### Output Validation
- Check that energy values are physically reasonable
- Verify network topology makes biological sense
- Confirm dashboard visualizations are correct
- Test that saved results can be reloaded

### Cross-Validation
- Compare results with published data when available
- Test same protein with different simulation parameters
- Verify reproducibility of results

## Resources for Testers
- [Test Systems Guide](test-systems.md) - Recommended proteins
- [Quick Simulation Guide](quick-simulation.md) - How to generate trajectories
- [Development Guide](development-guide.md) - For deeper testing

Remember: Every test is valuable, even (especially!) the ones that fail! 🧪
