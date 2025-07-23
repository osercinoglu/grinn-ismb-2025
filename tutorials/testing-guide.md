# Testing Guide for gRINN

## Overview
This guide helps you systematically test gRINN and report useful feedback to improve the tool for everyone.

## Quick Testing Checklist

### ✅ Basic Functionality Test
- [ ] Docker image builds successfully
- [ ] Test data runs without errors
- [ ] Dashboard launches and displays results
- [ ] Energy matrices look reasonable
- [ ] Network plots are generated

### ✅ Performance Test
- [ ] Note runtime for small protein (< 50 residues)
- [ ] Note runtime for medium protein (50-150 residues)
- [ ] Monitor memory usage during calculation
- [ ] Test with/without GPU acceleration

### ✅ Edge Cases
- [ ] Multi-chain proteins
- [ ] Proteins with missing residues
- [ ] Very small proteins (< 20 residues)
- [ ] Proteins with unusual amino acids
- [ ] Empty or corrupted input files

## Systematic Testing Protocol

### Step 1: Environment Setup
```bash
# Test Docker installation
docker --version
docker info

# Build gRINN image
git clone https://github.com/osercinoglu/grinn.git
cd grinn
docker build -t grinn .

# Verify build
docker images | grep grinn
```

### Step 2: Functional Testing
```bash
# Test with built-in data
mkdir test_results
docker run -v $(pwd)/test_results:/results grinn test

# Verify outputs exist
ls test_results/
# Should contain: energies_*.csv, system_dry.pdb, traj_dry.xtc

# Test dashboard
docker run -p 8051:8051 -v $(pwd)/test_results:/results grinn dashboard /results
# Check: http://localhost:8051
```

### Step 3: Performance Benchmarking
```bash
# Small protein test (Trp-cage)
wget https://files.rcsb.org/download/1L2Y.pdb
time docker run -v $(pwd):/data -v $(pwd)/results_small:/results grinn \
  workflow /data/1L2Y.pdb /results --gpu --solvate

# Medium protein test (Lysozyme)
wget https://files.rcsb.org/download/1LYS.pdb
time docker run -v $(pwd):/data -v $(pwd)/results_medium:/results grinn \
  workflow /data/1LYS.pdb /results --gpu --solvate

# Note runtimes and memory usage
```

### Step 4: Validation Testing
```bash
# Check energy value ranges
head -5 results_small/energies_intEnTotal.csv

# Verify network makes sense
# Strong interactions should be chemically reasonable
# (salt bridges, hydrogen bonds, hydrophobic contacts)
```

## What to Look For

### ✅ Success Indicators
- **Energy values**: Typically -10 to +5 kcal/mol range
- **Strong interactions**: < -2 kcal/mol, chemically reasonable
- **Network hubs**: Often at active sites or structural centers
- **Dashboard responsiveness**: Smooth interaction with plots

### 🚨 Warning Signs
- **Extreme energies**: > +10 or < -20 kcal/mol (possible errors)
- **Empty outputs**: Missing CSV files or zero-size files
- **Dashboard errors**: JavaScript errors in browser console
- **Memory issues**: Docker container killed due to OOM

### 🐛 Common Issues
- **Docker build failures**: Check Docker version and disk space
- **Permission errors**: Ensure output directories are writable
- **GROMACS errors**: Often due to malformed PDB files
- **Dashboard won't load**: Check port conflicts (8051)

## Reporting Results

### Success Report Template
```markdown
## Test Report: [Protein Name]

**System**: [OS, Docker version]
**Protein**: [PDB ID, size, description]
**Runtime**: [Total time from start to finish]
**Memory**: [Peak memory usage if known]

### Results
- Energy range: [min to max kcal/mol]
- Network hubs: [list key residues]
- Notable interactions: [interesting findings]

### Performance
- Build time: [X minutes]
- Simulation time: [X minutes] 
- Analysis time: [X minutes]

### Issues
- [None / list any problems]
```

### Bug Report Template
```markdown
## Bug Report: [Brief Description]

**Environment**:
- OS: [Windows/Mac/Linux + version]
- Docker: [version]
- System: [CPU, RAM, GPU if relevant]

**Steps to Reproduce**:
1. [Exact commands used]
2. [Input files used]
3. [When error occurred]

**Expected Behavior**:
[What should have happened]

**Actual Behavior**:
[What actually happened]

**Error Messages**:
```
[Full error output]
```

**Files**:
- Input PDB: [attach or link]
- Docker logs: [attach if helpful]
- Screenshots: [if dashboard related]
```

## Advanced Testing

### Custom Data Testing
- Test with your own MD trajectories
- Compare gRINN results with other tools
- Validate against known protein interactions

### Stress Testing
- Very large proteins (>500 residues)
- Long trajectories (>1000 frames)
- Multiple simultaneous runs
- Limited memory conditions

### Scientific Validation
- Compare with experimental data
- Test known allosteric networks
- Validate drug binding sites
- Cross-check with literature

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
