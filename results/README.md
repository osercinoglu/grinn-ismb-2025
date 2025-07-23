# gRINN Test Results

This directory contains test results, benchmarks, and validation data from community testing of gRINN.

## Directory Structure

```
results/
├── README.md           # This file
├── benchmarks/         # Performance testing results
├── validation/         # Scientific validation with known systems
├── bug-reports/        # Detailed bug reports with reproduction steps
└── success-stories/    # Working examples and case studies
```

## How to Contribute Results

### 1. Success Stories
Create a folder with your protein name and include:
- Input files used
- gRINN output files
- Screenshots of dashboard
- Brief description of findings

### 2. Benchmarks
Include:
- System specifications
- Protein size and simulation details
- Runtime measurements
- Memory usage
- Performance comparison data

### 3. Bug Reports
Provide:
- Full error messages
- Input files that caused the issue
- System information
- Steps to reproduce

### 4. Validation Studies
For known systems:
- Comparison with literature
- Expected vs. observed interactions
- Network analysis validation
- Scientific interpretation

## Example Structure

```
results/
└── lysozyme_1LYS/
    ├── input/
    │   ├── 1LYS.pdb
    │   └── simulation_params.txt
    ├── output/
    │   ├── energies_intEnTotal.csv
    │   ├── pen_networks.gml
    │   └── dashboard_screenshots/
    ├── analysis/
    │   └── findings.md
    └── metadata.json
```

## Sharing Guidelines

1. **Fork this repository**
2. **Add your results** to appropriate subfolder
3. **Include metadata** (protein, system size, runtime, etc.)
4. **Write a brief summary** of your findings
5. **Create pull request** with clear description

## Privacy

- Only share data you have permission to distribute
- Remove sensitive information from input files
- Consider using public PDB structures for reproducibility

---

**Questions?** Ask in [GitHub Discussions](https://github.com/osercinoglu/grinn-ismb-2025/discussions)
