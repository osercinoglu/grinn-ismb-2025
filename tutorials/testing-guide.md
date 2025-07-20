# Testing Guide - gRINN CollaborationFest

This guide will help you test gRINN effectively and contribute valuable feedback to improve the tool.

## Overview

As a tester, your role is crucial for:
- Validating gRINN works with different protein systems
- Identifying bugs and edge cases
- Providing feedback on usability
- Ensuring results are scientifically meaningful

## Prerequisites

- Basic understanding of molecular dynamics simulations
- Familiarity with protein structures (PDB format)
- GROMACS installation (optional, for your own data)

## Getting Started

### 1. Set Up gRINN

Follow the installation instructions in the [main repository](https://github.com/osercinoglu/grinn#installation--usage).

**Quick Docker Setup** (Recommended):
```bash
# Pull the Docker image
docker pull osercinoglu/grinn:latest

# Test with sample data
docker run -p 8051:8051 osercinoglu/grinn dashboard test
```

### 2. Test with Sample Data

We provide three sample datasets for testing:

#### Small Protein System (< 100 residues)
**Purpose**: Quick testing and validation
**What to check**:
- Tool runs without errors
- Results are generated in reasonable time
- Dashboard loads correctly

```bash
# Example command (adjust paths as needed)
python grinn_workflow.py sample_small.pdb results_small/ \
  --top sample_small.top \
  --traj sample_small.xtc \
  --nt 4
```

#### Protein-Ligand Complex
**Purpose**: Test protein-ligand interaction analysis
**What to check**:
- Ligand interactions are captured correctly
- Residue selection works for protein and ligand
- Energy calculations include all relevant interactions

```bash
# Example with residue selection
python grinn_workflow.py protein_ligand.pdb results_ligand/ \
  --top protein_ligand.top \
  --traj protein_ligand.xtc \
  --source_sel "protein" \
  --target_sel "resname LIG" \
  --nt 4
```

#### Large Protein Complex (> 500 residues)
**Purpose**: Test scalability and performance
**What to check**:
- Memory usage remains reasonable
- Processing time is acceptable
- Results are complete

```bash
# Example with frame skipping for efficiency
python grinn_workflow.py large_complex.pdb results_large/ \
  --top large_complex.top \
  --traj large_complex.xtc \
  --skip 10 \
  --nt 8
```

## Testing Checklist

### Basic Functionality
- [ ] Installation completes without errors
- [ ] Tool runs with sample data
- [ ] Output files are generated correctly
- [ ] Dashboard launches and displays data
- [ ] Docker container works on your system

### Output Validation
- [ ] Energy CSV files contain expected data
- [ ] Network files (GML) can be opened in network analysis tools
- [ ] Dashboard visualizations appear correct
- [ ] Log files show no critical errors

### Performance Testing
- [ ] Small systems complete in < 5 minutes
- [ ] Large systems complete in reasonable time
- [ ] Memory usage doesn't exceed system limits
- [ ] Frame skipping works as expected

### Edge Cases
- [ ] Empty or invalid input files are handled gracefully
- [ ] Missing trajectory files produce clear error messages
- [ ] Unusual residue names don't cause crashes
- [ ] Very small or very large systems work correctly

## Using Your Own Data

### Data Requirements
- **PDB file**: Protein structure (cleaned, no water molecules for dry analysis)
- **Topology file**: GROMACS topology (.top)
- **Trajectory file**: GROMACS trajectory (.xtc, .trr)

### Best Practices
1. **Start small**: Test with a short trajectory first
2. **Check topology**: Ensure topology matches your structure
3. **Use frame skipping**: For long trajectories, use `--skip` parameter
4. **Monitor resources**: Watch memory and CPU usage

### Common Issues and Solutions

| Issue | Possible Cause | Solution |
|-------|---------------|----------|
| Out of memory | Large system/trajectory | Use `--skip` parameter, reduce system size |
| Topology errors | Mismatched PDB/topology | Regenerate topology with GROMACS |
| No interactions found | High distance cutoff | Check `--initpairfiltercutoff` parameter |
| Dashboard won't load | Port conflict | Check if port 8051 is available |

## Reporting Issues

### What to Report
- **Bugs**: Unexpected crashes, incorrect results, error messages
- **Performance issues**: Slow execution, high memory usage
- **Usability problems**: Confusing interfaces, unclear documentation

### How to Report
1. **Check existing issues** on [GitHub](https://github.com/osercinoglu/grinn/issues)
2. **Create a new issue** with:
   - Clear title describing the problem
   - Steps to reproduce the issue
   - Your system information (OS, Python version, etc.)
   - Error messages or unexpected output
   - Sample data if relevant (small files only)

### Issue Template
```
**Problem Description**
Brief description of what went wrong

**Steps to Reproduce**
1. Command used: `python grinn_workflow.py ...`
2. Input files: [describe your data]
3. Expected behavior: [what should happen]
4. Actual behavior: [what actually happened]

**System Information**
- OS: [Windows/macOS/Linux distribution]
- Python version: [e.g., 3.10.2]
- gRINN version/commit: [if known]
- GROMACS version: [if applicable]

**Error Messages**
```
[paste any error messages here]
```

**Additional Context**
Any other relevant information
```

## Testing Results Documentation

### What to Document
- Systems tested (size, type, source)
- Commands used
- Processing times
- Memory usage
- Any issues encountered
- Suggestions for improvement

### Sharing Results
- Post successful analyses in GitHub discussions
- Share interesting findings or visualizations
- Contribute example workflows for others

## Advanced Testing

### Custom Residue Selections
Test different selection syntaxes:
```bash
# Chain-specific analysis
--source_sel "chain A" --target_sel "chain B"

# Residue type selection
--source_sel "resname ARG LYS" --target_sel "protein"

# Complex selections
--source_sel "chain A and resid 1:100" --target_sel "chain B and resname TRP PHE"
```

### Performance Benchmarking
For systematic performance testing:
```bash
# Time the execution
time python grinn_workflow.py [your_command]

# Monitor memory usage (Linux/macOS)
/usr/bin/time -v python grinn_workflow.py [your_command]
```

### Dashboard Testing
- Test with different browsers (Chrome, Firefox, Safari)
- Try different screen sizes and resolutions
- Upload different result sets
- Test all dashboard features and visualizations

## Next Steps

After testing:
1. **Report findings** via GitHub issues or discussions
2. **Share successful workflows** with the community
3. **Consider contributing** to documentation or development
4. **Help other testers** by answering questions

## Resources

- [Main gRINN Repository](https://github.com/osercinoglu/grinn)
- [GROMACS Documentation](https://manual.gromacs.org/)
- [ProDy Selection Syntax](http://prody.csb.pitt.edu/manual/reference/atomic/select.html)
- [GitHub Issues](https://github.com/osercinoglu/grinn/issues)

Happy testing! Your feedback helps make gRINN better for everyone. 🧪
