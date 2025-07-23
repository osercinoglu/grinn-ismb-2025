# gRINN Quickstart Guide

**Complete workflow from protein to results in ~5 minutes**

## Prerequisites
- Docker installed on your system
- Basic command line familiarity

## Step 1: Get gRINN Ready (1 min)

```bash
# Clone and build Docker image
git clone https://github.com/osercinoglu/grinn.git
cd grinn
docker build -t grinn .
```

## Step 2: Get Test Data (30 seconds)

```bash
# Download a small test protein (Trp-cage - 20 residues)
wget https://files.rcsb.org/download/1L2Y.pdb

# OR use built-in test data
docker run grinn test --help
```

## Step 3: Run gRINN (2-3 min)

```bash
# Create output directory
mkdir results

# Run full workflow with built-in test data
docker run -v $(pwd)/results:/results grinn test

# OR run with your downloaded protein
docker run -v $(pwd):/data -v $(pwd)/results:/results grinn \
  workflow /data/1L2Y.pdb /results --gpu --solvate --npt
```

## Step 4: View Results (30 seconds)

```bash
# Launch interactive dashboard
docker run -p 8051:8051 -v $(pwd)/results:/results grinn dashboard /results

# Open browser to: http://localhost:8051
```

## What You Should See

✅ **Energy matrices**: Pairwise interaction energies between residues  
✅ **Network plots**: Centrality measures for protein communication  
✅ **3D visualization**: Interactive molecular viewer with trajectories  

## Quick Validation

**Sanity checks for results:**
- Energy values typically range from -10 to +5 kcal/mol
- Strong interactions (< -2 kcal/mol) should make chemical sense
- Network hubs often correspond to structurally important residues

## Troubleshooting

**Container issues?**
```bash
# Check Docker installation
docker --version

# Check image was built
docker images | grep grinn
```

**Slow performance?**
- Use `--gpu` flag if NVIDIA GPU available
- Start with small proteins (< 100 residues)
- Reduce trajectory length with `--skip 10` (analyze every 10th frame)

**No results?**
- Check output directory permissions
- Verify input PDB file is valid
- Look for error messages in container logs

## Next Steps

**Working well?** 
- Try larger proteins from [test systems](test-systems.md)
- Experiment with different energy cutoffs in the dashboard
- Test with your own MD trajectories

**Found issues?**
- Report bugs in [CONTRIBUTING.md](../CONTRIBUTING.md#reporting)
- Share performance benchmarks
- Suggest improvements

**Want to develop?**
- Read the [development guide](development-guide.md)
- Check out the source code structure
- Try modifying dashboard features

---

**Total time: ~5 minutes** ⏱️  
**Questions?** Check [troubleshooting](testing-guide.md) or ask the community!
