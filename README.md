# gRINN CollaborationFest - Test & Improve Together! 🚀

**ISMB/ECCB 2025 Project: Making protein energy network analysis better for everyone**

## What You'll Do (Pick Your Speed)

### ⚡ Quick Test (5 minutes)
```bash
# Get gRINN
git clone https://github.com/osercinoglu/grinn.git
cd grinn && docker build -t grinn .

# Test it
mkdir results
docker run -v $(pwd)/results:/results grinn test

# View results
docker run -p 8051:8051 -v $(pwd)/results:/results grinn dashboard /results
# Open: http://localhost:8051
```

### 🧪 Real Protein Test (15-30 minutes)
```bash
# Get a protein (tiny one for speed)
wget https://files.rcsb.org/download/1L2Y.pdb

# Run full analysis
mkdir my_results
docker run -v $(pwd):/data -v $(pwd)/my_results:/results grinn \
  workflow /data/1L2Y.pdb /results --gpu --solvate

# View interactive results
docker run -p 8051:8051 -v $(pwd)/my_results:/results grinn dashboard /results
```

### 🎯 Your Own Data Test
```bash
# Use your protein/trajectory
docker run -v $(pwd):/data -v $(pwd)/results:/results grinn \
  workflow /data/your_protein.pdb /results --top /data/your_topology.top --traj /data/your_trajectory.xtc
```

## What to Look For

### ✅ Success = Report It!
- Reasonable energy values (-10 to +5 kcal/mol)
- Dashboard loads and is interactive
- Strong interactions make chemical sense
- Networks show expected connectivity

### 🐛 Problems = Report Them Too!
- Build failures
- Crashes or hangs
- Weird energy values
- Dashboard errors
- Performance issues

## Report Results

**Found bugs?** [File issue here](https://github.com/osercinoglu/grinn/issues/new)

**Everything worked?** [Share success here](https://github.com/osercinoglu/grinn-ismb-2025/issues/new/choose)

**Need help?** Ask during CollaborationFest or in [discussions](https://github.com/osercinoglu/grinn-ismb-2025/discussions)

## Test Proteins (Copy-Paste Ready)

```bash
# Tiny (< 1 min) - Good for debugging
wget https://files.rcsb.org/download/1L2Y.pdb  # Trp-cage, 20 residues

# Small (< 5 min) - Good for validation  
wget https://files.rcsb.org/download/1UBQ.pdb  # Ubiquitin, 76 residues

# Medium (< 30 min) - Good for stress testing
wget https://files.rcsb.org/download/1LYS.pdb  # Lysozyme, 129 residues
```

## Advanced Testing (Optional)

### Performance Testing
```bash
# Test memory usage
docker stats --no-stream grinn_container_name

# Test large protein
wget https://files.rcsb.org/download/1IGY.pdb  # 450 residues
time docker run -v $(pwd):/data -v $(pwd)/results_large:/results grinn \
  workflow /data/1IGY.pdb /results --gpu --skip 10
```

### Custom Analysis
```bash
# Test specific parameters
docker run -v $(pwd):/data -v $(pwd)/results:/results grinn \
  workflow /data/protein.pdb /results --initpairfiltercutoff 5.0 --create_pen --pen_cutoffs 1.0 2.0
```

### Scientific Validation
- Compare results with known protein interactions
- Validate against experimental binding data
- Cross-check with other network analysis tools

## That's It!

**No complex workflows, no tutorial maze, no setup confusion.**

Just test, report, and help make gRINN better! 🎉

---

### Team
- **Onur Serçinoğlu** - Original developer 
- **Tuğba Emine Eke** - Dashboard developer
- **You?** - Add your name when you contribute!

