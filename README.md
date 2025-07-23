# gRINN Collaborative Development - ISMB/ECCB 2025 CollaborationFest

**Advancing Protein Energy Network Analysis through Community Contributions**

---

## 🎯 Project Overview

This CollaborationFest project focuses on enhancing and expanding [**gRINN**](https://github.com/osercinoglu/grinn) (get Residue Interaction eNergies and Networks), a tool for calculating residue interaction energies and constructing protein energy networks from molecular dynamics simulations. We invite bioinformaticians, computational biologists, and MD simulation enthusiasts to join us in testing, improving, and documenting this resource for the structural biology community.

### About gRINN
gRINN is a comprehensive Python-based tool that:
- Calculates pairwise amino acid non-bonded interaction energies from GROMACS MD trajectories
- Constructs Protein Energy Networks (PENs) with customizable energy cutoffs
- Performs betweenness centrality analysis to identify key residues in protein communication pathways
- Provides an interactive web dashboard for visualization and analysis
- Supports Docker containerization for reproducible analysis

## 🚀 Goals

### Primary Objectives
1. **Tool Testing & Validation**: Test gRINN with diverse protein systems and simulation datasets
2. **Bug Identification & Fixes**: Identify and resolve issues encountered during real-world usage
3. **Tutorial Development**: Create tutorials for different use cases

## 🛠️ How to Contribute

### Everyone Starts Here:
1. **Get a protein system** - Use your own data OR pick from our [suggested systems](tutorials/test-systems.md)
2. **Run a quick simulation** - Generate a short trajectory (or use existing data)
3. **Test gRINN** - Run the tool and see what happens
4. **Report back** - Share results, bugs, or suggestions

### Then Choose Your Path:
- **🧪 Continue Testing**: Try more proteins, report bugs, validate edge cases
- **💻 Develop**: Fix issues you found, add features, improve performance
- **📖 Write Tutorials**: *After successful testing* - create guides for your specific use case

## 🛠️ Getting Started

### Step-by-Step Workflow
1. **[Setup](CONTRIBUTING.md#setup)** - Build gRINN Docker image
2. **[Get Test Data](tutorials/test-systems.md)** - Download a protein or use your own
3. **[Testing Guide](tutorials/testing-guide.md)** - Complete workflow from protein to results
4. **[Test Systems](tutorials/test-systems.md)** - Recommended proteins for testing
5. **[Report Results](CONTRIBUTING.md#reporting)** - Share what you found

### Start Here

#### ⚡ Setup
```bash
# Get gRINN (docker build will take time, so be patient...)
git clone https://github.com/osercinoglu/grinn.git
cd grinn && docker build -t grinn -f Dockerfile.dev .
```

#### 🧪 Real Protein Test with GROMACS Preparation (15-30 minutes)
```bash
# Get a protein (tiny one for speed)
mkdir test_1L2Y
cd test_1L2Y
wget https://files.rcsb.org/download/1L2Y.pdb

# Run all GROMACS steps inside container to avoid permission issues
docker run -it --rm -v $(pwd):/host_data grinn bash

# Then run the following commands:
cd /tmp
cp /host_data/1L2Y.pdb .

# Step 1: Prepare protein structure with AMBER force field
gmx pdb2gmx -f 1L2Y.pdb -o processed.gro -p topol.top -ff amber99sb-ildn -water tip3p -ignh

# Step 2: Define simulation box
gmx editconf -f processed.gro -o boxed.gro -c -d 1.0 -bt cubic

# Step 3: Add solvent
gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top

# Step 4: Energy minimization (adjust -nt for your device)
gmx grompp -f /app/mdp_files/minim.mdp -c solvated.gro -p topol.top -o em.tpr -maxwarn 5
gmx mdrun -v -deffnm em -nt 8

# Step 5: Short (toy) MD simulation
gmx grompp -f /app/mdp_files/prod.mdp -c em.gro -p topol.top -o md.tpr
gmx mdrun -v -deffnm md

# Copy results back to host
cp *.gro *.top *.tpr *.xtc *.edr *.log /host_data/
exit

# Step 6: Run gRINN analysis

docker run -w /data -v $(pwd):/data grinn workflow /data/em.gro /data/grinn_results --top /data/topol.top --traj /data/md.xtc

# Step 7: View interactive results
docker run -p 8051:8051 -v $(pwd):/data grinn dashboard /data/grinn_results
```

#### 🎯 Your Own Data Test
```bash

# Use your protein/trajectory
docker run -w /data -v $(pwd):/data grinn workflow /data/your_protein.pdb /data/results --top /data/your_system.top --traj /data/your_trajectory.xtc
```

### Skip Steps If:
- **Have existing trajectories?** Skip to Step 6 (gRINN analysis)
- **Want to develop?** First test the tool using Docker, then read [development notes](tutorials/development-guide.md)
- **Want to write tutorials?** First test successfully using Docker, then check [writing guide](tutorials/documentation-guide.md)

## Report Results

**Found bugs?** [File issue here](https://github.com/osercinoglu/grinn/issues/new)

## 📊 Current Status

### What's Working
- ✅ Core interaction energy calculations
- ✅ Protein Energy Network construction
- ✅ Interactive web dashboard
- ✅ Docker containerization
- ✅ Basic documentation
- ✅ GRO file support with proper CRYST1 record preservation (box information)

### Areas for Improvement
- 🎯 **Testing**: Validate with diverse protein systems and simulation setups
- 🎯 **Performance**: Optimize memory usage and speed for large systems
- 🎯 **Documentation**: Create comprehensive tutorials and troubleshooting guides
- 🎯 **Robustness**: Fix edge cases and improve error handling
- 🎯 **Usability**: Enhance dashboard features and user experience

## 📚 Resources & Links

### Essential Resources
- **[Main gRINN Repository](https://github.com/osercinoglu/grinn)** - Installation, source code, and basic documentation
- **[Scientific Publication](https://doi.org/10.1093/nar/gky381)** - Original method description and validation
- **[GROMACS Tutorials](http://www.mdtutorials.com/gmx/)** - Learn MD simulation basics if needed


### Custom Analysis
```bash

# Test specific parameters
docker run -w /data -v $(pwd):/data -v $(pwd)/results:/results grinn \
  workflow /data/protein.pdb /results --initpairfiltercutoff 5.0 --create_pen --pen_cutoffs 1.0 2.0
```

### Scientific Validation
- Compare results with known protein interactions
- Validate against experimental binding data
- Cross-check with other network analysis tools

## 👥 Team & Contributors

### Project Lead
- **Onur Serçinoğlu** - Original developer and maintainer

### Developer
- **Tuğba Emine Eke** - Dashboard developer and maintainer

### Contributors
- *Add your name here when you contribute*

