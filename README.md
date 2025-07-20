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
1. **[Setup](CONTRIBUTING.md#setup)** - Install gRINN and GROMACS
2. **[Get Test Data](tutorials/test-systems.md)** - Download a protein or use your own
3. **[Quick Simulation](tutorials/quick-simulation.md)** - Generate test trajectory (10-15 min)
4. **[Test gRINN](tutorials/test-grinn.md)** - Run the analysis
5. **[Report Results](CONTRIBUTING.md#reporting)** - Share what you found

### Skip Steps If:
- **Have existing trajectories?** Skip to Step 4 (testing)
- **Want to develop?** First test the tool, then read [development notes](tutorials/development.md)
- **Want to write tutorials?** First test successfully, then check [writing guide](tutorials/documentation.md)

## 📊 Current Status

### What's Working
- ✅ Core interaction energy calculations
- ✅ Protein Energy Network construction
- ✅ Interactive web dashboard
- ✅ Docker containerization
- ✅ Basic documentation

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

## 👥 Team & Contributors

### Project Lead
- **Onur Serçinoğlu** - Original developer and maintainer

### Developer
- **Tuğba Emine Eke** - Dashboard developer and maintainer

### Contributors
- *Add your name here when you contribute*

