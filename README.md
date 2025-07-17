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
4. **Feature Enhancement**: Implement new features based on community feedback and needs
5. **Performance Optimization**: Improve computational efficiency and scalability

## 🏗️ How to Contribute

### For Testers
- **Test gRINN** with your own MD simulation data or using the sample input trajectories that will be provided.
- **Report bugs** and unexpected behaviors via GitHub issues
- **Share your results** and analysis experiences
- **Suggest improvements** based on your workflow needs

### For Developers
- **Fix identified bugs** and improve code quality
- **Implement new features** from the feature wishlist
- **Optimize performance** for large-scale analyses
- **Enhance the web dashboard** with new visualization options

### For Documentation Contributors
- **Write tutorials** for specific use cases
- **Create example workflows** with sample data

## 🛠️ Getting Started

For detailed installation instructions and setup guides, please visit the [main gRINN repository](https://github.com/osercinoglu/grinn).

### Sample Data
We will provide sample datasets for testing:
- Small protein system (< 100 residues)
- Protein-ligand complex
- Membrane protein system
- Large protein complex (> 500 residues)

Contributors may also bring their own GROMACS-generated simulation trajectories.

## 📊 Current Status

### What's Working
- ✅ Core interaction energy calculations
- ✅ Protein Energy Network construction
- ✅ Interactive web dashboard
- ✅ Docker containerization
- ✅ Basic documentation (as found in this README)

### Known Issues
- 🔍 Memory usage with very large trajectories
- 🔍 Performance optimization needed for >1000 residue systems
- 🔍 Limited documentation for advanced use cases
- 🔍 Dashboard compatibility with some browser versions

## 📚 Resources

### Documentation
- [Main gRINN Repository](https://github.com/osercinoglu/grinn)
- [Publication (of the first version of the tool)](https://doi.org/10.1093/nar/gky381)
- [Installation Guide](https://github.com/osercinoglu/grinn#installation--usage)

## 👥 Team & Contributors

### Project Lead
- **Onur Serçinoğlu** - Original developer and maintainer

### Developer
- **Tuğba Emine Eke** - Dashboard developer and maintainer

### Contributors
- *Add your name here when you contribute*

