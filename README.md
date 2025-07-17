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
3. **Tutorial Development**: Create comprehensive tutorials and documentation for different use cases
4. **Feature Enhancement**: Implement new features based on community feedback and needs
5. **Performance Optimization**: Improve computational efficiency and scalability

### Expected Outcomes
- Enhanced stability and reliability of gRINN
- Comprehensive documentation and tutorial resources
- Expanded test coverage with diverse protein systems
- Performance improvements and optimizations
- Community-driven feature additions

## 🏗️ How to Contribute

### For Testers
- **Test gRINN** with your own MD simulation data
- **Report bugs** and unexpected behaviors via GitHub issues
- **Share your results** and analysis experiences
- **Suggest improvements** based on your workflow needs

### For Developers
- **Fix identified bugs** and improve code quality
- **Implement new features** from the feature wishlist
- **Optimize performance** for large-scale analyses
- **Enhance the web dashboard** with new visualization options

### For Documentation Contributors
- **Write tutorials** for specific use cases (protein-ligand, membrane proteins, etc.)
- **Create example workflows** with sample data

### For Data Contributors
- **Provide test datasets** with known expected results
- **Share interesting case studies** for tutorial development

## 📋 Project Tasks

### Testing & Validation
- [ ] Test gRINN with protein-ligand complexes
- [ ] Validate results with membrane protein systems
- [ ] Test scalability with large protein complexes
- [ ] Validate Docker containerization across different platforms
- [ ] Test the interactive dashboard with various browsers

### Bug Fixes & Improvements
- [ ] Address reported issues from the main repository
- [ ] Improve error handling and user feedback
- [ ] Fix memory management issues with large trajectories
- [ ] Resolve dependency conflicts and installation issues

### Documentation & Tutorials
- [ ] Create step-by-step tutorials for common workflows
- [ ] Document best practices for different protein systems
- [ ] Develop troubleshooting guides
- [ ] Write API documentation with examples
- [ ] Create video tutorials for the web dashboard

### Feature Development
- [ ] Implement additional network analysis metrics
- [ ] Add support for different force fields
- [ ] Develop batch processing capabilities
- [ ] Add export options for network visualization tools
- [ ] Implement statistical significance testing for interactions

### Performance Optimization
- [ ] Profile and optimize computational bottlenecks
- [ ] Implement parallel processing improvements
- [ ] Optimize memory usage for large systems
- [ ] Add progress tracking for long-running calculations

## 🛠️ Getting Started

### Prerequisites
- Python 3.10+
- GROMACS (for MD simulation processing)
- Basic knowledge of molecular dynamics simulations
- Familiarity with Python and/or Docker

### Quick Setup

#### Using Docker (Recommended)
```bash
# Clone this repository
git clone https://github.com/osercinoglu/grinn-ismb-2025.git
cd grinn-ismb-2025

# Build Docker image
docker build -t grinn .

# Run with sample data
docker run -p 8051:8051 grinn dashboard test
```

#### Using Conda
```bash
# Setup environment
conda create -n grinn python=3.10
conda activate grinn

# Install dependencies
pip install -r requirements.txt

# Test installation
python -m grinn.tests.test_installation
```

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

### Sample Workflows
- Basic protein analysis
- Protein-ligand interaction studies
- Membrane protein energy networks
- Large-scale comparative analysis

### Community Resources
- [GitHub Issues](https://github.com/osercinoglu/grinn/issues)
- [Discussion Forum](https://github.com/osercinoglu/grinn/discussions)
- [CollaborationFest Slack Channel](#grinn-collaboration)

## 👥 Team & Contributors

### Project Lead
- **Onur Serçinoğlu** - Original developer and maintainer

### Developer
- **Tuğba Emine Eke** - Dashboard developer and maintainer

### Contributors
- *Add your name here when you contribute*

