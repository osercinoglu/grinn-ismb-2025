# Contributing to gRINN CollaborationFest

## Quick Start

### Setup
1. **Install Docker**: Make sure Docker is installed and running
2. **Build gRINN Docker image**: 
   ```bash
   git clone https://github.com/osercinoglu/grinn.git
   cd grinn
   docker build -t grinn .
   ```
3. **Fork this repo**: For sharing your results and improvements

### The Workflow (Everyone Does This)
```bash
# 1. Create working directory and get a protein system
mkdir grinn_test && cd grinn_test
wget https://files.rcsb.org/download/1L2Y.pdb  # or use your own

# 2. Run quickstart workflow using Docker (see tutorials/quickstart.md)
docker build -f Dockerfile.dev -t grinn-dev . && cd path/to/test/data
# ... follow full workflow in quickstart.md

# 3. Test gRINN using Docker
docker run -v $(pwd):/data grinn workflow /data/em.gro /data/results --top /data/topol.top --traj /data/protein_traj.xtc

# 4. Check results and report back!
```

## What to Report

### ✅ Success
- "Tested protein X, worked fine, took Y minutes"
- Interesting findings or networks
- Performance observations

### 🐛 Problems
- Error messages (copy-paste the full error)
- What you were doing when it failed
- System details (protein size, OS, etc.)

### 💡 Ideas
- Feature suggestions
- Performance improvements
- Documentation clarifications

## Reporting Results

### GitHub Issues
- **Bugs**: Use bug report template
- **Features**: Use feature request template
- **Questions**: Use discussion forum

### Share Your Work
- Fork this repository
- Add your results to `results/` folder
- Create pull request with summary

## Project Structure
```
grinn-ismb-2025/
├── README.md           # Start here
├── CONTRIBUTING.md     # This file
├── tutorials/          # Short, focused guides
│   ├── test-systems.md      # Protein suggestions
│   ├── quickstart.md        # Complete 5-min workflow
│   ├── test-systems.md      # Recommended proteins
│   ├── development-guide.md # For developers
│   └── documentation-guide.md # For writers
└── results/            # Share your findings here
```

## Getting Help

- **Setup issues**: Ask in GitHub discussions
- **GROMACS help**: Check GROMACS documentation
- **gRINN bugs**: Report in main repository
- **During event**: Ask organizers

## Tips

- **Start small**: Use small proteins first (< 100 residues)
- **Short simulations**: 5-10 ns is plenty for testing
- **Document everything**: Note what works and what doesn't
- **Ask questions**: We're here to help!

Ready to contribute? Start with the [test systems guide](tutorials/test-systems.md)! 🚀
