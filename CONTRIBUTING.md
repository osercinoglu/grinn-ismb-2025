# Contributing to gRINN CollaborationFest

## Quick Start

### Setup
1. **Install gRINN**: Follow [main repository instructions](https://github.com/osercinoglu/grinn)
2. **Install GROMACS**: You'll need this for simulations (`conda install -c conda-forge gromacs`)
3. **Fork this repo**: For sharing your results and improvements

### The Workflow (Everyone Does This)
```bash
# 1. Get a protein system
wget https://files.rcsb.org/download/1L2Y.pdb  # or use your own

# 2. Run quick simulation (see tutorials/quick-simulation.md)
python prepare_system.py 1L2Y.pdb
gmx mdrun -deffnm protein -v -nt 4  # 5-10 ns is enough

# 3. Test gRINN
python grinn_workflow.py protein.pdb results/ --top protein.top --traj protein.xtc

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
│   ├── quick-simulation.md  # 10-min simulation setup
│   ├── test-grinn.md       # Basic testing
│   ├── development.md      # For developers
│   └── documentation.md    # For writers
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
