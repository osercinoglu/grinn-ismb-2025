# Development Guide for gRINN

## Quick Setup

### Option 1: Docker Development (Recommended)
```bash
# 1. Clone and build Docker image
git clone https://github.com/osercinoglu/grinn.git
cd grinn
docker build -f Dockerfile.dev -t grinn-dev .

# 2. Test Docker setup
docker run grinn-dev help

# 3. Run interactive session for development
docker run -it -v $(pwd):/app grinn-dev bash
```

### Option 2: Local Installation
```bash
# 1. Clone and setup
git clone https://github.com/osercinoglu/grinn.git
cd grinn
pip install -e .  # Install in development mode

# 2. Test installation
python -c "import grinn; print('Dev setup OK!')"

# 3. Run tests (if available)
python -m pytest tests/
```

## Prerequisites
- Docker installed and running
- Git for version control
- Text editor or IDE
- Basic knowledge of Python and molecular dynamics

## Setup Development Environment

### Docker-based Development (Recommended)
```bash
# 1. Clone the repository
git clone https://github.com/osercinoglu/grinn.git
cd grinn

# 2. Build development Docker image
docker build -f Dockerfile.dev -t grinn-dev .

# 3. Start interactive development session
docker run -it -v $(pwd):/app -w /app grinn-dev bash

# 4. Inside container, you can edit files and test immediately
```

### Local Development Setup
```bash
# 1. Clone repository
git clone https://github.com/osercinoglu/grinn.git
cd grinn

# 2. Create conda environment
conda create -n grinn-dev python=3.10
conda activate grinn-dev

# 3. Install dependencies
pip install -e .
conda install -c conda-forge gromacs

# 4. Test installation
python -c "import grinn; print('Setup complete!')"
```

## Development Workflow

### 1. Understanding the Codebase
- `grinn_workflow.py` - Main entry point
- `gRINN_Dashboard/` - Web interface components
- `mdp_files/` - GROMACS parameter files

### 2. Making Changes
1. Create a feature branch: `git checkout -b feature-name`
2. Make your changes
3. Test with small systems first
4. Test with larger systems
5. Update documentation if needed

### 3. Testing Your Changes
```bash
# Quick test with Docker
cd /path/to/test/data
docker run -v $(pwd):/data grinn-dev workflow /data/protein.gro /data/results \
  --top /data/topol.top --traj /data/traj.xtc

# Test dashboard
docker run -p 8051:8051 -v $(pwd):/data grinn-dev dashboard /data/results
```

## Common Development Tasks

### Adding New Features
1. **Plan**: Describe the feature in an issue first
2. **Implement**: Start with minimal working version
3. **Test**: Use test systems from tutorials/test-systems.md
4. **Document**: Update relevant documentation
5. **Submit**: Create pull request with clear description

### Fixing Bugs
1. **Reproduce**: Use the exact scenario that triggers the bug
2. **Debug**: Add logging or use debugger
3. **Fix**: Implement minimal fix
4. **Test**: Verify fix works and doesn't break other functionality
5. **Submit**: Create pull request referencing the issue

### Performance Optimization
1. **Profile**: Use tools like `cProfile` or `memory_profiler`
2. **Identify bottlenecks**: Focus on actual performance issues
3. **Optimize**: Make targeted improvements
4. **Benchmark**: Compare before/after performance
5. **Test**: Ensure optimization doesn't break functionality

## Code Style Guidelines
- Follow PEP 8 for Python code style
- Use meaningful variable and function names
- Add docstrings for new functions and classes
- Keep functions focused and relatively short
- Comment complex algorithms

## Testing Guidelines
- Test with small proteins first (< 50 residues)
- Use multiple test systems with different characteristics
- Test edge cases (very small/large proteins, unusual structures)
- Verify output format consistency (.gro files, proper trajectory format)
- Test both Docker and local installations when possible

## Debugging Tips
- Use Docker's interactive mode for debugging: `docker run -it grinn-dev bash`
- Check file permissions when using volume mounts
- Verify GROMACS installation inside container
- Use print statements or logging for tracking execution flow

## Performance Considerations
- Large proteins (>500 residues) may require significant memory
- Consider frame skipping for very long trajectories
- Profile memory usage for large systems
- Test on different hardware configurations when possible

## Getting Help
- Check existing issues on GitHub
- Ask questions in discussions
- Reach out to project maintainers
- Use community resources for GROMACS/MDAnalysis questions

## Resources
- [Main gRINN Repository](https://github.com/osercinoglu/grinn)
- [GROMACS Documentation](https://manual.gromacs.org/)
- [MDAnalysis Documentation](https://docs.mdanalysis.org/)
- [Plotly/Dash Documentation](https://dash.plotly.com/)

Ready to contribute? Start by testing the current version, then pick an issue to work on! 🚀
