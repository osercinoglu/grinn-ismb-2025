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
# First, build image with modified code.
docker build -t grinn -f Dockerfile.dev .

# Quick test with Docker
cd /path/to/test/data
docker run -v $(pwd):/data grinn-dev workflow /data/protein.gro /data/results \
  --top /data/topol.top --traj /data/traj.xtc

# Test dashboard
docker run -p 8051:8051 -v $(pwd):/data grinn-dev dashboard /data/results
```
