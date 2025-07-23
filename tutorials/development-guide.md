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

## 🔄 Contributing with Fork & Pull Request Model

**Important: Do NOT clone the main repository directly for development!**

### Step 1: Fork the Repository
1. Go to https://github.com/osercinoglu/grinn
2. Click the "Fork" button in the top-right corner
3. This creates a copy of gRINN in your GitHub account

### Step 2: Clone YOUR Fork
```bash
# Clone your fork (replace YOUR_USERNAME)
git clone https://github.com/YOUR_USERNAME/grinn.git
cd grinn

# Add the original repository as 'upstream'
git remote add upstream https://github.com/osercinoglu/grinn.git

# Verify remotes are set correctly
git remote -v
# Should show:
# origin    https://github.com/YOUR_USERNAME/grinn.git (fetch)
# origin    https://github.com/YOUR_USERNAME/grinn.git (push)
# upstream  https://github.com/osercinoglu/grinn.git (fetch)
# upstream  https://github.com/osercinoglu/grinn.git (push)
```

### Step 3: Create a Feature Branch
```bash
# Always start from the latest main branch
git checkout main
git pull upstream main

# Create and switch to a new feature branch
git checkout -b fix-dashboard-memory-leak
# or: git checkout -b add-gpu-acceleration
# or: git checkout -b improve-error-handling
```

### Step 4: Make Your Changes
1. Edit the code using your preferred editor
2. Test your changes thoroughly (see Testing section below)
3. Commit your changes with clear messages:

```bash
# Stage your changes
git add .

# Commit with a descriptive message
git commit -m "Fix memory leak in dashboard when loading large protein networks

- Reduce memory usage by 40% for networks >1000 nodes
- Add garbage collection after each visualization update
- Fixes issue #123"
```

### Step 5: Push to Your Fork
```bash
# Push your feature branch to your fork
git push origin fix-dashboard-memory-leak
```

### Step 6: Create a Pull Request
1. Go to your fork on GitHub: `https://github.com/YOUR_USERNAME/grinn`
2. Click "Compare & pull request" button
3. Fill out the PR template:
   ```markdown
   ## Description
   Fix memory leak in dashboard when loading large protein networks
   
   ## Changes Made
   - Added garbage collection after visualization updates
   - Optimized network data structure for large systems
   - Reduced memory usage by ~40% for networks >1000 nodes
   
   ## Testing
   - [x] Tested with 1L2Y (small system)
   - [x] Tested with 3HTB (large system, 2000+ residues)
   - [x] Verified no regression in functionality
   - [x] Memory usage confirmed improved
   
   ## Fixes
   Closes #123
   ```
4. Click "Create pull request"

### Step 7: Respond to Review Feedback
If maintainers request changes:

```bash
# Make the requested changes
# Edit files as needed...

# Commit the updates
git add .
git commit -m "Address review feedback: improve error handling"

# Push updates to the same branch
git push origin fix-dashboard-memory-leak
```

The pull request will automatically update with your new commits!

### Step 8: Keep Your Fork Up to Date
```bash
# Regularly sync your fork with the main repository
git checkout main
git pull upstream main
git push origin main

# Update your feature branch if needed
git checkout fix-dashboard-memory-leak
git merge main  # or: git rebase main
```

## Development Workflow

### 1. Understanding the Codebase
- `grinn_workflow.py` - Main entry point and workflow orchestration
- `gRINN_Dashboard/` - Web interface components and visualization
- `mdp_files/` - GROMACS parameter files for simulations
- `data/` - Sample data for testing
- `Dockerfile.dev` - Development container configuration

### 2. Making Quality Changes
1. **Start small**: Fix one issue or add one feature per PR
2. **Test thoroughly**: Use multiple protein systems
3. **Document changes**: Update code comments and user docs
4. **Follow conventions**: Match existing code style
5. **Performance matters**: Profile memory/speed impact

### 3. Testing Your Changes
```bash
# Build development image with your changes
docker build -t grinn-dev -f Dockerfile.dev .

# Test basic functionality
docker run grinn-dev help

# Test with sample data
cd /path/to/test/data
docker run -v $(pwd):/data grinn-dev workflow /data/protein.gro /data/results \
  --top /data/topol.top --traj /data/traj.xtc

# Test dashboard
docker run -p 8051:8051 -v $(pwd):/data grinn-dev dashboard /data/results

# Visit http://localhost:8051 to verify dashboard works
```

### 4. Submitting Quality Pull Requests

#### What Makes a Good PR?
- **Clear purpose**: Fix one issue or add one feature
- **Thorough testing**: Works with multiple protein systems
- **Good documentation**: Clear commit messages and PR description
- **No regressions**: Doesn't break existing functionality
- **Performance awareness**: Consider memory/speed impact

#### PR Template Example
```markdown
## Description
Brief summary of what this PR accomplishes

## Type of Change
- [ ] Bug fix (fixes an issue)
- [ ] New feature (adds functionality)
- [ ] Performance improvement
- [ ] Documentation update

## Changes Made
- Specific change 1
- Specific change 2
- Specific change 3

## Testing
- [ ] Tested with small protein (< 100 residues)
- [ ] Tested with large protein (> 500 residues)
- [ ] Dashboard loads and functions correctly
- [ ] No performance regression observed
- [ ] Added/updated tests if applicable

## Evidence
Paste screenshots, output logs, or performance metrics showing the fix works

## Issues
Fixes #(issue_number)
```

### 5. After Your PR is Merged
```bash
# Switch back to main branch
git checkout main

# Pull the latest changes (including your merged PR)
git pull upstream main

# Update your fork's main branch
git push origin main

# Clean up your feature branch
git branch -d fix-dashboard-memory-leak
git push origin --delete fix-dashboard-memory-leak
```

## 🎯 Development Best Practices

### Git Workflow
- **One feature per branch**: Keep PRs focused and small
- **Descriptive branch names**: `fix-memory-leak`, `add-gpu-support`
- **Clear commit messages**: Explain what and why, not just what
- **Regular syncing**: Keep your fork up to date with upstream

### Code Quality
- **Test before committing**: Always verify your changes work
- **Follow existing style**: Match the codebase conventions
- **Add comments**: Explain complex logic
- **Handle errors gracefully**: Don't let the tool crash on edge cases

### Performance Considerations
- **Memory efficiency**: Profile memory usage with large systems
- **Speed optimization**: Measure runtime impact
- **Scalability**: Test with various protein sizes
- **Resource cleanup**: Ensure proper cleanup of temporary files
