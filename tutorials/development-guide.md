# Development Guide - gRINN CollaborationFest

This guide will help developers contribute code improvements, bug fixes, and new features to gRINN.

## Overview

As a developer, you can contribute by:
- Fixing reported bugs and issues
- Improving code quality and performance
- Enhancing the web dashboard
- Optimizing algorithms for better scalability

## Prerequisites

- **Python 3.10+** with scientific computing experience
- **Git** for version control
- **Understanding of molecular dynamics** concepts (helpful)
- **Web development skills** (for dashboard work)
- **Docker knowledge** (optional)

## Setting Up Development Environment

### 1. Fork and Clone
```bash
# Fork the repository on GitHub, then clone your fork
git clone https://github.com/YOUR_USERNAME/grinn.git
cd grinn

# Add upstream remote
git remote add upstream https://github.com/osercinoglu/grinn.git
```

### 2. Create Development Environment
```bash
# Using conda (recommended)
conda env create -f environment.yml
conda activate grinn

# Or using pip
pip install -r requirements.txt
```

### 3. Verify Setup
```bash
# Run basic tests
python -c "import prody, numpy, scipy, pandas; print('Dependencies OK')"

# Test basic functionality
python grinn_workflow.py --help
```

## Code Structure Overview

### Key Files and Directories
```
grinn/
├── grinn_workflow.py          # Main workflow orchestration
├── gRINN_Dashboard/           # Web dashboard application
│   ├── grinn_dashboard.py     # Main dashboard script
│   ├── assets/                # CSS, JS, images
│   └── components/            # Dashboard components
├── test_data/                 # Sample datasets
├── environment.yml            # Conda environment
├── Dockerfile                 # Container configuration
└── README.md                  # Documentation
```

### Core Components

#### Main Workflow (`grinn_workflow.py`)
- **Entry point**: Command-line interface and argument parsing
- **GROMACS integration**: MD simulation processing
- **Energy calculations**: Core interaction energy algorithms
- **Network construction**: PEN building and analysis
- **File I/O**: Input/output handling

#### Dashboard (`gRINN_Dashboard/`)
- **Plotly Dash application**: Web interface
- **Data visualization**: Interactive plots and networks
- **File upload/management**: Result loading and processing
- **3D molecular viewer**: Protein structure display

## Development Workflow

### 1. Choose an Issue
- Browse [open issues](https://github.com/osercinoglu/grinn/issues)
- Look for labels: `good first issue`, `bug`, `enhancement`
- Comment on the issue to indicate you're working on it

### 2. Create a Branch
```bash
# Create feature branch
git checkout -b feature/your-feature-name

# Or for bug fixes
git checkout -b bugfix/issue-description
```

### 3. Make Changes
- Write clean, readable code
- Follow existing code style
- Add comments for complex logic
- Update documentation as needed

### 4. Test Your Changes
```bash
# Run with test data
python grinn_workflow.py test_data/sample.pdb output/ \
  --top test_data/sample.top \
  --traj test_data/sample.xtc

# Test dashboard (if modified)
python gRINN_Dashboard/grinn_dashboard.py output/
```

### 5. Commit and Push
```bash
# Stage changes
git add .

# Commit with descriptive message
git commit -m "Fix memory leak in energy calculation loop

- Properly close file handles in parse_interaction_energies
- Add garbage collection after processing large chunks
- Resolves issue #123"

# Push to your fork
git push origin feature/your-feature-name
```

### 6. Create Pull Request
- Go to GitHub and create a PR from your branch
- Fill out the PR template with details
- Link to related issues
- Request review from maintainers

## Common Development Tasks

### Bug Fixes

#### Memory Management Issues
```python
# Before: Memory not properly managed
def process_trajectory(trajectory_file):
    data = load_large_trajectory(trajectory_file)
    # ... processing ...
    return results

# After: Proper memory management
def process_trajectory(trajectory_file):
    try:
        data = load_large_trajectory(trajectory_file)
        # ... processing ...
        return results
    finally:
        if 'data' in locals():
            del data
        gc.collect()
```

#### Error Handling Improvements
```python
# Before: Generic error handling
try:
    result = calculate_energies(pairs)
except Exception as e:
    print(f"Error: {e}")

# After: Specific error handling
try:
    result = calculate_energies(pairs)
except FileNotFoundError:
    logger.error("Input files not found. Check file paths.")
    raise
except MemoryError:
    logger.error("Insufficient memory. Try using --skip parameter.")
    raise
except Exception as e:
    logger.error(f"Unexpected error in energy calculation: {e}")
    raise
```

### Performance Optimization

#### Parallel Processing
```python
# Example: Optimizing interaction calculations
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp

def optimize_energy_calculation(pairs_chunks, num_cores):
    with ProcessPoolExecutor(max_workers=num_cores) as executor:
        futures = [executor.submit(calculate_chunk, chunk) 
                  for chunk in pairs_chunks]
        results = [future.result() for future in futures]
    return results
```

#### Memory Optimization
```python
# Use generators for large datasets
def process_trajectory_frames(trajectory):
    for frame_idx, frame in enumerate(trajectory):
        yield process_frame(frame)
        # Memory is freed after each iteration

# Instead of loading everything at once
results = [process_frame(frame) for frame in trajectory]  # High memory
results = list(process_trajectory_frames(trajectory))     # Lower memory
```

### Dashboard Enhancements

#### Adding New Visualizations
```python
# In gRINN_Dashboard/grinn_dashboard.py

@app.callback(
    Output('new-plot', 'figure'),
    [Input('data-store', 'data')]
)
def update_new_plot(data):
    if not data:
        return {}
    
    df = pd.DataFrame(data)
    
    fig = px.scatter(df, x='x_col', y='y_col',
                     title='New Visualization')
    
    fig.update_layout(
        template='plotly_white',
        showlegend=True
    )
    
    return fig
```

#### Improving User Interface
```python
# Add new dashboard components
new_component = html.Div([
    html.H3("New Feature"),
    dcc.Graph(id='new-plot'),
    dbc.Button("Export Data", id="export-btn", 
               color="primary", className="mt-2")
])
```

## Testing Your Code

### Unit Testing
```python
# Create test files in tests/ directory
import unittest
from grinn_workflow import calculate_interaction_energies

class TestEnergyCalculations(unittest.TestCase):
    def test_basic_calculation(self):
        # Test with known input/output
        pairs = [(0, 1), (1, 2)]
        result = calculate_interaction_energies(pairs, test_data)
        self.assertIsNotNone(result)
        self.assertEqual(len(result), 2)

if __name__ == '__main__':
    unittest.main()
```

### Integration Testing
```bash
# Test complete workflows
python grinn_workflow.py test_data/small_protein.pdb test_output/ \
  --top test_data/small_protein.top \
  --traj test_data/small_protein.xtc \
  --test_only

# Check outputs exist and are valid
ls test_output/
head test_output/energies_intEnTotal.csv
```

## Code Style Guidelines

### Python Style
- Follow PEP 8 conventions
- Use meaningful variable names
- Add docstrings to functions
- Keep functions focused and small

```python
def calculate_residue_distances(structure, residue_indices, cutoff=10.0):
    """
    Calculate distances between specified residues.
    
    Parameters
    ----------
    structure : prody.AtomGroup
        Protein structure
    residue_indices : list
        List of residue indices to analyze
    cutoff : float, default 10.0
        Distance cutoff in Angstroms
        
    Returns
    -------
    dict
        Dictionary mapping residue pairs to distances
    """
    distances = {}
    # Implementation here...
    return distances
```

### Git Commit Messages
- Use present tense ("Add feature" not "Added feature")
- Keep first line under 50 characters
- Reference issue numbers when applicable

```
Fix memory leak in trajectory processing

- Close file handles properly in MDAnalysis reader
- Add explicit garbage collection for large frames
- Reduce memory usage by 30% for trajectories >1000 frames

Fixes #142
```

## Debugging Tips

### Common Issues

#### GROMACS Integration
```python
# Debug GROMACS command execution
import subprocess
import logging

def run_gromacs_command(cmd):
    logging.info(f"Running: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logging.error(f"GROMACS error: {result.stderr}")
        return result
    except FileNotFoundError:
        logging.error("GROMACS not found in PATH")
        raise
```

#### Dashboard Issues
```python
# Debug dashboard data flow
import json

@app.callback(...)
def debug_callback(input_data):
    print(f"Input data: {json.dumps(input_data, indent=2)}")
    # Your callback logic here
    result = process_data(input_data)
    print(f"Output data: {json.dumps(result, indent=2)}")
    return result
```

### Performance Profiling
```python
import cProfile
import pstats

# Profile your code
cProfile.run('your_function()', 'profile_output')

# Analyze results
stats = pstats.Stats('profile_output')
stats.sort_stats('cumulative').print_stats(10)
```

## Contributing Guidelines

### Before Submitting
- [ ] Code follows style guidelines
- [ ] New features have tests
- [ ] Documentation updated if needed
- [ ] Changes work with sample data
- [ ] No new dependencies without discussion

### Pull Request Template
```markdown
## Description
Brief description of changes

## Type of Change
- [ ] Bug fix
- [ ] New feature
- [ ] Performance improvement
- [ ] Documentation update

## Testing
- [ ] Tested with sample data
- [ ] Added/updated tests
- [ ] Manual testing completed

## Related Issues
Closes #123
```

## Advanced Topics

### Adding New Force Field Support
1. Identify force field parameters
2. Update topology parsing logic
3. Add parameter validation
4. Test with force field-specific data

### Implementing New Network Metrics
1. Research network analysis algorithms
2. Implement efficient calculation
3. Add to dashboard visualization
4. Validate with known test cases

### Docker Development
```dockerfile
# For testing Docker changes
FROM python:3.10-slim

# Install dependencies
COPY requirements.txt .
RUN pip install -r requirements.txt

# Copy source code
COPY . /app
WORKDIR /app

# Test build
RUN python -c "import grinn_workflow; print('Build successful')"
```

## Resources

- [Python Style Guide (PEP 8)](https://pep8.org/)
- [Plotly Dash Documentation](https://dash.plotly.com/)
- [ProDy Documentation](http://prody.csb.pitt.edu/)
- [GROMACS Manual](https://manual.gromacs.org/)
- [NetworkX Documentation](https://networkx.org/)

## Getting Help

- **GitHub Issues**: For bugs and feature requests
- **GitHub Discussions**: For development questions
- **Code Review**: Request reviews from maintainers
- **CollaborationFest**: Ask coordinators during the event

Ready to start developing? Pick an issue and start coding! 💻
