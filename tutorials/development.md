# Development Guide

## Quick Dev Setup
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

## Common Tasks

### Fix a Bug
1. **Reproduce**: Use the exact command/data that failed
2. **Debug**: Add print statements, use debugger
3. **Fix**: Make minimal changes
4. **Test**: Verify fix works, doesn't break other things
5. **Submit**: Create pull request with clear description

### Add a Feature
1. **Check issues**: See if already requested
2. **Discuss first**: Post idea before coding
3. **Start small**: Make minimal viable version
4. **Test thoroughly**: Try different inputs
5. **Document**: Update README/docstrings

### Improve Performance
- **Profile first**: Find actual bottlenecks (`python -m cProfile`)
- **Test with large systems**: Use proteins > 200 residues
- **Measure everything**: Time, memory, disk usage
- **Compare before/after**: Quantify improvements

## Code Areas

### Core Analysis (`grinn/analysis/`)
- Energy calculations
- Network construction
- Centrality analysis

### Dashboard (`grinn/dashboard/`)
- Plotly visualizations
- Web interface
- Data formatting

### I/O (`grinn/io/`)
- File reading/writing
- Format conversions
- Error handling

## Development Tips

### Testing Your Changes
```bash
# Test with small protein first
python grinn_workflow.py test_protein.pdb results/ --top test.top --traj test.xtc

# Then test with larger system
python grinn_workflow.py large_protein.pdb results_large/ --top large.top --traj large.xtc
```

### Debugging Performance
```bash
# Profile memory usage
python -m memory_profiler grinn_workflow.py [args]

# Profile execution time
python -m cProfile -o profile.stats grinn_workflow.py [args]
```

### Code Style
- Follow existing code style
- Add docstrings for new functions
- Use meaningful variable names
- Keep functions short and focused

## Useful Development Resources
- **Python debugging**: Use `pdb` or IDE debugger
- **Performance**: `cProfile`, `memory_profiler`
- **Plotting**: Plotly documentation
- **GROMACS**: Check trajectory reading with MDAnalysis

Ready to code? Check the [main repository issues](https://github.com/osercinoglu/grinn/issues)! 💻
