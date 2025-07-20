# Testing gRINN

## Basic Test Run

### Command
```bash
python grinn_workflow.py protein.gro results/ \
  --top topol.top \
  --traj protein_traj.xtc \
  --nt 4
```

### What Should Happen
1. **Processing**: Lots of output text showing energy calculations
2. **Files created**: Check `results/` folder for output files
3. **Dashboard**: Opens web browser with interactive plots
4. **Runtime**: Few minutes for small proteins, longer for big ones

### Check Your Results

**Files to look for:**
- `results/interaction_energies.csv` - Energy data
- `results/network_data.json` - Network structure
- `results/dashboard.html` - Interactive plots

**Quick validation:**
```bash
# Count interactions found
wc -l results/interaction_energies.csv

# Check for reasonable energy values
head results/interaction_energies.csv
```

### What to Report

#### ✅ Success
- "Worked! Protein X, Y interactions found, took Z minutes"
- Share interesting network features
- Note performance (speed, memory usage)

#### 🐛 Problems
Post the **full error message** plus:
- Protein used (PDB ID, size)
- Command you ran
- Your system (OS, Python version)
- Input file details

#### 🤔 Weird Results
- Unexpected energy values
- Strange network topology
- Dashboard not working
- Performance issues

### Common Issues

**"ModuleNotFoundError"**: gRINN not installed properly
**"File not found"**: Check file paths are correct
**"Memory error"**: Try smaller protein or skip frames (`--skip 10`)
**"Dashboard won't open"**: Try different browser

### Advanced Testing

Once basic test works:
```bash
# Test different cutoffs
python grinn_workflow.py protein.gro results_strict/ --top topol.top --traj protein_traj.xtc --cutoff -10

# Test frame skipping
python grinn_workflow.py protein.gro results_fast/ --top topol.top --traj protein_traj.xtc --skip 5

# Test specific residue ranges
python grinn_workflow.py protein.gro results_core/ --top topol.top --traj protein_traj.xtc --source_sel "resid 1-20"
```

### Share Your Results
1. Fork this repository
2. Add your results to `results/your_protein_name/`
3. Create pull request with summary
4. Or just post in GitHub discussions!

**Tip**: Even failed tests are valuable - report everything! 🧪
