# Quick Simulation Setup

## Goal: Generate a test trajectory in 10-15 minutes

### Prerequisites
- GROMACS installed (`conda install -c conda-forge gromacs`)
- A PDB file (from [test systems](test-systems.md))

### Step-by-Step
```bash
# 1. Download test protein
wget https://files.rcsb.org/download/1L2Y.pdb

# 2. Prepare system
gmx pdb2gmx -f 1L2Y.pdb -o protein.gro -p topol.top -ff amber99sb-ildn -water none
gmx editconf -f protein.gro -o boxed.gro -c -d 1.0 -bt cubic

# 3. Energy minimize (2 min)
gmx grompp -f minim.mdp -c boxed.gro -p topol.top -o em.tpr
gmx mdrun -deffnm em -v

# 4. Quick NVT (5 min simulation)
gmx grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt -v -nt 4

# 5. Extract trajectory for gRINN
echo "Protein" | gmx trjconv -s nvt.tpr -f nvt.xtc -o protein_traj.xtc -pbc mol
```

### MDP Files You Need

Create these simple .mdp files:

**minim.mdp** (energy minimization):
```
integrator = steep
nsteps = 5000
emtol = 1000.0
```

**nvt.mdp** (short simulation):
```
integrator = md
dt = 0.002
nsteps = 2500000    ; 5 ns
tcoupl = V-rescale
tc-grps = Protein
tau-t = 0.1
ref-t = 300
```

### Expected Runtime
- **Small protein** (< 50 residues): 5-10 minutes total
- **Medium protein** (100 residues): 10-15 minutes total
- **Large protein** (> 200 residues): 20+ minutes

### Troubleshooting
- **"Fatal error"**: Check your GROMACS installation
- **"Blowing up"**: Try shorter timestep (dt = 0.001)
- **Too slow**: Use fewer steps or smaller system

### Output Files for gRINN
You need:
- `protein.gro` or `protein.pdb` (structure)
- `topol.top` (topology)
- `protein_traj.xtc` (trajectory)

**Next**: Test these files with gRINN → [test-grinn.md](test-grinn.md) 🚀
