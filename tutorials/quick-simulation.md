# Quick Simulation Setup

## Goal: Generate a minimal test trajectory in 3-5 minutes using gRINN Docker

This tutorial covers a streamlined MD workflow: energy minimization followed by a very short production simulation (just 1000 steps) for quick gRINN testing.

### Prerequisites
- Docker installed and running
- gRINN Docker image built:
  - **For development/testing**: `docker build -f Dockerfile.dev -t grinn-dev .` (faster)
  - **For production**: `docker build -t grinn .` (with tests)
- A PDB file (from [test systems](test-systems.md))

**Note**: Use `Dockerfile.dev` for faster builds during development. It skips tests and builds quicker for iteration.

### Step-by-Step
```bash
# 1. Create a working directory and download test protein
mkdir -p grinn_sim && cd grinn_sim
wget https://files.rcsb.org/download/1L2Y.pdb

# 2. Prepare system using Docker gmx mode
docker run -v $(pwd):/data -w /data grinn-dev gmx pdb2gmx -f 1L2Y.pdb -o protein.gro -p topol.top -ff amber99sb-ildn -water tip3p -ignh
docker run -v $(pwd):/data -w /data grinn-dev gmx editconf -f protein.gro -o boxed.gro -c -d 1.0 -bt cubic

# 3. Add water and ions
docker run -v $(pwd):/data -w /data grinn-dev gmx solvate -cp boxed.gro -cs spc216.gro -o solvated.gro -p topol.top
docker run -v $(pwd):/data -w /data grinn-dev gmx grompp -f /app/mdp_files/ions.mdp -c solvated.gro -p topol.top -o ions.tpr
echo "SOL" | docker run -i -v $(pwd):/data -w /data grinn-dev gmx genion -s ions.tpr -o system.gro -p topol.top -pname NA -nname CL -neutral

# 4. Energy minimize (1-2 min)
docker run -v $(pwd):/data -w /data grinn-dev gmx grompp -f /app/mdp_files/minim.mdp -c system.gro -p topol.top -o em.tpr
# You may adjust -nt (number of threads) to a level more suitable for your device.
docker run -v $(pwd):/data -w /data grinn-dev gmx mdrun -deffnm em -v -nt 4 

# 5. Quick production run (30 seconds - just 1000 steps)
docker run -v $(pwd):/data -w /data grinn-dev gmx grompp -f /app/mdp_files/prod.mdp -c em.gro -p topol.top -o prod.tpr
docker run -v $(pwd):/data -w /data grinn-dev gmx mdrun -deffnm prod -v -nt 4

# 6. Extract trajectory for gRINN
echo "Protein" | docker run -i -v $(pwd):/data -w /data grinn-dev gmx trjconv -s prod.tpr -f prod.xtc -o protein_traj.xtc -pbc mol
```

### MDP Files Used

The Docker image includes pre-configured MDP files in `/app/mdp_files/` that are used automatically in the commands above. Here's what they contain:

**ions.mdp** (for adding ions):
```
; ions.mdp - used as input into grompp to generate ions.tpr
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = 50000         ; Maximum number of (minimization) steps to perform

nstlist         = 1         ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet    ; Buffered neighbor searching 
ns_type         = grid      ; Method to determine neighbor list (simple, grid)
coulombtype     = cutoff    ; Treatment of long range electrostatic interactions
rcoulomb        = 1.0       ; Short-range electrostatic cut-off
rvdw            = 1.0       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
```

**minim.mdp** (energy minimization):
```
; minim.mdp - used as input into grompp to generate em.tpr
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = 50000         ; Maximum number of (minimization) steps to perform

nstlist         = 1         ; Frequency to update the neighbor list
cutoff-scheme   = Verlet    ; Buffered neighbor searching
coulombtype     = PME       ; Treatment of long range electrostatic interactions
rcoulomb        = 1.0       ; Short-range electrostatic cut-off
rvdw            = 1.0       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
```

**prod.mdp** (quick production run):
```
; prod.mdp - used for quick production run
title                   = Quick Production MD
integrator              = md        ; leap-frog integrator
nsteps                  = 1000      ; Just 1000 steps (2 ps)
dt                      = 0.002     ; 2 fs
nstxout-compressed      = 10        ; save coordinates every 20 fs
nstenergy               = 10        ; save energies every 20 fs
continuation            = no        ; Starting from minimization
cutoff-scheme           = Verlet    ; Buffered neighbor searching
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = Protein Non-Protein   ; two coupling groups - more accurate
tau_t                   = 0.1     0.1           ; time constant, in ps
ref_t                   = 300     300           ; reference temperature, one for each group, in K
pcoupl                  = Parrinello-Rahman     ; Pressure coupling on
pcoupltype              = isotropic             ; uniform scaling of box vectors
tau_p                   = 5.5                   ; time constant, in ps
ref_p                   = 1.0                   ; reference pressure, in bar
compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1
pbc                     = xyz       ; 3-D PBC
gen_vel                 = yes       ; Generate velocities from Maxwell distribution
gen_temp                = 300       ; Temperature for Maxwell distribution
```

### Expected Runtime
- **Small protein** (< 50 residues): 2-4 minutes total
- **Medium protein** (100 residues): 3-6 minutes total  
- **Large protein** (> 200 residues): 5-10 minutes total

### Simulation Protocol
This streamlined workflow is designed for quick gRINN testing:
1. **Energy minimization**: Remove bad contacts
2. **Quick production**: Generate minimal trajectory for testing (1000 steps = 2 ps)

**Note**: This produces a very short trajectory suitable only for testing gRINN functionality. For research purposes, longer simulations with proper equilibration are recommended.

### Customizing MDP Files (Optional)
The workflow above uses optimized MDP files from `/app/mdp_files/` automatically. If you want to customize them:

```bash
# Copy MDP files from Docker image to current directory for editing
docker run -v $(pwd):/data grinn-dev bash -c "cp /app/mdp_files/*.mdp /data/"

# Edit the files as needed, then use local files instead:
docker run -v $(pwd):/data -w /data grinn-dev gmx grompp -f ./minim.mdp -c system.gro -p topol.top -o em.tpr
```

**Note**: The pre-configured files in the Docker image are optimized for most protein systems and should work well for testing gRINN.

### Troubleshooting
- **"Fatal error"**: Check Docker is running and image is built correctly
- **"Blowing up"**: Try using pre-configured MDP files with more conservative settings
- **Too slow**: Use fewer steps or smaller system
- **Permission issues**: Ensure your current directory is mounted correctly with `-v $(pwd):/data`

### Output Files for gRINN
You need these files (all in .gro format for consistency):
- `em.gro` or `prod.gro` (structure from minimization or production)
- `topol.top` (topology) 
- `protein_traj.xtc` (trajectory from production run - very short!)

### Quality Check
Before running gRINN, verify your simulation generated output:
```bash
# Check trajectory length (should show ~100 frames from 1000 steps)
docker run -v $(pwd):/data -w /data grinn-dev gmx check -f protein_traj.xtc

# Quick energy check
docker run -v $(pwd):/data -w /data grinn-dev gmx energy -f prod.edr -o energy.xvg << EOF
Potential
0
EOF
```

**Next**: Test these files with gRINN → [testing-guide.md](testing-guide.md) 🚀
