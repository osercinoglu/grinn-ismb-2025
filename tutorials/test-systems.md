# Test Systems for gRINN

## Quick Pick List

### Small & Fast (Good for First Test)
- **Trp-cage**: `1L2Y` (20 residues) - Ultra-fast folding
- **Villin**: `2F4K` (35 residues) - Classic test system
- **WW domain**: `1PIN` (37 residues) - Beta-sheet protein

### Medium Size (Realistic Testing)
- **Lysozyme**: `1LYS` (129 residues) - Well-studied enzyme
- **Ubiquitin**: `1UBQ` (76 residues) - Stable, compact
- **Protein G**: `1PGB` (56 residues) - Mixed structure

### With Ligands (Test Binding Sites)
- **Lysozyme + inhibitor**: `1LYZ` - Clear binding site
- **Trypsin + inhibitor**: `1PPH` - Strong binding
- **HIV Protease**: `1HTM` - Drug target

### Large (Performance Testing)
- **Antibody Fab**: `1IGY` (~450 residues) - Multi-chain
- **T4 Lysozyme**: `1L63` (164 residues) - Larger enzyme

## How to Use
1. Pick a protein from the list above
2. Download: `wget https://files.rcsb.org/download/PDBID.pdb`
3. Follow the [quick simulation guide](quick-simulation.md)

## What to Test
- **Does it run?** Basic functionality
- **Are results reasonable?** Strong interactions make sense
- **How fast is it?** Note runtime and memory usage
- **Any errors?** Report everything that breaks

**Tip**: Start with Trp-cage (1L2Y) - it's tiny and runs in minutes! 🚀
