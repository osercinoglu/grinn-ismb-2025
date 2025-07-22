# Test Systems for gRINN

## Quick Pick List

### Medium Size (Realistic Testing)
- **Lysozyme**: `1LYS` (129 residues) - Well-studied enzyme
- **Ubiquitin**: `1UBQ` (76 residues) - Stable, compact
- **Protein G**: `1PGB` (56 residues) - Mixed structure

### Large (Performance Testing)
- **Antibody Fab**: `1IGY` (~450 residues) - Multi-chain
- **T4 Lysozyme**: `1L63` (164 residues) - Larger enzyme

## How to Use
1. Pick a protein from the list above
2. Download: `wget https://files.rcsb.org/download/PDBID.pdb`
3. Follow the [quickstart guide](quickstart.md)

## What to Test
- **Does it run?** Basic functionality
- **Are results reasonable?** Strong interactions make sense
- **How fast is it?** Note runtime and memory usage
- **Any errors?** Report everything that breaks

**Tip**: Start with Trp-cage (1L2Y) - it's tiny and runs in minutes! 🚀
