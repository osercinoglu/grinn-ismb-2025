# Documentation Guide - gRINN CollaborationFest

This guide will help you contribute to gRINN's documentation, create tutorials, and improve the overall user experience.

## Overview

Documentation contributors help by:
- Writing step-by-step tutorials for different use cases
- Creating example workflows with sample data
- Improving existing documentation clarity
- Developing troubleshooting guides

## Prerequisites

- **Strong technical writing skills**
- **Understanding of protein analysis workflows**
- **Basic knowledge of molecular dynamics**
- **Markdown formatting experience**
- **Familiarity with gRINN functionality** (can be learned)

## Types of Documentation Contributions

### 1. User Tutorials
Step-by-step guides for specific analysis scenarios

### 2. API Documentation
Code documentation and function references

### 3. Troubleshooting Guides
Common problems and their solutions

### 4. Example Workflows
Complete analysis examples with sample data

### 5. Best Practices
Guidelines for effective use of gRINN

## Getting Started

### 1. Set Up Documentation Environment
```bash
# Clone the repository
git clone https://github.com/osercinoglu/grinn.git
cd grinn

# Install gRINN to understand functionality
conda env create -f environment.yml
conda activate grinn

# Test basic functionality
python grinn_workflow.py --help
```

### 2. Explore Existing Documentation
- Read the main README
- Try existing examples
- Identify gaps or unclear sections
- Note common user questions

### 3. Choose Your Focus
Pick a documentation area that matches your interests and expertise.

## Tutorial Writing Guidelines

### Structure Template
```markdown
# Tutorial Title

## Overview
Brief description of what this tutorial covers and who it's for.

## Prerequisites
- Required knowledge
- Software requirements
- Data requirements

## Learning Objectives
By the end of this tutorial, you will be able to:
- Objective 1
- Objective 2
- Objective 3

## Step-by-Step Instructions

### Step 1: [Clear Action Title]
Explanation of what we're doing and why.

[Code block or command]

Expected output or what to look for.

### Step 2: [Next Action]
Continue with clear steps...

## Troubleshooting
Common issues and solutions.

## Summary
Recap of what was accomplished.

## Next Steps
Where to go from here.
```

## Priority Tutorial Topics

### 1. Basic Protein Analysis
**Target audience**: New users
**Content**:
- Installing gRINN
- Running first analysis
- Understanding output files
- Viewing results in dashboard

```markdown
# Basic Protein Energy Network Analysis

## Overview
This tutorial guides you through your first protein energy network analysis using gRINN with a small example protein.

## Prerequisites
- gRINN installed and working
- Basic understanding of protein structure

## What You'll Learn
- How to run a basic gRINN analysis
- How to interpret energy output files
- How to visualize results in the dashboard

## Step 1: Prepare Your Data
We'll use the provided sample data for a small protein system.

Download the sample files:
- `small_protein.pdb` - Protein structure
- `small_protein.top` - GROMACS topology
- `small_protein.xtc` - MD trajectory

## Step 2: Run Basic Analysis
```bash
python grinn_workflow.py small_protein.pdb results/ \
  --top small_protein.top \
  --traj small_protein.xtc \
  --nt 4
```

**What this command does**:
- Analyzes `small_protein.pdb` structure
- Uses topology file for force field parameters
- Processes trajectory to calculate energies
- Uses 4 CPU cores for faster processing

Expected runtime: 2-5 minutes

## Step 3: Examine Output Files
After completion, check the `results/` directory:
...
```

### 2. Protein-Ligand Analysis
**Target audience**: Drug discovery researchers
**Content**:
- Setting up protein-ligand systems
- Using residue selections
- Analyzing binding site interactions
- Interpreting ligand-specific results

### 3. Large System Optimization
**Target audience**: Advanced users
**Content**:
- Performance optimization strategies
- Memory management
- Frame skipping techniques
- Parallel processing options

### 4. Dashboard Usage
**Target audience**: All users
**Content**:
- Launching the dashboard
- Uploading results
- Navigating visualizations
- Exporting data and figures

### 5. Custom Analysis Workflows
**Target audience**: Researchers with specific needs
**Content**:
- Advanced residue selections
- Custom energy cutoffs
- Batch processing multiple systems
- Integration with other tools

## Writing Best Practices

### Clear and Concise Language
- Use simple, direct sentences
- Define technical terms when first used
- Avoid jargon when possible
- Write for your target audience level

### Code Examples
- Always test code before including it
- Use realistic examples, not toy data
- Include expected outputs
- Explain what each parameter does

```markdown
# Good example:
```bash
python grinn_workflow.py protein.pdb output/ \
  --top protein.top \          # GROMACS topology file
  --traj trajectory.xtc \      # MD trajectory
  --nt 8 \                     # Use 8 CPU cores
  --skip 10                    # Analyze every 10th frame
```

Expected output:
```
INFO: Starting gRINN workflow
INFO: Processing 1000 frames (with skip=10: 100 frames)
INFO: Calculating interaction energies...
```
```

### Visual Aids
- Include screenshots when helpful
- Use diagrams for complex concepts
- Show expected output examples
- Highlight important information

### Error Handling
- Include common error messages
- Provide clear solutions
- Explain what causes errors
- Offer alternative approaches

## API Documentation

### Function Documentation Template
```python
def calculate_interaction_energies(pairs, trajectory, topology, cutoff=10.0):
    """
    Calculate pairwise residue interaction energies from MD trajectory.
    
    This function processes a molecular dynamics trajectory to compute
    non-bonded interaction energies between specified residue pairs.
    
    Parameters
    ----------
    pairs : list of tuples
        Residue index pairs to analyze, e.g., [(0, 1), (1, 2)]
    trajectory : str
        Path to GROMACS trajectory file (.xtc or .trr)
    topology : str
        Path to GROMACS topology file (.top)
    cutoff : float, default 10.0
        Distance cutoff in Angstroms for initial pair filtering
        
    Returns
    -------
    dict
        Dictionary mapping pair tuples to energy arrays
        Keys: (res1_idx, res2_idx)
        Values: numpy arrays of energies per frame
        
    Raises
    ------
    FileNotFoundError
        If trajectory or topology files don't exist
    ValueError
        If residue indices are invalid
        
    Examples
    --------
    >>> pairs = [(0, 1), (5, 10)]
    >>> energies = calculate_interaction_energies(
    ...     pairs, 'traj.xtc', 'topol.top', cutoff=12.0
    ... )
    >>> print(energies[(0, 1)].mean())
    -15.4
    
    Notes
    -----
    Energy calculations include both electrostatic and van der Waals
    components. Positive values indicate repulsive interactions,
    negative values indicate attractive interactions.
    
    The cutoff parameter is used for initial geometric filtering but
    does not affect the final energy calculations.
    """
```

## Troubleshooting Documentation

### Common Issues Template
```markdown
## Issue: Out of Memory Error

### Symptoms
```
MemoryError: Unable to allocate 8.5 GiB for an array
```

### Causes
- Trajectory file too large for available RAM
- System has too many residues
- No frame skipping used

### Solutions

#### Option 1: Use Frame Skipping
```bash
python grinn_workflow.py protein.pdb output/ \
  --traj trajectory.xtc \
  --skip 10  # Analyze every 10th frame
```

#### Option 2: Reduce System Size
- Remove solvent molecules
- Focus on specific protein regions
- Use shorter trajectory segments

#### Option 3: Increase Available Memory
- Close other applications
- Use system with more RAM
- Consider cloud computing resources

### Prevention
- Always start with small test cases
- Monitor memory usage during analysis
- Use appropriate skip values for large systems
```

## Documentation Testing

### Content Validation
- Test all code examples
- Verify file paths and names
- Check that outputs match descriptions
- Ensure links work correctly

### User Testing
- Have someone else follow your tutorial
- Note where they get confused
- Revise based on feedback
- Test with different experience levels

## Contribution Workflow

### 1. Plan Your Documentation
- Identify the need (user question, missing guide)
- Define target audience
- Outline content structure
- Gather necessary examples and data

### 2. Create Draft
- Write in Markdown format
- Include code examples and outputs
- Add screenshots if helpful
- Follow style guidelines

### 3. Test Your Documentation
- Follow your own instructions step-by-step
- Test on a fresh system if possible
- Verify all examples work
- Check for typos and formatting

### 4. Submit for Review
```bash
# Create branch for documentation
git checkout -b docs/tutorial-topic

# Add your files
git add tutorials/new-tutorial.md
git commit -m "Add tutorial for protein-ligand analysis

- Step-by-step guide for drug discovery workflows
- Includes residue selection examples
- Added troubleshooting section"

# Push and create PR
git push origin docs/tutorial-topic
```

## Documentation Tools

### Markdown Editors
- **Typora**: WYSIWYG Markdown editor
- **Mark Text**: Real-time preview editor
- **VS Code**: With Markdown extensions
- **GitHub's editor**: For quick edits

### Screenshot Tools
- **macOS**: Built-in Screenshot (Cmd+Shift+4)
- **Windows**: Snipping Tool or Snip & Sketch
- **Linux**: GNOME Screenshot or Spectacle

### Diagram Creation
- **Draw.io**: Free online diagram tool
- **Lucidchart**: Professional diagramming
- **ASCII diagrams**: For simple text-based visuals

## Style Guidelines

### Formatting
- Use clear headings hierarchy (H1, H2, H3)
- Include table of contents for long documents
- Use code blocks with syntax highlighting
- Use tables for structured information

### Tone and Voice
- Write in second person ("you will...")
- Use active voice
- Be encouraging and supportive
- Acknowledge when tasks are difficult

### File Organization
```
docs/
├── tutorials/
│   ├── basic-analysis.md
│   ├── protein-ligand.md
│   ├── advanced-options.md
│   └── troubleshooting.md
├── api/
│   ├── functions.md
│   └── classes.md
├── examples/
│   ├── workflows/
│   └── datasets/
└── images/
    ├── screenshots/
    └── diagrams/
```

## Measuring Documentation Success

### Metrics to Track
- User questions reduced
- Tutorial completion rates
- GitHub issues related to documentation
- Community feedback

### Continuous Improvement
- Regular review of documentation
- Update examples with new features
- Incorporate user feedback
- Keep information current

## Resources

### Documentation Tools
- [Markdown Guide](https://www.markdownguide.org/)
- [GitHub Markdown](https://guides.github.com/features/mastering-markdown/)
- [Technical Writing Courses](https://developers.google.com/tech-writing)

### Science Communication
- [Writing for Scientists](https://www.nature.com/articles/d41586-018-02404-4)
- [Clear Scientific Writing](https://www.americanscientist.org/article/the-science-of-scientific-writing)

### gRINN-Specific
- [Main Repository](https://github.com/osercinoglu/grinn)
- [Scientific Paper](https://doi.org/10.1093/nar/gky381)
- [Existing Issues](https://github.com/osercinoglu/grinn/issues)

## Getting Help

- **Content questions**: Ask maintainers or experienced users
- **Technical writing**: Use online resources and style guides
- **gRINN functionality**: Test the software hands-on
- **Community feedback**: Engage with users to understand needs

Ready to improve gRINN's documentation? Start with a topic you're passionate about! 📚
