# Documentation Guide

## Prerequisites for Documentation Contributors

**You should write tutorials ONLY if:**
- ✅ You've successfully tested gRINN with your own protein system
- ✅ You found the tool useful for a specific use case
- ✅ You have an idea for a tutorial that would help others
- ✅ You encountered and solved specific problems during testing

**Why this matters**: The best tutorials come from real experience using the tool!

## What We Need

### Tutorials
- **Step-by-step guides** for specific use cases
- **Troubleshooting sections** for common problems
- **Example workflows** with real data

### API Documentation
- **Function descriptions** with clear examples
- **Parameter explanations** with default values
- **Return value descriptions** with expected formats

### User Guides
- **Installation troubleshooting** for different systems
- **Performance tips** for large datasets
- **Best practices** for MD simulation setup

## Writing Guidelines

### Keep It Simple
- **Short sentences** and paragraphs
- **Clear headings** and structure
- **Code examples** that actually work
- **Expected outputs** for each step

### Test Everything
- **Run all commands** you write about
- **Check all links** actually work
- **Verify examples** produce expected results
- **Test on different systems** if possible

### Template for Tutorials
```markdown
# Title: What This Does

## Goal
One sentence describing what users will accomplish.

## Prerequisites
- List what they need installed
- Link to installation guides

## Step-by-Step
1. **Do this**: `command here`
2. **Check that**: Expected output
3. **If problems**: Common fixes

## Troubleshooting
- **Error X**: Solution Y
- **Problem Z**: Fix A

## Next Steps
Link to related guides.
```

## Quick Tasks

### Improve Existing Docs
1. **Test the examples**: Do they work?
2. **Fill gaps**: What's missing?
3. **Fix errors**: Typos, broken links, wrong commands
4. **Add screenshots**: Especially for dashboard features

### Write New Content
- **Use case tutorials**: "How to analyze protein-protein interfaces"
- **Troubleshooting guides**: "gRINN won't start"
- **Comparison guides**: "gRINN vs other tools"
- **Performance guides**: "Optimizing for large systems"

### Tools You Can Use
- **Markdown**: For most documentation
- **Screenshots**: For dashboard features
- **Diagrams**: For workflow illustrations
- **GIFs**: For interactive features

## Where to Contribute
- **Main repo**: Update README, add examples
- **This repo**: Improve tutorials, add guides
- **Wiki pages**: Create in-depth explanations
- **Video tutorials**: Record walkthroughs

## Getting Started

### Before You Write
1. **Test gRINN thoroughly** with your protein system
2. **Note what worked well** and what was confusing
3. **Identify your use case**: What specific problem did gRINN solve for you?
4. **Think about your audience**: Who else would benefit from your experience?

### Good Tutorial Ideas (Based on Real Experience)
- **"Analyzing Protein-Ligand Binding Sites with gRINN"** (after testing with ligand complexes)
- **"gRINN for Large Protein Complexes: Performance Tips"** (after testing >300 residue systems)
- **"Troubleshooting gRINN with Membrane Proteins"** (after solving specific membrane protein issues)
- **"Comparing Force Fields using gRINN Networks"** (after testing different force fields)

## Getting Started
1. **Pick something confusing**: What did you struggle with during testing?
2. **Write it down**: Make it clearer based on your experience
3. **Test with others**: Does it help people who haven't used gRINN?
4. **Submit**: Create pull request

**Remember**: Your real testing experience is what makes tutorials valuable! 📖
