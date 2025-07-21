# Documentation Guide for gRINN

## Overview
This guide helps contributors create clear, useful documentation for the gRINN project. Good documentation makes tools accessible to the wider community.

## Prerequisites for Documentation Writing
- Successfully tested gRINN with at least one protein system
- Understanding of the Docker-based workflow
- Experience with the specific use case you want to document

## Documentation Principles

### Write for Your Audience
- **Beginners**: Assume minimal MD simulation experience
- **Experienced users**: Focus on gRINN-specific details
- **Developers**: Include technical implementation details

### Be Practical
- Provide complete, working examples
- Include expected outputs and runtimes
- Mention common pitfalls and solutions
- Test all commands before publishing

### Keep It Current
- Use Docker-based commands consistently
- Reference current file formats (.gro consistently)
- Update examples with recent protein structures

## Types of Documentation Needed

### 1. Tutorial Documentation
**Format**: Step-by-step guides
**Audience**: New users
**Examples**: Current guides in `tutorials/` folder

### 2. How-To Guides
**Format**: Problem-focused instructions
**Audience**: Users with specific goals
**Examples**: "How to analyze membrane proteins", "How to handle large trajectories"

### 3. Reference Documentation
**Format**: Comprehensive parameter lists
**Audience**: All users
**Examples**: Complete command-line reference, file format specifications

### 4. Troubleshooting Guides
**Format**: Problem → Solution mapping
**Audience**: Users encountering issues
**Examples**: Error message database, performance optimization guide

## Documentation Standards

### File Naming Convention
- Use kebab-case: `protein-ligand-analysis.md`
- Be descriptive: `membrane-protein-simulation.md` not `membrane.md`
- Include version if needed: `advanced-analysis-v2.md`

### Structure Template
```markdown
# Clear, Descriptive Title

## Overview
Brief description of what this guide covers and who it's for.

## Prerequisites
- Required software/setup
- Background knowledge needed
- Previous tutorials to complete

## Step-by-Step Instructions
### Step 1: Setup
Clear, numbered steps with code blocks

### Step 2: Analysis
More steps...

## Expected Results
What users should see if successful

## Troubleshooting
Common issues and solutions

## Next Steps
Where to go from here
```

### Code Block Standards
Always use Docker-based commands:
```bash
# ✅ Good - Docker command with volume mount
docker run -v $(pwd):/data grinn-dev workflow /data/protein.gro /data/results \
  --top /data/topol.top --traj /data/traj.xtc

# ❌ Avoid - Local installation commands
python grinn_workflow.py protein.gro results/ --top topol.top --traj traj.xtc
```

### File Format Consistency
- Use `.gro` format for structure files consistently
- Always mention topology (.top) and trajectory (.xtc) requirements
- Specify full file paths in examples

## Documentation Workflow

### 1. Planning Phase
- **Identify the gap**: What documentation is missing?
- **Define scope**: What specific problem will this solve?
- **Check audience**: Who needs this information?
- **Research**: Test the workflow thoroughly yourself

### 2. Writing Phase
- **Start with outline**: Plan the structure first
- **Write iteratively**: Get basics down, then refine
- **Include examples**: Real commands, real outputs
- **Test everything**: Every command should work

### 3. Review Phase
- **Self-review**: Read from a new user's perspective
- **Test with others**: Have someone else follow your guide
- **Update based on feedback**: Fix unclear parts
- **Finalize**: Polish language and formatting

## Specific Documentation Tasks

### High-Priority Needs
1. **Membrane protein analysis guide**
   - Special considerations for membrane systems
   - Appropriate simulation parameters
   - Interpretation of results

2. **Large system optimization guide**
   - Memory management strategies
   - Performance tuning
   - Parallelization options

3. **Error troubleshooting database**
   - Common error messages and solutions
   - System-specific issues
   - Recovery strategies

4. **Advanced analysis examples**
   - Custom energy cutoffs
   - Specific residue selection
   - Time-dependent analysis

### Medium-Priority Needs
1. **Comparison with other tools**
   - How gRINN differs from alternatives
   - When to use gRINN vs other approaches
   - Migration guides from other tools

2. **Integration guides**
   - Using gRINN with existing workflows
   - Scripting and automation
   - Data export/import

3. **Visualization customization**
   - Dashboard customization
   - Export options
   - Publication-ready figures

## Writing Tips

### Language and Style
- **Be concise**: Get to the point quickly
- **Use active voice**: "Run the command" not "The command should be run"
- **Include context**: Explain why steps are necessary
- **Define terms**: Explain technical concepts briefly

### Examples and Code
- **Complete examples**: Show full commands, not fragments
- **Real data**: Use actual proteins and realistic parameters
- **Expected timing**: Mention how long things should take
- **Error handling**: Show what to do when things go wrong

### Visual Elements
- **Screenshots**: Include dashboard screenshots when relevant
- **Diagrams**: Use simple diagrams for complex workflows
- **Tables**: Organize parameter lists and comparisons
- **Callout boxes**: Highlight important warnings or tips

## Testing Your Documentation

### Self-Testing
1. **Fresh environment**: Test in clean directory
2. **Follow exactly**: Don't skip steps or add knowledge
3. **Time it**: Note how long each step takes
4. **Document issues**: Note any problems encountered

### User Testing
1. **Find a volunteer**: Someone unfamiliar with your specific workflow
2. **Watch them work**: Don't help unless they're truly stuck
3. **Note confusion points**: Where do they hesitate or make mistakes?
4. **Gather feedback**: What was unclear or missing?

### Quality Checklist
- [ ] All commands tested and working
- [ ] File paths and names consistent
- [ ] Prerequisites clearly listed
- [ ] Expected outputs described
- [ ] Common errors addressed
- [ ] Timing estimates provided
- [ ] Next steps suggested

## Documentation Tools and Resources

### Markdown Editing
- Use any text editor with Markdown preview
- VS Code with Markdown extensions
- Typora for WYSIWYG editing

### Screen Capture
- Screenshot tools for dashboard images
- Screen recording for complex procedures
- Image optimization for web display

### Collaboration
- GitHub for version control and review
- Issues for planning and discussion
- Pull requests for review process

## Example Documentation Projects

### Beginner: "Analyzing Your First Protein"
- Target audience: Complete beginners
- Scope: Single small protein from PDB to dashboard
- Length: ~1000 words
- Effort: 2-4 hours

### Intermediate: "Membrane Protein Analysis Workflow"
- Target audience: Users with MD experience
- Scope: Special considerations for membrane systems
- Length: ~2000 words
- Effort: 1-2 days

### Advanced: "Performance Optimization for Large Systems"
- Target audience: Power users and developers
- Scope: Technical optimization strategies
- Length: ~3000 words
- Effort: 3-5 days

## Contributing Your Documentation

### Submission Process
1. **Create branch**: `git checkout -b docs/your-topic`
2. **Write documentation**: Follow standards above
3. **Test thoroughly**: Verify all examples work
4. **Submit PR**: Include description of audience and scope
5. **Respond to feedback**: Address reviewer comments

### Review Criteria
- **Accuracy**: All technical content correct
- **Completeness**: Covers stated scope adequately
- **Clarity**: Accessible to target audience
- **Consistency**: Follows project standards
- **Utility**: Solves real user problems

Remember: Good documentation is just as valuable as good code! 📚
