# Contributing to gRINN - CollaborationFest Guide

Welcome to the gRINN CollaborationFest! This guide will help you get started with contributing to the project, whether you're a tester, developer, or documentation contributor.

## Table of Contents
- [Quick Start](#quick-start)
- [Contribution Types](#contribution-types)
- [Setting Up Your Environment](#setting-up-your-environment)
- [Project Structure](#project-structure)
- [Workflow Guidelines](#workflow-guidelines)
- [Getting Help](#getting-help)

## Quick Start

1. **Fork and Clone**: Fork the [main gRINN repository](https://github.com/osercinoglu/grinn) and clone it locally
2. **Set up Environment**: Follow the installation guide in the main repository
3. **Choose Your Focus**: Pick from testing, development, or documentation tasks
4. **Start Contributing**: Follow the specific guides linked below

## Contribution Types

### 🧪 For Testers
**Goal**: Validate gRINN with diverse datasets and identify issues

**What you'll do**:
- Test gRINN with different protein systems
- Run the tool with provided sample data
- Report bugs and unexpected behaviors
- Share results and experiences

**Skills needed**: 
- Basic knowledge of molecular dynamics
- Familiarity with protein structures
- Experience with GROMACS (helpful but not required)

**Start here**: [Testing Guide](./tutorials/testing-guide.md)

### 💻 For Developers
**Goal**: Fix bugs and enhance the codebase

**What you'll do**:
- Fix identified bugs and issues
- Improve code quality and performance
- Enhance the web dashboard
- Optimize algorithms for large systems

**Skills needed**:
- Python programming
- Experience with scientific computing
- Knowledge of web development (for dashboard work)
- Familiarity with Docker (helpful)

**Start here**: [Development Guide](./tutorials/development-guide.md)

### 📚 For Documentation Contributors
**Goal**: Create tutorials and improve documentation

**What you'll do**:
- Write step-by-step tutorials
- Create example workflows
- Document best practices
- Improve code documentation

**Skills needed**:
- Technical writing skills
- Understanding of protein analysis workflows
- Markdown/documentation tools experience

**Start here**: [Documentation Guide](./tutorials/documentation-guide.md)

## Setting Up Your Environment

### Prerequisites
- Python 3.10+
- Git
- GROMACS (for testing with your own data)
- Docker (optional but recommended)

### Quick Setup
```bash
# Clone the repository
git clone https://github.com/osercinoglu/grinn.git
cd grinn

# Follow installation instructions in the main repository
# See: https://github.com/osercinoglu/grinn#installation--usage
```

## Project Structure

```
grinn/
├── grinn_workflow.py          # Main workflow script
├── gRINN_Dashboard/           # Web dashboard code
├── test_data/                 # Sample datasets
├── docs/                      # Documentation
├── tests/                     # Test suite
├── environment.yml            # Conda environment
├── Dockerfile                 # Docker configuration
└── README.md                  # Main documentation
```

## Workflow Guidelines

### Issue Tracking
- Check existing [GitHub issues](https://github.com/osercinoglu/grinn/issues) before starting
- Create new issues for bugs or feature requests
- Use clear, descriptive titles and include relevant details

### Branch Naming
- `bugfix/issue-description`
- `feature/new-feature-name`
- `docs/tutorial-topic`
- `test/system-type`

### Pull Requests
- Create focused PRs that address specific issues
- Include clear descriptions of changes
- Add tests for new features when applicable
- Update documentation as needed

### Communication
- Use GitHub issues for bug reports and feature requests
- Use GitHub discussions for questions and general discussion
- Tag maintainers (@osercinoglu) for urgent issues

## Getting Help

### During CollaborationFest
- Ask questions in the dedicated Slack channel
- Find team coordinators for immediate help
- Join daily check-ins and progress reviews

### Resources
- [Main gRINN Repository](https://github.com/osercinoglu/grinn)
- [Scientific Publication](https://doi.org/10.1093/nar/gky381)
- [Tutorial Files](./tutorials/)
- [GitHub Discussions](https://github.com/osercinoglu/grinn/discussions)

### Common Issues
- **Installation problems**: Check the main repository's installation guide
- **GROMACS errors**: Ensure GROMACS is properly installed and accessible
- **Memory issues**: Use smaller test systems or enable frame skipping
- **Dashboard problems**: Check browser compatibility and port availability

## Next Steps

Choose your contribution path:
- **Testing**: Go to [Testing Guide](./tutorials/testing-guide.md)
- **Development**: Go to [Development Guide](./tutorials/development-guide.md)  
- **Documentation**: Go to [Documentation Guide](./tutorials/documentation-guide.md)

Ready to contribute? Let's make gRINN better together! 🚀
