# Contributing to MCP Drug Discovery Pipeline

Thank you for your interest in contributing to the MCP Drug Discovery Pipeline! This document provides guidelines and information for contributors.

## 🚀 Quick Start for Contributors

### Prerequisites
- Python 3.11+
- Docker & Docker Compose
- Git
- Basic understanding of drug discovery and bioinformatics

### Setup Development Environment

1. **Fork and Clone**
   ```bash
   git clone https://github.com/YOUR_USERNAME/mcp-drug-discovery.git
   cd mcp-drug-discovery
   ```

2. **Set up Environment**
   ```bash
   # Automated setup
   chmod +x conda_env_setup.sh
   ./conda_env_setup.sh
   source activate_env.sh
   
   # Or manual setup
   conda env create -f environment.yml
   conda activate mcp-drug-discovery
   ```

3. **Start MCP Services**
   ```bash
   docker compose -f docker/docker-compose.yaml up -d
   ```

4. **Run Tests**
   ```bash
   python scripts/test_functionality.py
   ```

## 🏗️ Project Structure

```
mcp-drug-discovery/
├── orchestrator/           # Main pipeline code
│   ├── adapters/          # ML model adapters
│   ├── mcp_clients/       # MCP server clients
│   ├── schemas/           # Data models
│   ├── scoring/           # Scoring algorithms
│   └── pipeline.py        # Main orchestration
├── services/              # MCP server implementations
├── configs/               # Configuration files
├── scripts/               # Utility scripts
├── data/                  # Data files
├── docker/                # Docker configurations
└── docs/                  # Documentation
```

## 🔧 Development Guidelines

### Code Style
- Follow PEP 8 for Python code
- Use type hints for all functions
- Document all public methods and classes
- Use meaningful variable and function names

### Testing
- Write tests for new functionality
- Ensure all tests pass before submitting PR
- Test with real data when possible
- Use the existing test framework

### Documentation
- Update relevant documentation for new features
- Add docstrings to new functions and classes
- Update README.md if adding new features
- Include examples in documentation

## 🐛 Bug Reports

When reporting bugs, please include:

1. **Environment Information**
   - OS and version
   - Python version
   - Docker version
   - Conda environment details

2. **Reproduction Steps**
   - Clear steps to reproduce the issue
   - Configuration files used
   - Input data (if applicable)

3. **Error Information**
   - Complete error messages
   - Stack traces
   - Log files

4. **Expected vs Actual Behavior**
   - What you expected to happen
   - What actually happened

## ✨ Feature Requests

When requesting features, please include:

1. **Use Case**
   - What drug discovery problem does this solve?
   - How would this improve the pipeline?

2. **Implementation Details**
   - Which components need modification?
   - Any dependencies or external services?
   - Performance considerations?

3. **Examples**
   - Sample input/output
   - Configuration examples
   - Usage scenarios

## 🔄 Pull Request Process

1. **Create Feature Branch**
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. **Make Changes**
   - Write code following project guidelines
   - Add tests for new functionality
   - Update documentation

3. **Test Changes**
   ```bash
   # Run all tests
   python scripts/test_functionality.py
   
   # Test specific components
   python scripts/test_component.py adapters all
   python scripts/test_component.py mcp all
   ```

4. **Commit Changes**
   ```bash
   git add .
   git commit -m "Add: brief description of changes"
   ```

5. **Push and Create PR**
   ```bash
   git push origin feature/your-feature-name
   ```

## 🧪 Testing Guidelines

### Running Tests
```bash
# Full functionality test
python scripts/test_functionality.py

# Individual component tests
python scripts/test_component.py adapters all
python scripts/test_component.py mcp all
python scripts/test_component.py scoring all

# MCP connection tests
docker compose -f docker/docker-compose.yaml up -d
python -c "
import asyncio
from orchestrator.mcp_clients import *
async def test(): 
    clients = [ReactomeClient(), ProteinAtlasClient(), STRINGClient(), UniProtClient(), PDBClient(), ChEMBLClient()]
    for client in clients: 
        print(f'{client.__class__.__name__}: {\"✅\" if await client.test_connection() else \"❌\"}')
        await client.close()
asyncio.run(test())
"
```

### Writing Tests
- Test both success and failure cases
- Use real data when possible
- Mock external services appropriately
- Ensure tests are deterministic

## 📚 Component Development

### Adding New MCP Clients
1. Create client in `orchestrator/mcp_clients/`
2. Inherit from `BaseMCPClient`
3. Implement required methods
4. Add to `__init__.py`
5. Create corresponding service in `services/`

### Adding New ML Adapters
1. Create adapter in `orchestrator/adapters/`
2. Implement `setup()` and prediction methods
3. Add to `__init__.py`
4. Update pipeline integration

### Adding New Scoring Methods
1. Create scorer in `orchestrator/scoring/`
2. Implement normalization and fusion methods
3. Add to pipeline scoring workflow
4. Update configuration schemas

## 🐳 Docker Development

### Building Services
```bash
# Build all services
docker compose -f docker/docker-compose.yaml build

# Build specific service
docker compose -f docker/docker-compose.yaml build kegg-mcp
```

### Testing Services
```bash
# Start all services
docker compose -f docker/docker-compose.yaml up -d

# Check service logs
docker compose -f docker/docker-compose.yaml logs service-name

# Stop services
docker compose -f docker/docker-compose.yaml down
```

## 📖 Documentation

### Updating Documentation
- Update relevant `.md` files for new features
- Keep examples current and working
- Update API documentation if applicable
- Maintain consistency across all docs

### Documentation Structure
- `README.md`: Overview and quick start
- `INSTALLATION.md`: Detailed installation guide
- `TESTING.md`: Testing strategies and procedures
- `STATUS.md`: Current project status
- `CONTRIBUTING.md`: This file

## 🤝 Community Guidelines

- Be respectful and inclusive
- Help others learn and grow
- Provide constructive feedback
- Follow the code of conduct
- Ask questions when unsure

## 📞 Getting Help

- **Issues**: Use GitHub issues for bugs and feature requests
- **Discussions**: Use GitHub discussions for questions and ideas
- **Documentation**: Check existing docs first
- **Examples**: Look at existing configurations and scripts

## 🎯 Areas for Contribution

### High Priority
- Real DeepDTA model integration
- Additional MCP server implementations
- Performance optimizations
- Extended test coverage

### Medium Priority
- Additional ML model adapters
- New scoring algorithms
- API improvements
- Documentation enhancements

### Low Priority
- UI/Web interface
- Additional data formats
- Advanced visualization
- Integration with other tools

## 📋 Checklist for Contributors

Before submitting a PR, ensure:

- [ ] Code follows project style guidelines
- [ ] All tests pass
- [ ] Documentation is updated
- [ ] No hardcoded data or credentials
- [ ] Error handling is appropriate
- [ ] Performance is acceptable
- [ ] Security considerations addressed

Thank you for contributing to the MCP Drug Discovery Pipeline! 🧬
