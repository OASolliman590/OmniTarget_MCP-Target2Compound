# Installation Guide - MCP Drug Discovery Pipeline

This guide covers different installation methods for the MCP Drug Discovery Pipeline.

## 📊 Installation Status: **FULLY FUNCTIONAL WITH MCP SERVICES** ✅

### ✅ **Successfully Installed**
- **Conda Environment**: `mcp-drug-discovery` with Python 3.11
- **AutoDock Vina**: Built from source (v1.2.7-20-g93cdc3d-mod)
- **DeepDTA**: Repository cloned and adapter configured
- **Core Dependencies**: RDKit, FastAPI, Pydantic, NumPy, Pandas, etc.
- **Real Compounds**: 11 sertraline compounds ready for testing
- **MCP Services**: All 6 services running and responding correctly
- **End-to-End Pipeline**: Fully functional with MCP integration

### 🔧 **Technical Notes**
- Vina binary: `vina_custom` (available in PATH)
- DeepDTA: Using placeholder predictions (real model integration pending)
- NumPy 2.x compatibility warnings with RDKit (non-blocking)
- MCP Services: All operational with live data integration

## 🐍 Method 1: Conda Environment (Recommended)

### Prerequisites
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/products/distribution) installed

### Quick Setup
```bash
# Clone the repository
git clone <repository-url>
cd mcp-drug-discovery

# Run the automated setup script
chmod +x conda_env_setup.sh
./conda_env_setup.sh

# Activate the environment
source activate_env.sh
```

### Manual Conda Setup
```bash
# Create environment from environment.yml
conda env create -f environment.yml

# Activate environment
conda activate mcp-drug-discovery
```

### ML Dependencies (Already Installed)

**AutoDock Vina**
```bash
# Vina is already built from source and available as 'vina_custom'
# Binary location: $CONDA_PREFIX/bin/vina_custom
# Version: v1.2.7-20-g93cdc3d-mod
# Built with Boost libraries from conda-forge
```

**DeepDTA**
```bash
# DeepDTA repository is already cloned in third_party/DeepDTA/
# Adapter is configured and working with placeholder predictions
# Ready for real model integration when needed
```

# Verify installation
python -c "import orchestrator; print('✅ Installation successful')"
```

## 🐍 Method 2: Pip Installation

### Prerequisites
- Python 3.11+
- pip

### Setup
```bash
# Clone the repository
git clone <repository-url>
cd mcp-drug-discovery

# Create virtual environment (recommended)
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Install in development mode
pip install -e .
```

## 🐳 Method 3: Docker (Production)

### Prerequisites
- Docker and Docker Compose

### Setup
```bash
# Clone the repository
git clone <repository-url>
cd mcp-drug-discovery

# Build all services
docker compose build

# Start services
docker compose up -d

# Test installation
docker compose exec orchestrator python -c "import orchestrator; print('✅ Installation successful')"
```

## 🔧 Development Setup

### For Contributors
```bash
# Clone with submodules
git clone --recursive <repository-url>
cd mcp-drug-discovery

# Set up conda environment
conda env create -f environment.yml
conda activate mcp-drug-discovery

# Install pre-commit hooks
pre-commit install

# Install development dependencies
pip install -e ".[dev]"

# Run tests
pytest orchestrator/tests/ -v
```

## 📦 Environment Details

### Conda Environment (`environment.yml`)
- **Python**: 3.11
- **Channels**: conda-forge, bioconda, defaults
- **Key packages**:
  - Scientific: numpy, pandas, scipy, scikit-learn
  - Web: fastapi, uvicorn, httpx
  - Validation: pydantic, pydantic-settings
  - CLI: typer
  - Bioinformatics: biopython, rdkit
  - Visualization: plotly, matplotlib
  - Development: pytest, black, mypy

### Pip Requirements (`requirements.txt`)
- Core framework dependencies
- Scientific computing stack
- Development and testing tools
- Optional ML dependencies (commented out)

## 🧪 Verify Installation

### Quick Test
```bash
# Activate environment (if using conda)
conda activate mcp-drug-discovery

# Run quick functionality test
./scripts/quick_test.sh

# Test core imports
python -c "
from orchestrator.pipeline import DrugDiscoveryPipeline
from orchestrator.mcp_clients import KEGGClient
from orchestrator.adapters import DeepDTAAdapter
print('✅ All core components imported successfully')
"
```

### Comprehensive Test
```bash
# Run full test suite
python scripts/test_functionality.py

# Test specific components
python scripts/test_component.py mcp kegg
python scripts/test_component.py adapter deepdta
```

## 🐛 Troubleshooting

### Common Issues

**1. Conda Environment Creation Fails**
```bash
# Update conda
conda update conda

# Clear conda cache
conda clean --all

# Try with mamba (faster)
conda install mamba
mamba env create -f environment.yml
```

**2. RDKit Installation Issues**
```bash
# Install RDKit from conda-forge
conda install -c conda-forge rdkit

# Or use rdkit-pypi
pip install rdkit-pypi
```

**3. Import Errors**
```bash
# Check Python path
python -c "import sys; print(sys.path)"

# Reinstall in development mode
pip install -e .

# Check for conflicting packages
pip list | grep -E "(pydantic|fastapi|numpy)"
```

**4. Docker Build Issues**
```bash
# Clean Docker cache
docker system prune -a

# Rebuild without cache
docker compose build --no-cache

# Check Docker resources
docker system df
```

### Platform-Specific Notes

**macOS (Apple Silicon)**
```bash
# Some packages may need specific channels
conda install -c conda-forge rdkit
conda install -c conda-forge openbabel
```

**Windows**
```bash
# Use conda instead of pip for system dependencies
conda install -c conda-forge rdkit
conda install -c conda-forge openbabel

# For WSL2, follow Linux instructions
```

**Linux**
```bash
# Install system dependencies
sudo apt-get update
sudo apt-get install build-essential

# For Ubuntu/Debian
sudo apt-get install libopenbabel-dev
```

## 🔄 Environment Management

### Updating Environment
```bash
# Update conda environment
conda env update -f environment.yml

# Update pip requirements
pip install -r requirements.txt --upgrade

# Update specific packages
conda update numpy pandas
pip install --upgrade fastapi pydantic
```

### Environment Information
```bash
# List environments
conda env list

# Show environment details
conda list -n mcp-drug-discovery

# Export environment
conda env export -n mcp-drug-discovery > environment_backup.yml
```

### Cleanup
```bash
# Remove environment
conda env remove -n mcp-drug-discovery

# Clean conda cache
conda clean --all

# Remove Docker images
docker compose down --rmi all
```

## 🚀 Next Steps

After successful installation:

1. **Start MCP Services**:
   ```bash
   docker compose -f docker/docker-compose.yaml up -d
   ```

2. **Test MCP Connections**:
   ```bash
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

3. **Test Pipeline**:
   ```bash
   python scripts/test_functionality.py
   python -m orchestrator.cli run configs/test_sertraline.yaml
   ```

4. **Use API**:
   ```bash
   curl http://localhost:8000/health
   ```

5. **Read Documentation**:
   - [README.md](README.md) - Overview and usage
   - [TESTING.md](TESTING.md) - Testing guide
   - [API Documentation](http://localhost:8000/docs) - When running

## 📞 Support

If you encounter issues:

1. Check the [troubleshooting section](#-troubleshooting)
2. Review [TESTING.md](TESTING.md) for testing procedures
3. Check [GitHub Issues](https://github.com/your-org/mcp-drug-discovery/issues)
4. Join our [Discord community](https://discord.gg/mcp-drug-discovery)

---

**Happy drug discovery! 🧬💊**
