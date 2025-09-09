#!/bin/bash
set -e

echo "🐍 Setting up MCP Drug Discovery Conda Environment"
echo "=================================================="

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Function to print status
print_status() {
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✅ $1 - SUCCESS${NC}"
    else
        echo -e "${RED}❌ $1 - FAILED${NC}"
    fi
}

print_info() {
    echo -e "${BLUE}ℹ️  $1${NC}"
}

print_warning() {
    echo -e "${YELLOW}⚠️  $1${NC}"
}

# Check if conda is installed
if ! command -v conda &> /dev/null; then
    echo -e "${RED}❌ Conda is not installed or not in PATH${NC}"
    echo "Please install Miniconda or Anaconda first:"
    echo "  https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

print_info "Conda found: $(conda --version)"

# Check if environment already exists
ENV_NAME="mcp-drug-discovery"
if conda env list | grep -q "^${ENV_NAME} "; then
    print_warning "Environment '${ENV_NAME}' already exists"
    read -p "Do you want to remove and recreate it? (y/N): " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        print_info "Removing existing environment..."
        conda env remove -n ${ENV_NAME} -y
        print_status "Environment removal"
    else
        print_info "Using existing environment"
    fi
fi

# Create environment from environment.yml
print_info "Creating conda environment from environment.yml..."
conda env create -f environment.yml
print_status "Environment creation"

# Activate environment
print_info "Activating environment..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate ${ENV_NAME}
print_status "Environment activation"

# Verify installation
print_info "Verifying installation..."

# Test Python version
PYTHON_VERSION=$(python --version 2>&1 | cut -d' ' -f2)
print_info "Python version: ${PYTHON_VERSION}"

# Test key imports
python -c "
import sys
print(f'Python executable: {sys.executable}')

# Test core imports
try:
    import numpy as np
    print(f'✅ NumPy {np.__version__}')
except ImportError as e:
    print(f'❌ NumPy import failed: {e}')

try:
    import pandas as pd
    print(f'✅ Pandas {pd.__version__}')
except ImportError as e:
    print(f'❌ Pandas import failed: {e}')

try:
    import pydantic
    print(f'✅ Pydantic {pydantic.__version__}')
except ImportError as e:
    print(f'❌ Pydantic import failed: {e}')

try:
    import fastapi
    print(f'✅ FastAPI {fastapi.__version__}')
except ImportError as e:
    print(f'❌ FastAPI import failed: {e}')

try:
    import rdkit
    print(f'✅ RDKit {rdkit.__version__}')
except ImportError as e:
    print(f'❌ RDKit import failed: {e}')

try:
    import biopython
    print(f'✅ BioPython {biopython.__version__}')
except ImportError as e:
    print(f'❌ BioPython import failed: {e}')
"

print_status "Import verification"

# Test orchestrator imports
print_info "Testing orchestrator imports..."
python -c "
try:
    from orchestrator.settings import settings
    print('✅ Orchestrator settings imported')
    
    from orchestrator.schemas.config import RunConfig
    print('✅ Configuration schemas imported')
    
    from orchestrator.pipeline import DrugDiscoveryPipeline
    print('✅ Pipeline imported')
    
    from orchestrator.mcp_clients import KEGGClient
    print('✅ MCP clients imported')
    
    from orchestrator.adapters import DeepDTAAdapter
    print('✅ ML adapters imported')
    
    print('🎉 All orchestrator components imported successfully!')
    
except ImportError as e:
    print(f'❌ Orchestrator import failed: {e}')
    print('This is expected if you haven\'t set up the project structure yet')
"

print_status "Orchestrator import test"

# Create activation script
print_info "Creating activation script..."
cat > activate_env.sh << 'EOF'
#!/bin/bash
# MCP Drug Discovery Environment Activation Script

echo "🐍 Activating MCP Drug Discovery Environment"
echo "==========================================="

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "❌ Conda not found. Please install Miniconda/Anaconda first."
    exit 1
fi

# Initialize conda for bash
eval "$(conda shell.bash hook)"

# Activate environment
conda activate mcp-drug-discovery

if [ $? -eq 0 ]; then
    echo "✅ Environment activated successfully"
    echo "📁 Current directory: $(pwd)"
    echo "🐍 Python: $(which python)"
    echo "📦 Environment: $(conda info --envs | grep '*' | awk '{print $1}')"
    echo ""
    echo "🚀 Ready to use the MCP Drug Discovery Pipeline!"
    echo ""
    echo "Quick commands:"
    echo "  python -m orchestrator.cli --help"
    echo "  python scripts/quick_test.sh"
    echo "  python -m orchestrator.cli status"
    echo ""
else
    echo "❌ Failed to activate environment"
    exit 1
fi
EOF

chmod +x activate_env.sh
print_status "Activation script creation"

# Create deactivation reminder
print_info "Creating deactivation reminder..."
cat > deactivate_env.sh << 'EOF'
#!/bin/bash
echo "🐍 Deactivating MCP Drug Discovery Environment"
conda deactivate
echo "✅ Environment deactivated"
EOF

chmod +x deactivate_env.sh
print_status "Deactivation script creation"

# Summary
echo ""
echo -e "${GREEN}🎉 Conda Environment Setup Complete!${NC}"
echo "=================================="
echo ""
echo -e "${BLUE}Environment Details:${NC}"
echo "  Name: ${ENV_NAME}"
echo "  Python: ${PYTHON_VERSION}"
echo "  Location: $(conda info --envs | grep ${ENV_NAME} | awk '{print $2}')"
echo ""
echo -e "${BLUE}Usage:${NC}"
echo "  To activate:   source activate_env.sh"
echo "  To deactivate: conda deactivate"
echo "  Or manually:   conda activate ${ENV_NAME}"
echo ""
echo -e "${BLUE}Next Steps:${NC}"
echo "  1. Activate environment: source activate_env.sh"
echo "  2. Run quick test: ./scripts/quick_test.sh"
echo "  3. Start MCP services: docker compose up -d"
echo "  4. Run pipeline: python -m orchestrator.cli run configs/run.example.yaml"
echo ""
echo -e "${YELLOW}Note:${NC} The environment is now active in this shell session."
echo "For new terminal sessions, run: source activate_env.sh"
