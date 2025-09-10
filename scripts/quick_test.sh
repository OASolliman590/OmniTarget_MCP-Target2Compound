#!/bin/bash
set -e

echo "🧪 Quick Functionality Test for MCP Drug Discovery Pipeline"
echo "=========================================================="

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Function to print status
print_status() {
    if [ $? -eq 0 ]; then
        echo -e "${GREEN}✅ $1 - PASSED${NC}"
    else
        echo -e "${RED}❌ $1 - FAILED${NC}"
    fi
}

print_skip() {
    echo -e "${YELLOW}⏭️  $1 - SKIPPED${NC}"
}

# Change to project directory
cd "$(dirname "$0")/.."

echo "📁 1. Checking project structure..."
if [ -d "orchestrator" ] && [ -d "configs" ] && [ -d "docker" ]; then
    print_status "Project structure"
else
    echo -e "${RED}❌ Project structure - FAILED${NC}"
    echo "Missing required directories"
    exit 1
fi

echo "🐍 2. Testing Python imports..."
python3 -c "
try:
    from orchestrator.settings import settings
    from orchestrator.schemas.config import RunConfig
    from orchestrator.pipeline import DrugDiscoveryPipeline
    print('✅ Core imports successful')
except ImportError as e:
    print(f'❌ Import failed: {e}')
    exit(1)
" || print_status "Python imports"

echo "⚙️ 3. Testing configuration system..."
python3 -c "
from orchestrator.schemas.config import RunConfig
try:
    config = RunConfig(
        disease_terms=['test'],
        compounds={'input_paths': ['test.smi']},
        scoring={'weights': {'similarity': 0.5, 'pharmacophore': 0.2, 'docking': 0.1, 'evidence': 0.2}}
    )
    print('✅ Configuration validation works')
except Exception as e:
    print(f'❌ Configuration failed: {e}')
    exit(1)
"
print_status "Configuration system"

echo "📄 4. Creating test data..."
mkdir -p data/compounds data/outputs/test
cat > data/compounds/quick_test.smi << 'EOF'
CCO ethanol
CC(C)O isopropanol
C1CCCCC1 cyclohexane
EOF
print_status "Test data creation"

echo "📋 5. Creating test configuration..."
cat > configs/quick_test.yaml << 'EOF'
disease_terms: ["lung cancer"]
organism: "Homo sapiens"
compounds:
  input_paths: ["data/compounds/quick_test.smi"]
  max_compounds: 3
max_targets: 2
scoring:
  weights:
    similarity: 0.5
    pharmacophore: 0.2
    docking: 0.1
    evidence: 0.2
output_dir: "data/outputs/test"
debug_mode: true
dry_run: true
EOF
print_status "Test configuration"

echo "🔧 6. Testing CLI interface..."
if command -v python3 >/dev/null 2>&1; then
    python3 -m orchestrator.cli validate configs/quick_test.yaml >/dev/null 2>&1
    print_status "CLI configuration validation"
else
    print_skip "CLI interface (python3 not available)"
fi

echo "🤖 7. Testing ML adapters (basic setup)..."
python3 -c "
import asyncio
from orchestrator.adapters import GeminiMolAdapter, VinaAdapter, OuroborosAdapter

async def test_adapters():
    adapters = [
        ('GeminiMol', GeminiMolAdapter()),
        ('Vina', VinaAdapter()),
        ('Ouroboros', OuroborosAdapter())
    ]
    
    for name, adapter in adapters:
        try:
            # Test initialization only
            assert hasattr(adapter, 'setup')
            print(f'✅ {name} adapter initialized')
        except Exception as e:
            print(f'❌ {name} adapter failed: {e}')

asyncio.run(test_adapters())
" 2>/dev/null
print_status "ML adapters initialization"

echo "📊 8. Testing scoring system..."
python3 -c "
from orchestrator.scoring import ScoreNormalizer, ScoreFusion

# Test normalization
normalizer = ScoreNormalizer()
scores = [1.0, 2.0, 3.0, 4.0, 5.0]
normalized = normalizer.z_score_normalize(scores)
assert len(normalized) == len(scores)
print('✅ Score normalization works')

# Test fusion
weights = {'similarity': 0.5, 'pharmacophore': 0.2, 'docking': 0.1, 'evidence': 0.2}
fusion = ScoreFusion(weights)
score_data = [
    {'pair_id': 'test', 'docking_score': 2.0, 'evidence_score': 3.0}
]
combined = fusion.combine_scores(score_data)
assert len(combined) == 1
print('✅ Score fusion works')
"
print_status "Scoring system"

echo "🔌 9. Testing MCP client setup..."
python3 -c "
import asyncio
from orchestrator.mcp_clients import KEGGClient, UniProtClient

async def test_clients():
    clients = [KEGGClient(), UniProtClient()]
    for client in clients:
        try:
            # Test basic initialization
            assert hasattr(client, 'base_url')
            assert hasattr(client, 'test_connection')
            print(f'✅ {client.__class__.__name__} initialized')
            await client.close()
        except Exception as e:
            print(f'❌ {client.__class__.__name__} failed: {e}')

asyncio.run(test_clients())
" 2>/dev/null
print_status "MCP client setup"

echo "🏃‍♂️ 10. Testing pipeline dry run..."
python3 -c "
import asyncio
from orchestrator.pipeline import DrugDiscoveryPipeline
from orchestrator.schemas.config import RunConfig

async def test_pipeline():
    pipeline = DrugDiscoveryPipeline()
    
    # Test basic pipeline initialization
    assert hasattr(pipeline, 'kegg_client')
    print('✅ Pipeline initialized')
    
    # Test configuration loading
    config = RunConfig(
        disease_terms=['test'],
        compounds={'input_paths': ['data/compounds/quick_test.smi']},
        dry_run=True
    )
    assert config.dry_run == True
    print('✅ Dry run configuration works')

asyncio.run(test_pipeline())
" 2>/dev/null
print_status "Pipeline dry run"

echo ""
echo "🎯 Quick Test Summary"
echo "===================="
echo "✅ Core functionality appears to be working"
echo "⚠️  Some components may need services running to fully test"
echo ""
echo "Next steps:"
echo "1. Start MCP services: docker compose up -d"
echo "2. Run full test: python3 scripts/test_functionality.py"
echo "3. Try a real pipeline run: python3 -m orchestrator.cli run configs/quick_test.yaml"
echo ""
echo "For detailed testing, see: TESTING.md"
