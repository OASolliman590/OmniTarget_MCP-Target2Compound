# Testing Summary - How to Test Each Functionality

## 📊 Current Testing Status: **FULLY FUNCTIONAL WITH MCP SERVICES** ✅

### ✅ **Completed Tests**
- **Environment Setup**: Conda environment with Python 3.11 and all dependencies
- **AutoDock Vina**: Successfully built from source (v1.2.7-20-g93cdc3d-mod) and tested
- **DeepDTA**: Repository cloned, adapter tested with real sertraline compounds
- **ML Adapters**: GeminiMol, DeepDTA, Vina, Ouroboros all working with real data
- **Real Compounds**: 11 sertraline compounds loaded and validated
- **Scoring System**: Normalization and fusion tested with realistic compound data
- **Configuration**: Test config validated for sertraline compounds
- **MCP Servers**: All 6 services running and responding correctly
- **End-to-End Pipeline**: Fully functional with real compounds and MCP integration

### 🎯 **Fully Operational**
- **MCP Services**: All services operational with live data integration
- **End-to-End Pipeline**: Fully functional with real compounds
- **Target Discovery**: Working with disease terms (depression, anxiety, OCD)

## 🎯 Quick Testing Guide

Here are the different ways to test each component of the MCP Drug Discovery Pipeline:

## 📋 Test Categories

### 1. 🚀 Quick Health Check
```bash
# Run basic functionality check (5 minutes)
./scripts/quick_test.sh

# Check if all core components load
python3 -c "
from orchestrator.pipeline import DrugDiscoveryPipeline
from orchestrator.mcp_clients import KEGGClient
from orchestrator.adapters import DeepDTAAdapter
print('✅ All imports successful')
"
```

### 2. 🔌 MCP Server Testing (All Services Operational)

**Test All MCP Connections:**
```bash
# Start all MCP services
docker compose -f docker/docker-compose.yaml up -d

# Test all MCP connections
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

# Test specific MCP clients
python scripts/test_component.py mcp kegg
python scripts/test_component.py mcp uniprot
python scripts/test_component.py mcp string
"
```

**Check MCP Service Status:**
```bash
# Using CLI
python -m orchestrator.cli status

# Using Docker
docker compose ps
docker compose logs kegg-mcp
```

### 3. 🤖 ML Adapter Testing

**Individual Adapter Testing:**
```bash
# Test ML adapters individually
python scripts/test_component.py adapter deepdta
python scripts/test_component.py adapter geminimol
python scripts/test_component.py adapter vina
python scripts/test_component.py adapter ouroboros
```

**Test Adapter Setup:**
```python
import asyncio
from orchestrator.adapters import DeepDTAAdapter

async def test_deepdta():
    adapter = DeepDTAAdapter()
    setup_ok = await adapter.setup()
    print(f"Setup: {'✅' if setup_ok else '❌'}")
    
    if setup_ok:
        score = await adapter.predict_affinity("CCO", "MVLSPADKTNVKAAW")
        print(f"Prediction: {score.predicted_affinity:.3f}")

asyncio.run(test_deepdta())
```

### 4. 🔄 Pipeline Stage Testing

**Test Individual Pipeline Stages:**
```bash
# Test compound loading
python scripts/test_component.py pipeline compound_loading

# Test target discovery  
python scripts/test_component.py pipeline target_discovery
```

**Test Full Pipeline Components:**
```python
import asyncio
from orchestrator.pipeline import DrugDiscoveryPipeline
from orchestrator.schemas.config import RunConfig

async def test_pipeline():
    pipeline = DrugDiscoveryPipeline()
    setup_ok = await pipeline.setup()
    print(f"Pipeline setup: {'✅' if setup_ok else '❌'}")
    
    # Test with small dataset
    config = RunConfig(
        disease_terms=["lung cancer"],
        compounds={"input_paths": ["data/compounds/test.smi"]},
        dry_run=True
    )
    
    # Run dry run
    result = await pipeline.run(config)
    print(f"Dry run: {result.run_id}")

asyncio.run(test_pipeline())
```

### 5. 📊 Scoring System Testing

**Test Score Normalization:**
```bash
python scripts/test_component.py scoring normalization
```

**Test Score Fusion:**
```bash
python scripts/test_component.py scoring fusion
```

**Manual Scoring Test:**
```python
from orchestrator.scoring import ScoreNormalizer, ScoreFusion

# Test normalization
normalizer = ScoreNormalizer()
scores = [1.0, 2.5, 3.0, 4.5, 5.0]
normalized = normalizer.z_score_normalize(scores)
print(f"Original: {scores}")
print(f"Normalized: {[f'{x:.2f}' for x in normalized]}")

# Test fusion
weights = {"deepdta": 0.6, "docking": 0.3, "evidence": 0.1}
fusion = ScoreFusion(weights)
score_data = [
    {"pair_id": "test", "deepdta_score": 8.5, "docking_score": -7.2, "evidence_score": 6.0}
]
combined = fusion.combine_scores(score_data)
print(f"Combined score: {combined[0]['combined_score']:.3f}")
```

### 6. ⚙️ Configuration Testing

**Test Configuration Validation:**
```bash
# Validate configuration file
python -m orchestrator.cli validate configs/run.example.yaml

# Test programmatically
python -c "
from orchestrator.schemas.config import RunConfig
config = RunConfig(
    disease_terms=['lung cancer'],
    compounds={'input_paths': ['test.smi']},
    scoring={'weights': {'deepdta': 0.6, 'docking': 0.3, 'evidence': 0.1}}
)
print('✅ Configuration valid')
"
```

### 7. 🌐 API Testing

**Test API Endpoints:**
```bash
# Start API server
docker compose up orchestrator &

# Test health endpoint
curl http://localhost:8000/health | jq '.'

# Test configuration validation
curl -X POST http://localhost:8000/validate \
  -H "Content-Type: application/json" \
  -d '{
    "config": {
      "disease_terms": ["lung cancer"],
      "compounds": {"input_paths": ["data/compounds/test.smi"]},
      "scoring": {"weights": {"deepdta": 0.6, "docking": 0.3, "evidence": 0.1}}
    }
  }' | jq '.'

# Submit pipeline run
curl -X POST http://localhost:8000/runs \
  -H "Content-Type: application/json" \
  -d @configs/api_request.json
```

### 8. 🧪 Integration Testing

**Run Comprehensive Test Suite:**
```bash
# Run all integration tests
python -m pytest orchestrator/tests/ -v

# Run specific test categories
python -m pytest orchestrator/tests/test_mcp_clients.py -v
python -m pytest orchestrator/tests/test_adapters.py -v
python -m pytest orchestrator/tests/test_integration.py -v
```

**Run Custom Test Suite:**
```bash
# Comprehensive functionality testing
python scripts/test_functionality.py

# Interactive demo
python scripts/demo_testing.py
```

### 9. 🏃‍♂️ End-to-End Testing

**Small Dataset Test:**
```bash
# Create test data
echo "CCO ethanol" > data/compounds/test.smi
echo "CC(C)O isopropanol" >> data/compounds/test.smi

# Create test config
cat > configs/test.yaml << EOF
disease_terms: ["lung cancer"]
compounds:
  input_paths: ["data/compounds/test.smi"]
  max_compounds: 2
max_targets: 2
scoring:
  weights: {deepdta: 0.6, docking: 0.3, evidence: 0.1}
output_dir: "data/outputs/test"
EOF

# Run pipeline
python -m orchestrator.cli run configs/test.yaml
```

**Dry Run Test:**
```bash
# Test without actually running expensive computations
python -m orchestrator.cli run configs/test.yaml --dry-run
```

### 10. 🔍 Debugging and Troubleshooting

**Debug Individual Components:**
```bash
# Enable debug logging
export LOG_LEVEL=DEBUG

# Test with debug mode
python -c "
import logging
logging.basicConfig(level=logging.DEBUG)

import asyncio
from orchestrator.mcp_clients import KEGGClient

async def debug_test():
    client = KEGGClient()
    result = await client.test_connection()
    print(f'Result: {result}')
    await client.close()

asyncio.run(debug_test())
"
```

**Check Docker Services:**
```bash
# Check service status
docker compose ps

# View logs
docker compose logs orchestrator
docker compose logs kegg-mcp

# Restart specific service
docker compose restart kegg-mcp

# Check resource usage
docker stats
```

**Performance Testing:**
```bash
# Monitor memory usage
python -c "
import psutil
import time
import subprocess

process = subprocess.Popen(['python', '-m', 'orchestrator.cli', 'run', 'configs/test.yaml'])
peak_memory = 0

while process.poll() is None:
    try:
        proc = psutil.Process(process.pid)
        memory = proc.memory_info().rss / 1024 / 1024  # MB
        peak_memory = max(peak_memory, memory)
        print(f'Memory: {memory:.1f} MB')
        time.sleep(5)
    except psutil.NoSuchProcess:
        break

print(f'Peak memory: {peak_memory:.1f} MB')
"
```

## 📊 Test Results Interpretation

### ✅ Success Indicators
- MCP clients return `True` for `test_connection()`
- ML adapters complete `setup()` without errors
- Pipeline stages process data and return expected objects
- Scoring produces normalized, ranked results
- API endpoints return proper HTTP status codes
- Configuration validation catches invalid inputs

### ⚠️ Expected Warnings
- MCP services may be unavailable (return `False` for connections)
- ML models may not have weights downloaded (setup fails)
- Some tests may be skipped if dependencies missing
- Docker services may take time to start

### ❌ Failure Investigation
- Check Docker service logs: `docker compose logs [service]`
- Verify file paths in configurations
- Ensure Python dependencies installed: `pip install -e .`
- Check network connectivity for MCP services
- Verify sufficient disk space for temporary files

## 🔄 Continuous Testing

**Regular Health Checks:**
```bash
# Daily health check script
#!/bin/bash
echo "$(date): Running health checks..."
./scripts/quick_test.sh > health_check.log 2>&1
if [ $? -eq 0 ]; then
    echo "✅ Health check passed"
else
    echo "❌ Health check failed - see health_check.log"
fi
```

**Automated Testing:**
```bash
# Set up automated testing with cron
# Add to crontab: 0 */6 * * * /path/to/health_check.sh
```

This comprehensive testing approach ensures all components work individually and together, providing confidence in the pipeline's functionality.
