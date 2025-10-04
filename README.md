# MCP Drug Discovery Pipeline

A comprehensive drug discovery pipeline that orchestrates target discovery, molecular docking, affinity prediction, and bioactivity validation using MCP servers and state-of-the-art ML models.

## 🎯 Overview

This pipeline integrates multiple biological databases and ML models to identify promising drug-target pairs:

- **Target Discovery**: KEGG pathways, Reactome enrichment
- **Target Characterization**: Protein Atlas expression, STRING PPI, UniProt sequences, PDB structures  
- **Molecular Representation**: GeminiMol embeddings, Ouroboros features
- **Ligand Similarity & Features**: GeminiMol/Ouroboros, RDKit fingerprints
- **Molecular Docking**: AutoDock Vina structure-based screening
- **Evidence Validation**: ChEMBL bioactivity data

Non-docking target prediction is the default workflow. Comparator evidence (ChEMBL) plus similarity/ML features drive ranking. Docking is optional and gated: it runs only when a co-crystal pocket is validated (BioLiP/PDBe + CCD additive filter + SIFTS mapping) and required tools are available (Open Babel + AutoDock Vina + RDKit). See `docs/Documentation.md` and `docs/Methods.md` for details.

**Key Principles:**
- ✅ **Live API integration** - real data from MCP servers, no mocks or fixtures
- ✅ **Reproducible** - complete provenance tracking and versioning
- ✅ **Scalable** - async processing, caching, and containerized services

## 📊 Current Status: **ENHANCED WITH E3FP, GEMINIMOL & OUROBOROS** ✅

### ✅ **Ready Components**
- **E3FP 3D Fingerprints**: Fully installed and operational (v1.2.7)
- **Basic Pipeline**: Non-docking with Morgan 2D, pharmacophore, evidence scoring
- **MCP Services**: All 6 services running and responding correctly
- **Soft Failure Design**: Graceful handling of missing dependencies

### ✅ **Model & Config Status**
- **GeminiMol**: Repo present, model weights placed, env vars added; example config enables it
- **Ouroboros**: Repo present, model placed, Chemical* scripts detected; example config enables it

### 🎯 **Enhanced Capabilities**
- **Multi-channel Similarity**: Morgan 2D + E3FP 3D + GeminiMol embeddings
- **Target Identification**: PharmProfiler integration for ligand-based prediction
- **Chemical Generation**: Ouroboros Chemical modes for directed evolution
- **Configurable Weights**: Balanced fusion scoring across all channels

📋 **For setup instructions and model placement, see [docs/Setup_Guide.md](docs/Setup_Guide.md)**

## 🏗️ Architecture

```
┌─────────────────┐    ┌─────────────────┐    ┌─────────────────┐
│   MCP Servers   │    │   Orchestrator  │    │   ML Adapters   │
│                 │    │                 │    │                 │
│ • KEGG          │────│ • Pipeline      │────│ • Similarity    │
│ • Reactome      │    │ • FastAPI       │    │ • GeminiMol     │
│ • ProteinAtlas  │    │ • CLI           │    │ • AutoDock Vina │
│ • STRING        │    │ • Scoring       │    │ • Ouroboros     │
│ • UniProt       │    │ • Caching       │    │                 │
│ • PDB           │    │                 │    │                 │
│ • ChEMBL        │    │                 │    │                 │
└─────────────────┘    └─────────────────┘    └─────────────────┘
```

## 🚀 Quick Start

### Prerequisites

- Docker & Docker Compose
- Python 3.11+ (or Miniconda/Anaconda)
- Git with submodules support

### 1. Clone and Setup

**Option A: Conda Environment (Recommended)**
```bash
git clone <repository-url>
cd mcp-drug-discovery

# Automated conda setup
chmod +x conda_env_setup.sh
./conda_env_setup.sh

# Activate environment
source activate_env.sh
```

### Non-docking quick run

```
python -m orchestrator.cli example_nd_config --output-path run.nd.yaml
python -m orchestrator.cli run_nd run.nd.yaml
```

### Optional docking (gated)

1) Ensure dependencies:
   - RDKit import works
   - `obabel` on PATH (Open Babel)
   - AutoDock Vina binary configured in settings (`models.vina_binary_path`)

2) Edit `run.nd.yaml`:

```
docking:
  enabled: true
  receptor_pdb_path: "data/structures/target.pdb"  # or .cif
  cocrystal:
    ligand_resname: "LIG"    # non-additive residue ID
    chain_id: "A"            # optional
    uniprot_id: "P12345"     # optional; enforces SIFTS mapping when pdb_id present
  pdb_id: "1ABC"             # optional; used if auto_validate is true
  auto_validate: true         # optionally fetch/select ligand via PDBe
  box_margin: 12.0

scoring:
  weights:
    docking: 0.10             # include docking in fused score
```

3) Run: `python -m orchestrator.cli run run.nd.yaml`

Outputs will include `docking_results.csv` and `results.csv` with a `docking_score` column.

Tip: For receptor preparation and pocket analysis, consider evaluating the PDB Prepare Wizard (see docs below) to improve cleanup and pocket quality.

### Reporting & Visualization

Generate a static HTML report after a run (non-docking default):

```
python scripts/generate_report.py --results data/outputs/run_*/results.csv --manifest data/outputs/run_*/manifest.json
```

Or enable reporting in the config:

```
visualization:
  enabled: true
  title: "My Pipeline Report"
```

Compare two runs:

```
python scripts/compare_reports.py \
  --run1-results data/outputs/run_A/results.csv --run1-manifest data/outputs/run_A/manifest.json \
  --run2-results data/outputs/run_B/results.csv --run2-manifest data/outputs/run_B/manifest.json \
  --out data/outputs/compare --title "Run A vs Run B"
```

### 2. Install ML Dependencies

**AutoDock Vina (Already Built)**
```bash
# Vina is already built and available as 'vina_custom'
# Binary location: $CONDA_PREFIX/bin/vina_custom
# Version: v1.2.7-20-g93cdc3d-mod
```

```

**Option B: Docker Only**
```bash
git clone <repository-url>
cd mcp-drug-discovery
```

**Option C: Pip Installation**
```bash
git clone <repository-url>
cd mcp-drug-discovery
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate
pip install -r requirements.txt
```

For detailed installation instructions, see [INSTALLATION.md](INSTALLATION.md).

```bash
# Initialize submodules for third-party tools
git submodule update --init --recursive

# Add MCP servers as submodules
git submodule add https://github.com/Augmented-Nature/KEGG-MCP-Server.git services/kegg-mcp
git submodule add https://github.com/Augmented-Nature/Reactome-MCP-Server.git services/reactome-mcp
git submodule add https://github.com/Augmented-Nature/ProteinAtlas-MCP-Server.git services/proteinatlas-mcp
git submodule add https://github.com/Augmented-Nature/STRING-db-MCP-Server.git services/string-mcp
git submodule add https://github.com/TakumiY235/uniprot-mcp-server.git services/uniprot-mcp
git submodule add https://github.com/Augmented-Nature/PDB-MCP-Server.git services/pdb-mcp
git submodule add https://github.com/Augmented-Nature/ChEMBL-MCP-Server.git services/chembl-mcp

# Add ML model repositories
git submodule add https://github.com/Wang-Lin-boop/GeminiMol.git third_party/GeminiMol
git submodule add https://github.com/ccsb-scripps/AutoDock-Vina.git third_party/autodock-vina
```

### 3. Prepare Your Data

```bash
# Create a compound file (SMILES format)
mkdir -p data/compounds
cat > data/compounds/example.smi << EOF
CCO ethanol
CC(C)O isopropanol
C1CCCCC1 cyclohexane
EOF
```

### 4. Configure Environment

```bash
# An `.env` is provided with sensible defaults incl. GeminiMol/Ouroboros paths
# Edit it to adjust paths or add API keys
${EDITOR:-vi} .env
```

### 5. Start MCP Services and Run Pipeline

```bash
# Start all MCP servers
docker compose -f docker/docker-compose.yaml up -d

# Test MCP connections
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

# Run functionality tests
python scripts/test_functionality.py

# Run pipeline with sertraline compounds
python -m orchestrator.cli run configs/test_sertraline.yaml

# Run dry-run validation
python -m orchestrator.cli run configs/test_sertraline.yaml --dry-run
```

## 📋 Configuration

The pipeline is configured via YAML files. See `configs/run.example.yaml` for a complete example:

```yaml
disease_terms: ["lung cancer"]
organism: "Homo sapiens"

compounds:
  input_paths:
    - "data/compounds/my_compounds.smi"
  max_compounds: 1000

scoring:
  weights:
    similarity: 0.5
    pharmacophore: 0.2
    docking: 0.1
    docking: 0.3
    evidence: 0.1

output_dir: "data/outputs"
```

### Configuration Sections

- **`disease_terms`**: Disease names for target discovery
- **`ppi`**: Protein-protein interaction filtering
- **`expression`**: Tissue expression filtering  
- **`structures`**: 3D structure requirements
- **`compounds`**: Input compound files and filters
- **`scoring`**: Scoring weights and parameters

## 🔧 API Usage

The pipeline also provides a REST API:

```bash
# Start API server
docker compose up orchestrator

# Submit a run
curl -X POST "http://localhost:8000/runs" \
  -H "Content-Type: application/json" \
  -d @configs/api_request.json

# Check run status
curl "http://localhost:8000/runs/{run_id}"

# Download results
curl "http://localhost:8000/runs/{run_id}/results"
```

## 📊 Output

Each pipeline run generates:

```
data/outputs/run_20240115_143022_abc123/
├── results.csv              # Ranked compound-target pairs
├── manifest.json            # Complete provenance record
├── targets.json             # Discovered and characterized targets
├── compounds.json           # Processed compound library
├── docking_results.csv      # Molecular docking scores (if docking enabled and validated)
├── evidence_records.csv     # Bioactivity evidence
└── plots/                   # Visualization plots
    ├── score_distribution.html
    ├── target_network.html
    └── compound_similarity.html
```

### Results Format

The main `results.csv` contains (non-docking default):

| Column | Description |
|--------|-------------|
| `compound_id` | Compound identifier |
| `compound_smiles` | Canonical SMILES |
| `target_id` | Target identifier |
| `target_gene` | Gene symbol |
| `cosine_max` | Similarity (GeminiMol cosine or fallback) |
| `tanimoto_max` | RDKit Morgan Tanimoto (fallback) |
| `ph4_best` | Pharmacophore best overlap score |
| `evidence_strength` | Comparator evidence strength (0–1) |
| `docking_score` | Normalized docking score (if enabled) |
| `fused_score` | Final integrated score |
| `rank` | Rank index |

## 🧪 Testing

The pipeline includes integration tests that use real data:

```bash
# Run tests (requires test compounds in data/compounds/test_set.smi)
docker compose exec orchestrator python -m pytest orchestrator/tests/

# Run specific test
docker compose exec orchestrator python -m pytest orchestrator/tests/test_integration.py::test_full_pipeline
```

**Note**: Tests require real compound files - no synthetic or mock data is used.

### Environment check

Run `python scripts/check_env.py` to verify RDKit, Open Babel, Vina, and optional PDBe connectivity before enabling docking.

## 🔄 Development

### Adding New MCP Servers

1. Add as git submodule in `services/`
2. Update `docker-compose.yaml` with new service
3. Create client in `orchestrator/mcp_clients/`
4. Update pipeline to use new client

### Adding New ML Models

1. Create adapter in `orchestrator/adapters/`
2. Implement required interface methods
3. Update pipeline to integrate new scores
4. Add configuration options

### Extending Scoring

1. Add new score types to `schemas/scoring.py`
2. Update normalization in `scoring/normalize.py`
3. Modify fusion weights in `scoring/fuse.py`

## 📁 Project Structure

```
mcp-drug-discovery/
├── orchestrator/           # Main Python application
│   ├── mcp_clients/       # MCP server clients
│   ├── adapters/          # ML model adapters
│   ├── scoring/           # Score normalization and fusion
│   ├── schemas/           # Pydantic data models
│   ├── tests/             # Integration tests
│   ├── cli.py             # Command-line interface
│   ├── api.py             # FastAPI application
│   ├── pipeline.py        # Main orchestration logic
│   └── settings.py        # Configuration management
├── services/              # MCP server submodules
├── third_party/           # ML model submodules
├── docker/                # Docker configuration
├── configs/               # Configuration examples
├── data/                  # User data directory
│   ├── compounds/         # Input compound files
│   └── outputs/           # Pipeline results
└── README.md
```

## 🛠️ Troubleshooting

### Common Issues

**MCP Server Connection Failures**
```bash
# Check server status
docker compose ps
docker compose logs kegg-mcp

# Restart specific service
docker compose restart kegg-mcp
```

**Memory Issues with Large Compound Sets**
```bash
# Increase Docker memory limits
# Or reduce compounds.max_compounds in config
```

**Missing Structures for Docking**
```bash
# Check structure availability
# Set structures.require_structure: false to skip
```

### Performance Optimization

- **Caching**: Enable `cache_results: true` in config
- **Parallelization**: Increase `parallel_jobs` for CPU-bound tasks
- **Filtering**: Use tissue/expression filters to reduce target count
- **Limits**: Set `max_compounds` and `max_targets` for large datasets

## 📚 References

1. **KEGG**: [github.com/Augmented-Nature/KEGG-MCP-Server](https://github.com/Augmented-Nature/KEGG-MCP-Server)
2. **Reactome**: [github.com/Augmented-Nature/Reactome-MCP-Server](https://github.com/Augmented-Nature/Reactome-MCP-Server)
3. **GeminiMol**: [github.com/Wang-Lin-boop/GeminiMol](https://github.com/Wang-Lin-boop/GeminiMol)
4. **AutoDock Vina**: [github.com/ccsb-scripps/AutoDock-Vina](https://github.com/ccsb-scripps/AutoDock-Vina)
5. **PDBe APIs (SIFTS & Ligands)**: https://www.ebi.ac.uk/pdbe/api
6. **PDB Prepare Wizard**: [github.com/OASolliman590/pdb-prepare-wizard](https://github.com/OASolliman590/pdb-prepare-wizard)

## 📄 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests for new functionality
5. Submit a pull request


---

**Built with ❤️ using MCP servers, modern ML models, and reproducible science principles.**
