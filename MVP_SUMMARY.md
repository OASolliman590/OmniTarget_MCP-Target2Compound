# MVP Summary - MCP Drug Discovery Pipeline

## 🎯 **WHAT IS THE MVP?**

The Minimum Viable Product is a **complete non-docking drug discovery workflow** that:

1. Takes disease terms as input
2. Discovers relevant target proteins from pathways
3. Loads your compound library
4. Scores compound-target pairs using similarity and evidence
5. Ranks and outputs prioritized hits with provenance

**Status**: ✅ **FULLY OPERATIONAL** (tested end-to-end with real data)

---

## 🔧 **MVP PIPELINE STAGES**

### 1. **Ligand Preparation**
- Load SMILES/SDF files
- Standardize structures (pH 7.4)
- Remove salts, select largest fragment
- Generate 3D conformers
- Export to SDF + SMILES

### 2. **Target Discovery**
- Query KEGG pathways for disease terms
- Query Reactome pathways for enrichment
- Extract unique gene targets
- Limit to `max_targets` if specified

### 3. **ChEMBL Comparators**
- For each target, fetch ChEMBL bioactivity data (IC50)
- Filter by `min_pchembl` threshold (default: 6.0)
- Build quality summary (median pChEMBL, assay consistency)
- Apply gates: `min_comparators`, `require_assay_consistency`
- Cache results to avoid repeated API calls

### 4. **Similarity Scoring**
- Compute Morgan 2D fingerprints (RDKit Tanimoto)
- Compute top-k neighbors against ChEMBL actives
- Extract features: `tanimoto_max`, `tanimoto_mean_topk`

### 5. **E3FP 3D Similarity (Optional)**
- Generate 3D conformer-based fingerprints
- Compute 3D Tanimoto against ChEMBL actives
- Extract: `e3fp_tanimoto_max`, `e3fp_tanimoto_mean_topk`

### 6. **Pharmacophore Scoring (Optional)**
- RDKit feature pharmacophore (donors, acceptors, aromatics)
- Best overlap score against comparators

### 7. **Score Fusion**
- Normalize all scores (z-score or min-max)
- Apply configurable weights:
  - `tanimoto_max`: 30% (default)
  - `e3fp_tanimoto_max`: 10% (if enabled)
  - `ph4_best`: 25% (if enabled)
  - `evidence_strength`: 10%
  - Custom weights via config
- Compute `fused_score` and rank

### 8. **Output Generation**
- Write `results.csv` with ranked compound-target pairs
- Write `manifest.json` with full provenance
- Optional: docking results, visualization reports

---

## 📊 **MVP OUTPUT**

### `results.csv` Columns

| Column | Description |
|--------|-------------|
| `compound_id` | Compound identifier |
| `compound_smiles` | Canonical SMILES |
| `compound_inchikey` | InChI key for deduplication |
| `compound_source` | Source file path |
| `target_uniprot` | UniProt accession |
| `target_chembl` | ChEMBL target ID |
| `target_gene` | Gene symbol |
| `evidence_strength` | ChEMBL comparator evidence (0-1) |
| `n_comparators` | Number of ChEMBL actives |
| `median_pchembl` | Median pChEMBL value of actives |
| `assay_consistency` | Assay type consistency (0-1) |
| `gate_reason` | Why target was rejected (if any) |
| `tanimoto_max` | Best 2D Tanimoto similarity |
| `e3fp_tanimoto_max` | Best 3D Tanimoto similarity |
| `e3fp_tanimoto_mean_topk` | Mean top-k 3D similarity |
| `gm_cosine_max` | GeminiMol cosine (if enabled) |
| `gm_cosine_mean_topk` | GeminiMol mean top-k |
| `gm_profile_score` | PharmProfiler target score |
| `ph4_best` | Best pharmacophore overlap |
| `qsar_score` | QSAR prediction (if available) |
| `docking_score` | Vina docking (if enabled) |
| `fused_score` | Final weighted score |
| `rank` | Rank index (1 = best) |

### `manifest.json` Contents
- Run ID, timestamp, config hash
- Stage results (n_targets, n_comparators, n_results)
- Output file paths with SHA256 hashes
- E3FP/GeminiMol/Ouroboros summaries
- Duration and performance metrics

---

## 🚀 **HOW TO RUN THE MVP**

### Step 1: Setup
```bash
# Activate environment
conda activate mcp-drug-discovery

# Start MCP services
docker compose -f docker/docker-compose.yaml up -d

# Verify services
python -c "
import asyncio
from orchestrator.mcp_clients import *
async def test(): 
    clients = [KEGGClient(), ReactomeClient(), UniProtClient(), ChEMBLClient()]
    for client in clients: 
        print(f'{client.__class__.__name__}: {\"✅\" if await client.test_connection() else \"❌\"}')
        await client.close()
asyncio.run(test())
"
```

### Step 2: Prepare Data
```bash
# Create compound file (SMILES format)
cat > data/compounds/my_compounds.smi << EOF
CCO ethanol
CC(C)O isopropanol
C1CCCCC1 cyclohexane
EOF
```

### Step 3: Create Config
```bash
# Generate example config
python -m orchestrator.cli example_nd_config --output-path my_run.yaml

# Edit config
vi my_run.yaml
```

Example minimal config:
```yaml
disease_terms:
  - lung cancer
  - breast cancer

organism: "Homo sapiens"
max_targets: 20

ligand_prep:
  input_paths:
    - "data/compounds/my_compounds.smi"
  ph: 7.4
  conformers:
    enable: true

chembl:
  min_pchembl: 6.0
  min_comparators: 5
  require_assay_consistency: true
  limit: 200

similarity:
  top_k: 50
  enable_e3fp: true  # Enable 3D fingerprints
  e3fp:
    radius: 2
    shells: 5
    top_k: 25

pharmacophore:
  method: "rdkit_features"

scoring:
  weights:
    tanimoto_max: 0.30
    e3fp_tanimoto_max: 0.10
    ph4_best: 0.25
    evidence_strength: 0.10

output_dir: "data/outputs"
```

### Step 4: Run Pipeline
```bash
# Run with your config
python -m orchestrator.cli run_nd my_run.yaml

# Or dry-run to validate
python -m orchestrator.cli run_nd my_run.yaml --dry-run
```

### Step 5: Check Results
```bash
# Find your run directory
ls -la data/outputs/

# View results
head -20 data/outputs/run_*/results.csv

# View manifest
cat data/outputs/run_*/manifest.json | jq .
```

---

## 📈 **EXPECTED PERFORMANCE**

### Timing (Typical)
- **Small dataset** (10 compounds, 5 targets): 2-5 minutes
- **Medium dataset** (100 compounds, 20 targets): 10-20 minutes
- **Large dataset** (1000 compounds, 50 targets): 30-60 minutes

### Bottlenecks
- ChEMBL API calls (cached after first run)
- 3D conformer generation (if enabled)
- E3FP 3D fingerprints (computationally intensive)

### Optimization Tips
- Use caching (enabled by default)
- Limit `max_targets` and `max_compounds`
- Disable 3D features for initial screening
- Use `dry_run: true` to validate config without running

---

## 🎯 **SUCCESS CRITERIA**

### ✅ MVP is successful if:
1. Pipeline completes without errors
2. `results.csv` contains ranked compound-target pairs
3. All scores are within [0, 1] range (normalized)
4. `manifest.json` includes all stage results
5. No synthetic or mock data in outputs
6. Provenance is complete (input hashes, output hashes)

### 🚫 MVP fails if:
- MCP services are down (pre-check with status command)
- No ChEMBL data for any target (increase target pool)
- No compounds pass standardization (check SMILES validity)
- All targets fail quality gates (adjust thresholds)

---

## 🔧 **TROUBLESHOOTING**

### Issue: "No targets discovered"
**Solution**: 
- Check disease terms are specific (e.g., "lung cancer" not "cancer")
- Increase `max_targets` limit
- Check KEGG/Reactome MCP services are running

### Issue: "No ChEMBL comparators found"
**Solution**:
- Lower `min_pchembl` threshold (try 5.5 or 5.0)
- Lower `min_comparators` requirement (try 3)
- Check target ChEMBL IDs are valid
- Disable `require_assay_consistency` temporarily

### Issue: "Pipeline takes too long"
**Solution**:
- Disable E3FP: `enable_e3fp: false`
- Reduce `max_targets` and `max_compounds`
- Use cached ChEMBL data (automatic)
- Disable conformer generation initially

### Issue: "MCP services not responding"
**Solution**:
```bash
# Check service status
docker compose ps

# View logs
docker compose logs kegg-mcp
docker compose logs reactome-mcp

# Restart all services
docker compose down
docker compose up -d

# Wait 30 seconds for startup
sleep 30

# Test connections
python -m orchestrator.cli status
```

---

## 📚 **WHAT'S BEYOND THE MVP?**

### Enhanced Features (Integrated, Need Setup)
- **GeminiMol**: Deep learning embeddings (1024-D vectors)
- **Ouroboros**: Chemical generation and optimization
- **PharmProfiler**: Ligand-based target identification
- **Docking**: Structure-based scoring with Vina

### Future Extensions
- QSAR models per target
- Multi-objective optimization
- Active learning loops
- Virtual screening at scale
- Multi-target drug design

---

## ✅ **CONCLUSION**

**The MVP is a production-ready, end-to-end drug discovery pipeline.**

**You can run it right now with:**
- Your disease terms (any KEGG/Reactome pathway-related disease)
- Your compound library (SMILES or SDF format)
- Real MCP service data (no mocks)
- Configurable scoring weights
- Full provenance tracking

**Next steps:**
1. Prepare your compound library
2. Choose your disease terms
3. Run the pipeline
4. Analyze the results
5. Iterate on parameters

**The pipeline will output ranked compound-target pairs backed by real ChEMBL evidence and molecular similarity.**

---

**For detailed documentation, see:**
- Full status: `STATUS.md`
- Setup guide: `docs/Setup_Guide.md`
- Methods: `docs/Methods.md`
- Testing: `TESTING.md`

