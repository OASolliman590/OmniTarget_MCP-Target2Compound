# Setup Guide for Enhanced Drug Discovery Pipeline

This guide covers the setup requirements for all optional components in the enhanced drug discovery pipeline.

## Current Status

### ✅ Ready Components
- **E3FP 3D Fingerprints**: Fully installed and operational
- **Basic Pipeline**: Non-docking pipeline with Morgan 2D fingerprints, pharmacophore, and evidence scoring

### ✅ Updated Status
- **GeminiMol**: Repo present, model files placed, env vars added
- **Ouroboros**: Repo present, model placed, Chemical* scripts detected

---

## E3FP 3D Fingerprints

### Status: ✅ READY
E3FP is already installed and working correctly.

```bash
# Verify installation
python -c "import e3fp; print(f'E3FP v{e3fp.__version__} installed')"
```

### Usage
```bash
# Enable E3FP in pipeline
python -m orchestrator.cli run-nd config.yaml --enable-e3fp
```

---

## GeminiMol Integration

### Status: ✅ MODELS PRESENT

#### What’s already done
- Repo cloned to `third_party/GeminiMol`
- Models placed in `third_party/GeminiMol/models/GeminiMol` (e.g., `GeminiMol.pt`)
- `.env` updated with `GeminiMol`, `geminimol_app`, `geminimol_lib`, `geminimol_data`

#### Fetch/Place models (if needed)
```bash
# Option A: gdown (Google Drive)
python -m pip install -U gdown
mkdir -p third_party/GeminiMol/models/GeminiMol
cd third_party/GeminiMol/models/GeminiMol
python -m gdown https://drive.google.com/uc?id=1IgpkabylSJ0aIwUjdIWOC5ID51pcucw5 -O GeminiMol.zip
unzip -o GeminiMol.zip -d . && rm -f GeminiMol.zip
[ -d GeminiMol ] && (shopt -s dotglob; mv -f GeminiMol/* .; shopt -u dotglob; rmdir GeminiMol || true)
cd -

# Option B: Hugging Face
python -m pip install -U huggingface_hub
huggingface-cli login
huggingface-cli download AlphaMWang/GeminiMol --repo-type model --local-dir third_party/GeminiMol/models/GeminiMol --local-dir-use-symlinks False

# Helper script
bash scripts/fetch_geminimol_models.sh
```

#### Verify
```bash
ls -la third_party/GeminiMol/models/GeminiMol
```

---

## Ouroboros Chemical Modes

### Status: ✅ MODELS PRESENT

#### What’s already done
- Repo cloned to `third_party/Ouroboros`
- Chemical* scripts available under `third_party/Ouroboros/ouroboros/`
- Model placed and normalized in `third_party/Ouroboros/models/Ouroboros_M1c`
- `.env` updated with `Ouroboros`, `ouroboros_app`, `ouroboros_lib`, `ouroboros_dataset`

#### Place models (if needed)
```bash
mkdir -p third_party/Ouroboros/models/Ouroboros_M1c
# Download from ZhangLab page and unzip here, or use a direct link via curl:
cd third_party/Ouroboros/models
curl -L "<DIRECT_URL_TO_OUROBOROS_MODEL_ZIP>" -o Ouroboros_M1c.zip
unzip -o Ouroboros_M1c.zip -d Ouroboros_M1c && rm -f Ouroboros_M1c.zip
cd -

# Helper scaffold
bash scripts/fetch_ouroboros_models.sh
```

#### Normalize structure (if nested)
```bash
cd third_party/Ouroboros/models/Ouroboros_M1c
shopt -s dotglob; mv -f Ouroboros_V0.1.0/* .; shopt -u dotglob
rmdir Ouroboros_V0.1.0 || true
cd -
```

#### Verify
```bash
ls -la third_party/Ouroboros/models/Ouroboros_M1c
```

#### Example Ouroboros Configuration:
```yaml
ouroboros:
  enabled: true
  repo_path: third_party/Ouroboros
  model_name: "Ouroboros_M1c"
  jobs:
    - type: "check"
      start_smiles: "CCO"
      model_name: "Ouroboros_M1c"
      job_name: "ethanol_check"
    - type: "exploration"
      start_smiles: "CCO"
      mode: "scaffold_hopping"
      model_name: "Ouroboros_M1c"
      job_name: "ethanol_exploration"
      replica_num: 1
      steps: 100
      step_interval: 10
      temperature: 1.0
      learning_rate: 0.01
```

---

## Soft Failure Behavior

All optional components implement **soft failure** - they gracefully handle missing dependencies without crashing the pipeline:

### GeminiMol Soft Failure
- Returns `None` values for all GeminiMol features
- Logs warnings about missing setup
- Pipeline continues with other similarity methods

### Ouroboros Soft Failure
- Skips Chemical mode execution
- Logs warnings about missing setup
- Pipeline continues with baseline target identification

### E3FP Soft Failure
- Returns `None` values for E3FP features
- Logs warnings about missing library
- Pipeline continues with 2D fingerprints

---

## Testing Setup

### Test Current Functionality
```bash
# Test with only E3FP (ready)
python -m orchestrator.cli run-nd configs/run.nd.example.yaml --enable-e3fp

# Test with missing dependencies (should work gracefully)
python -m orchestrator.cli run-nd configs/run.nd.example.yaml --enable-geminimol --enable-ouroboros

# Example config already enables both; simply run
python -m orchestrator.cli run-nd configs/run.nd.example.yaml
```

### Verify Soft Failure
```bash
python -c "
from orchestrator.adapters.geminimol_adapter import compute_geminimol_features
result = compute_geminimol_features(['CCO'], ['CCCO'], ['test'], 'third_party/GeminiMol', 'models/GeminiMol', None, True)
print('GeminiMol soft failure:', result)
"
```

---

## Recommended Workflow

### Phase 1: Current Capabilities
- Use E3FP 3D fingerprints (ready)
- Use Morgan 2D fingerprints (ready)
- Use pharmacophore features (ready)
- Use evidence scoring (ready)

### Phase 2: Enhanced Capabilities (when setup complete)
- Add GeminiMol embeddings for molecular representation
- Add PharmProfiler for target identification
- Add Ouroboros Chemical modes for molecular generation

### Phase 3: Full Integration
- Combine all similarity channels
- Use generated molecules from Ouroboros (if explicitly enabled)
- Comprehensive multi-channel drug discovery

---

## Troubleshooting

### Common Issues

1. **Import Errors**
   - Ensure all dependencies are installed
   - Check Python path and environment

2. **Missing Model Files**
   - Verify model files are in correct directories
   - Check file permissions

3. **Environment Variables**
   - Ensure all required environment variables are set
   - Check variable paths are correct

4. **Script Not Found**
   - Verify Chemical* scripts exist in scripts/ directory
   - Check script permissions and shebang lines

### Getting Help

- Check individual repository READMEs for detailed setup instructions
- Review pipeline logs for specific error messages
- Test individual components before full integration

---

## Next Steps

1. **Immediate**: Use current E3FP-enabled pipeline for enhanced 3D similarity
2. **Short-term**: Set up GeminiMol for molecular embeddings and target identification
3. **Long-term**: Set up Ouroboros for chemical generation and directed evolution

The pipeline is designed to work incrementally - you can add components as they become available without breaking existing functionality.



