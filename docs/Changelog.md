# Changelog

## 0.4.1 (Models placed, env + scripts)
- Cloned upstream Ouroboros into `third_party/Ouroboros` and normalized adapter to detect Chemical* under `ouroboros/`.
- Placed GeminiMol weights in `third_party/GeminiMol/models/GeminiMol` and Ouroboros model in `third_party/Ouroboros/models/Ouroboros_M1c` (flattened).
- Added helpers: `scripts/fetch_geminimol_models.sh` (gdown/HF) and `scripts/fetch_ouroboros_models.sh` (placement guide).
- Added `.env` entries for GeminiMol (`GeminiMol`, `geminimol_app`, `geminimol_lib`, `geminimol_data`) and Ouroboros (`Ouroboros`, `ouroboros_app`, `ouroboros_lib`, `ouroboros_dataset`).
- Enabled GeminiMol and Ouroboros in `configs/run.nd.example.yaml` by default.
- Improved soft-failure: GeminiMol PharmProfiler path resolution via env and multiple locations.

## 0.4.0 (GeminiMol & Ouroboros Integration)
- **Added GeminiMol integration**: Molecular embeddings and PharmProfiler for ligand-based target identification
- **Added Ouroboros Chemical modes**: ChemicalCheck, ChemicalExploration, ChemicalMigration, and ChemicalFusion for directed chemical evolution
- **Enhanced similarity channels**: Now supports Morgan 2D, E3FP 3D, and GeminiMol embeddings with configurable weights
- **Updated fusion scoring**: Added GeminiMol weight (default 25%) and rebalanced weights (Morgan 30%, E3FP 10%, GeminiMol 25%, pharmacophore 25%, evidence 10%)
- **New CLI flags**: Added --enable-geminimol and --enable-ouroboros flags for optional channel activation
- **Enhanced results schema**: Added gm_cosine_max, gm_cosine_mean_topk, gm_profile_score columns
- **Third-party integration**: Added third_party directory structure for GeminiMol and Ouroboros repositories
- **Comprehensive documentation**: Updated Methods.md, Validation.md with new channels and ablation strategies
- **Soft failure design**: All new channels fail gracefully when dependencies are unavailable (no fabricated data)
- **Setup guide**: Added comprehensive Setup_Guide.md with current readiness status and setup instructions
- **Current status**: E3FP ready, GeminiMol and Ouroboros require repository setup and model files

## 0.3.0 (E3FP Integration & Enhanced Reporting)
- **Added E3FP 3D fingerprints**: Extended 3-Dimensional FingerPrint integration for 3D molecular similarity analysis
- **Enhanced configuration**: Relaxed 2D similarity thresholds (min_tanimoto_morgan: 0.25) and pharmacophore support (min_ph4_support: 0.40) for discovery runs
- **Enabled conformers by default**: 3D conformer generation now enabled in ligand preparation for pharmacophore and E3FP analysis
- **Exposed comparator diagnostics**: Added n_comparators, median_pchembl, assay_consistency, and gate_reason columns to results.csv
- **Updated fusion scoring**: Added E3FP weight (default 10%) to score fusion with configurable weights
- **Enhanced manifest logging**: E3FP parameters and computation status now logged in manifest.json
- **Improved documentation**: Updated Methods.md, Validation.md with E3FP details and new ablation strategies

## 0.2.1 (Remove sequence-based predictor)
- Removed the previous sequence-based affinity predictor and its configuration; default pipeline remains non-docking (similarity + pharmacophore + evidence). Optional docking stays gated and off by default.

## 0.2.0 (Non-Docking Default)
- Added non-docking pipeline (`run_nd`) with ChEMBL comparators, GeminiMol/Ouroboros stubs, and RDKit pharmacophore.
- Implemented ligand IO supporting .mol/.mol2/.sdf/.smi with 3D conformers.
- Added provenance manifest and docs scaffolding.
