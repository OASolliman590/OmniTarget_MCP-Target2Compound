# Validation

Report the following checks:
- **Ablations**: 
  - 2D-only (Morgan fingerprints)
  - 2D+E3FP (Morgan + 3D fingerprints)
  - GeminiMol-only (embeddings + PharmProfiler)
  - Pharmacophore-only
  - Fused (2D+E3FP+GeminiMol+PH4+evidence)
  - (Optional) Docking ablation when Vina is enabled
- **Comparator gating**: skip targets without sufficient comparators or consistent assays; log reasons.
- **External checks** (when available): known actives should score above decoys.
- **E3FP validation**: verify 3D fingerprint computation when E3FP is installed and SDF files are available.
- **GeminiMol validation**: verify embedding computation and PharmProfiler scoring when GeminiMol is available.
- **Ouroboros validation**: verify Chemical mode execution and output generation when Ouroboros is enabled.

Recommended artifacts:
- CSVs for each ablation run with comparator diagnostics columns.
- `manifest.json` capturing thresholds, weights, counts, and all channel parameters (E3FP, GeminiMol, Ouroboros).
- E3FP computation logs showing fingerprint generation success/failure.
- GeminiMol computation logs showing embedding and PharmProfiler success/failure.
- Ouroboros job execution logs and output CSV files in `outputs/run_*/third_party/`.
