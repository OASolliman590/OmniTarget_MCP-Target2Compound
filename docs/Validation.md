# Validation

Report the following checks:
- Ablations: similarity-only vs pharmacophore-only vs fused.
- Comparator gating: skip targets without sufficient comparators or consistent assays; log reasons.
- External checks (when available): known actives should score above decoys.

Recommended artifacts:
- CSVs for each ablation run.
- `manifest.json` capturing thresholds, weights, and counts.
