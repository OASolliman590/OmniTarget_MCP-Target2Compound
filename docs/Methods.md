# Methods

## Comparators (ChEMBL)
- Filter by `min_pchembl` (config).
- Summarize set quality: size, median pChEMBL, and assay consistency score (0–1).
- Gate targets with insufficient comparators or inconsistent assays.

## Similarity
- GeminiMol embeddings with cosine similarity (adapter interface provided). When the model is not available, no embeddings are fabricated.
- RDKit Morgan fingerprints (Tanimoto) used as a cross-check or fallback.
- Aggregate features per query: `cosine_max`, `cosine_mean_top5`, `tanimoto_max`, `n_cos_ge_tau`.

## QSAR (Ouroboros)
- Optional; attempted only when labeled data suffice (>= 20 actives).
- Adapter raises if infeasible; pipeline treats QSAR as disabled (weight=0).

## (Removed) Sequence-based predictor
Removed: no sequence–ligand affinity predictor is included in the default pipeline.

## Pharmacophore
- RDKit FeatureFactory + ETKDG conformers. Computes a simple overlap score based on feature family proximity.

## Fusion
- Z-score each component per run.
- Combined score: `S = w_sim*z(sim) + w_qsar*z(qsar) + w_dta*z(dta) + w_ph4*z(ph4) + w_dock*z(docking) + w_evd*evidence`.
- Weights are config-driven; set QSAR weight to 0 if disabled.

## Docking (gated and optional)
- Receptor/ligand to PDBQT via Open Babel (`obabel`).
- Docking via AutoDock Vina with a pocket box derived from a validated co-crystal ligand.
- Validation gates: CCD additive filter, SIFTS mapping (PDBe API), optional BioLiP/PDBe ligand selection (auto-validate).
- Docking scores contribute to fusion as `docking_score = -Vina_energy` (more negative energy => higher score).
- For structure preparation and pocket analysis, optionally integrate PDB Prepare Wizard for chain selection, HETATM handling, and PLIP-based pocket characterization.
