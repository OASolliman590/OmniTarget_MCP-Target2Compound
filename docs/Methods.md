# Methods

## Comparators (ChEMBL)
- Filter by `min_pchembl` (config).
- Summarize set quality: size, median pChEMBL, and assay consistency score (0–1).
- Gate targets with insufficient comparators or inconsistent assays.

## Similarity
- **RDKit Morgan fingerprints (Tanimoto)**: 2D molecular similarity using Morgan fingerprints.
- **E3FP 3D fingerprints**: Extended 3-Dimensional FingerPrint method for 3D molecular similarity analysis.
  - Computed from 3D conformers (SDF files) using configurable radius and shells parameters.
  - Tanimoto similarity between E3FP bitvectors provides 3D structural similarity.
  - Features: `e3fp_tanimoto_max`, `e3fp_tanimoto_mean_topk`.
  - Fails soft if E3FP library not installed (no fabricated scores).
- **GeminiMol embeddings**: Molecular embeddings for ligand-based target identification.
  - Computed from SMILES using GeminiMol model (requires model setup and environment variables).
  - Cosine similarity between embeddings provides molecular representation similarity.
  - Features: `gm_cosine_max`, `gm_cosine_mean_topk`, `gm_profile_score`.
  - PharmProfiler integration for target identification scoring.
  - **Status**: Requires GeminiMol repository setup, model files, and environment variables.
  - Fails soft if GeminiMol not available (no fabricated scores).
- Aggregate features per query: `cosine_max`, `cosine_mean_top5`, `tanimoto_max`, `n_cos_ge_tau`.

## Ouroboros Chemical Modes
- **ChemicalCheck**: Reconstruction similarity testing to validate model performance.
  - Tests reconstruction similarity near seed molecules (similarity >0.9 indicates well-modeled local space).
  - Outputs: generation CSV with similarity scores and molecular properties.
- **ChemicalExploration**: Scaffold hopping and optimization modes.
  - Modes: scaffold_hopping, directional_optimization, directional_scaffold_hopping.
  - Configurable hyperparameters: replica_num, steps, step_interval, temperature, learning_rate.
  - Outputs: exploration CSV with generated molecules and properties.
- **ChemicalMigration**: Directed evolution between molecules.
  - Evolves from start_smiles toward ref_smiles using directed chemical evolution.
  - Outputs: migration CSV with evolutionary path and intermediate molecules.
- **ChemicalFusion**: Chimera generation for multi-target compounds.
  - Generates compounds targeting multiple fusion targets (e.g., "AURKA:PI3Kg").
  - Uses probe dataset and configurable fusion temperature.
  - Outputs: fusion CSV with chimera molecules and target profiles.
- **Default behavior**: Disabled by default (`ouroboros.enabled: false`).
- **Integration**: Runs after baseline target identification; does not modify fused ranking unless explicitly enabled.
- **Provenance**: All commands, parameters, and outputs logged in manifest.json.
- **Status**: Requires Ouroboros repository setup, model files, and Chemical* scripts.

## (Removed) Sequence-based predictor
Removed: no sequence–ligand affinity predictor is included in the default pipeline.

## Pharmacophore
- RDKit FeatureFactory + ETKDG conformers. Computes a simple overlap score based on feature family proximity.
- Requires 3D conformers (enabled by default in `ligand_prep.conformers`).
- Configurable minimum pharmacophore support threshold (default: 0.40 for discovery runs).

## Fusion
- Z-score each component per run.
- Combined score: `S = w_morgan*z(morgan) + w_e3fp*z(e3fp) + w_gm*z(geminimol) + w_qsar*z(qsar) + w_ph4*z(ph4) + w_dock*z(docking) + w_evd*evidence`.
- Weights are config-driven; set component weight to 0 if disabled.
- Default weights: Morgan 2D (30%), E3FP (10%), GeminiMol (25%), pharmacophore (25%), evidence (10%), docking (0%), QSAR (0%).

## Docking (gated and optional)
- Receptor/ligand to PDBQT via Open Babel (`obabel`).
- Docking via AutoDock Vina with a pocket box derived from a validated co-crystal ligand.
- Validation gates: CCD additive filter, SIFTS mapping (PDBe API), optional BioLiP/PDBe ligand selection (auto-validate).
- Docking scores contribute to fusion as `docking_score = -Vina_energy` (more negative energy => higher score).
- For structure preparation and pocket analysis, optionally integrate PDB Prepare Wizard for chain selection, HETATM handling, and PLIP-based pocket characterization.
