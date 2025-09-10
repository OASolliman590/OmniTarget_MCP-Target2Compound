# OmniTarget MCP-Target2Compound (Non-Docking Default)

This pipeline performs ligand-to-target prediction without docking by default. It leverages curated ChEMBL comparators, learned molecular representations, optional QSAR, and pharmacophore features. A validated co-crystal pocket can be used later to enable docking, but it is not required for the primary workflow.

Key stages:
- Target discovery and annotation using MCPs (UniProt/STRING/ProteinAtlas/PDB).
- Comparator construction from ChEMBL with pChEMBL thresholds and assay consistency gating.
- Similarity features from GeminiMol embeddings (when available) and RDKit Morgan Tanimoto as a cross-check.
- Optional QSAR training per target when data suffices.
- RDKit pharmacophore feature overlap against top-k comparators.
- Score fusion on z-scored components with configurable weights; results ranked and written to CSV with a manifest for provenance.

Optional docking (gated):
- If a validated co-crystal pocket exists (BioLiP + CCD additive filter + SIFTS mapping) and tools are available (Open Babel + AutoDock Vina), the pipeline can run docking and produce a separate `docking_results.csv`. Docking is disabled by default and does not affect fused scores unless you later enable it in config and adjust fusion.

Config keys for docking:
- `docking.enabled`: boolean; defaults to false.
- `docking.receptor_pdb_path`: path to receptor PDB/mmCIF.
- `docking.cocrystal.ligand_resname`: ligand residue name to define pocket; skipped if common additive.
- `docking.cocrystal.chain_id`, `docking.cocrystal.uniprot_id`: if provided, SIFTS mapping is enforced.
- `docking.auto_validate` + `docking.pdb_id`: when true and a PDB ID is provided, the pipeline attempts to fetch ligands (PDBe) and select a relevant ligand automatically; falls back gracefully if offline.

Inputs:
- Ligands in .mol/.mol2/.sdf/.smi (mix allowed). The pipeline canonicalizes SMILES for ML; generates 3D SDF for pharmacophore.
- Disease terms (for target discovery) and organism.

Outputs:
- `results.csv` with fused ranking and component features.
- `manifest.json` with provenance: versions, endpoints/queries, comparator counts, thresholds, weights, and file hashes.

See `docs/Methods.md` for implementation details and `configs/run.nd.example.yaml` for a complete configuration.
- PDB Prepare Wizard (optional):
  You may leverage the PDB Prepare Wizard to download/clean structures, extract ligands, and analyze pockets (distance/PLIP) with druggability estimates. This can improve receptor preparation and pocket definition. See: https://github.com/OASolliman590/pdb-prepare-wizard
