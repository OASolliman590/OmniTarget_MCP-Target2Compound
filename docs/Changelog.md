# Changelog

## 0.2.1 (Remove sequence-based predictor)
- Removed the previous sequence-based affinity predictor and its configuration; default pipeline remains non-docking (similarity + pharmacophore + evidence). Optional docking stays gated and off by default.

## 0.2.0 (Non-Docking Default)
- Added non-docking pipeline (`run_nd`) with ChEMBL comparators, GeminiMol/Ouroboros stubs, and RDKit pharmacophore.
- Implemented ligand IO supporting .mol/.mol2/.sdf/.smi with 3D conformers.
- Added provenance manifest and docs scaffolding.
