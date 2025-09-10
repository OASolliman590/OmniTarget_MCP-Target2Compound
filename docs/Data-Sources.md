# Data Sources

- KEGG/Reactome: pathway associations used during target discovery.
- STRING: protein–protein interactions supporting target context.
- ProteinAtlas: tissue expression filters and annotations.
- UniProt: canonical target identifiers and sequences.
- PDB: structure availability and co-crystal context.
- PDBe/SIFTS: PDB–UniProt mappings and ligand annotations for pocket validation.
- ChEMBL: active ligand comparators per target (pChEMBL, assays).

Each data source is accessed via an MCP service in production. The pipeline logs endpoints, parameters, and counts to `manifest.json`.
