import os
import json
from pathlib import Path

import pytest

from orchestrator.pipeline_nondocking import NonDockingPipeline


@pytest.mark.skipif(not Path("data/compounds").exists(), reason="no real inputs present")
def test_e2e_nondocking_minimal(tmp_path):
    # Use any existing .smi/.sdf files from data/compounds if available
    inputs = []
    for pattern in ("*.smi", "*.sdf"):
        inputs.extend(str(p) for p in Path("data/compounds").glob(pattern))

    if not inputs:
        pytest.skip("no ligand inputs available")

    cfg = {
        "chembl": {"min_pchembl": 6.0, "min_comparators": 5, "require_assay_consistency": True},
        "similarity": {"top_k": 10, "min_tanimoto_morgan": 0.4},
        "pharmacophore": {"method": "rdkit_features"},
        "scoring": {"weights": {"similarity": 0.6, "qsar": 0.0, "pharmacophore": 0.0, "docking": 0.0, "evidence": 0.4}},
        "ligand_prep": {"input_paths": inputs, "ph": 7.4, "conformers": {"enable": False}},
        "output": {"dir": str(tmp_path)},
        "disease_terms": ["example"],
        "organism": "Homo sapiens",
    }

    pipeline = NonDockingPipeline(cfg)
    out = pipeline.run()

    csv_path = Path(out["results_csv"])
    man_path = Path(out["manifest"])
    assert csv_path.exists()
    assert man_path.exists()

    manifest = json.loads(man_path.read_text())
    assert "stages" in manifest
