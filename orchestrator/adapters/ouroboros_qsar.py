"""
QSAR adapter using Ouroboros (scaffold for integration).

The functions here define a narrow interface expected by the pipeline. They do
not fabricate models. If sufficient labeled data are not available or the
environment lacks dependencies, they raise clear errors so the pipeline can
optionally skip QSAR.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional


@dataclass
class ModelCard:
    model: Any
    metadata: Dict


def fit_qsar_for_target(records: List[Dict]) -> ModelCard:
    """Fit a QSAR model for a target if sufficient labeled data are provided.

    Expect `records` with fields including at least `smiles` and `pchembl_value`.
    If the dataset is too small or environment is missing dependencies, raise
    a RuntimeError; do not fabricate models or labels.
    """
    if len(records) < 20:
        raise RuntimeError("Insufficient labeled data for QSAR (need >=20)")
    raise RuntimeError("QSAR backend not integrated in this environment")


def predict_qsar(model_card: ModelCard, query_smiles: str) -> float:
    """Predict QSAR score for a query SMILES using a fitted model.

    Raises a RuntimeError in this scaffold since no real model is loaded.
    """
    raise RuntimeError("QSAR prediction unavailable (no model loaded)")

