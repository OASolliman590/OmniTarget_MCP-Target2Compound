"""
GeminiMol embedding and nearest neighbor utilities.

When GeminiMol is unavailable, this module does not fabricate embeddings.
It can optionally compute RDKit Morgan fingerprints and Tanimoto similarity
as a cross-check or fallback for similarity features.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, List, Optional, Tuple

import numpy as np
from loguru import logger


def _try_import_rdkit():
    try:
        from rdkit import Chem  # type: ignore
        from rdkit.Chem import AllChem, DataStructs  # type: ignore
        return Chem, AllChem, DataStructs
    except Exception:
        return None, None, None


def embed_smiles(smiles_list: List[str]) -> np.ndarray:
    """Compute GeminiMol embeddings for SMILES.

    Returns
    -------
    np.ndarray
        Array of shape (n, d). If the GeminiMol model is not available in the
        environment, raises a RuntimeError. Do not fabricate embeddings.
    """
    # Placeholder: no model integration here to avoid fabricating data.
    raise RuntimeError("GeminiMol embeddings unavailable (model not loaded in this environment)")


def build_index(embeddings: np.ndarray) -> Any:
    """Build an ANN index for embeddings (e.g., FAISS/Annoy).

    This function is a stub to clarify the expected interface. Implement with
    FAISS/Annoy in production. Offline, raise to signal unavailability.
    """
    raise RuntimeError("ANN index build unavailable (no backend in this environment)")


@dataclass
class Neighbor:
    idx: int
    similarity: float
    tanimoto: Optional[float] = None


def _morgan_fp(smiles: str, radius: int = 2, nbits: int = 2048):
    Chem, AllChem, _ = _try_import_rdkit()
    if Chem is None:
        return None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=nbits)


def nearest_neighbors(
    query_smiles: str,
    index: Any,
    embeddings: np.ndarray,
    k: int,
) -> List[Neighbor]:
    """Find k nearest neighbors using cosine similarity on GeminiMol embeddings.

    As a cross-check, compute RDKit Morgan Tanimoto between query and each
    neighbor (if RDKit is available). If the ANN index is unavailable, raise.
    """
    raise RuntimeError("Nearest neighbor search unavailable (no ANN backend in this environment)")


def tanimoto_neighbors_only(
    query_smiles: str, smiles_list: List[str], k: int
) -> List[Neighbor]:
    """RDKit-only fallback to compute Tanimoto neighbors among SMILES.

    This does not fabricate data: if RDKit is not present, returns an empty
    list and logs a warning.
    """
    Chem, AllChem, DataStructs = _try_import_rdkit()
    if Chem is None:
        logger.warning("RDKit not available; cannot compute Tanimoto neighbors")
        return []
    qfp = _morgan_fp(query_smiles)
    if qfp is None:
        return []
    sims: List[Tuple[int, float]] = []
    for idx, smi in enumerate(smiles_list):
        fp = _morgan_fp(smi)
        if fp is None:
            continue
        sim = DataStructs.TanimotoSimilarity(qfp, fp)
        sims.append((idx, float(sim)))
    sims.sort(key=lambda x: x[1], reverse=True)
    top = sims[: max(0, k)]
    return [Neighbor(idx=i, similarity=s, tanimoto=s) for i, s in top]

