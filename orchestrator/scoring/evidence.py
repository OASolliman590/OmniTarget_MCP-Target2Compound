"""
Evidence strength calculations for comparator sets.
"""

from __future__ import annotations

from typing import Dict, List


def evidence_strength(records: List[dict]) -> float:
    """Combine comparator characteristics into a strength score in [0,1].

    Inputs are records per-target for actives (e.g., from ChEMBL). The exact
    functional form can be simple and monotonic without fabricating values.

    Heuristic:
    - Size: saturating to 1.0 around N=50
    - Median pChEMBL: scaled from 5.0–8.0 to 0–1
    - Assay consistency: already 0–1
    Final strength is the mean of the three components.
    """
    if not records:
        return 0.0

    # Expect one summary record or a dict
    rec = records[0] if isinstance(records[0], dict) else {}
    size = float(rec.get("size", 0.0))
    median_p = float(rec.get("median_pchembl", 0.0))
    assay_c = float(rec.get("assay_consistency_score", 0.0))

    size_term = min(size / 50.0, 1.0)
    p_term = max(0.0, min((median_p - 5.0) / 3.0, 1.0))
    c_term = max(0.0, min(assay_c, 1.0))

    return float((size_term + p_term + c_term) / 3.0)


def similarity_features(query: str, neighbors: List[dict]) -> dict:
    """Aggregate similarity features from neighbor results.

    neighbors should include keys: 'cosine' (optional), 'tanimoto' (optional).
    Returns a dict with cosine_max, cosine_mean_top5, tanimoto_max, n_cos_ge_tau.
    """
    cos = [n.get("cosine") for n in neighbors if n.get("cosine") is not None]
    tan = [n.get("tanimoto") for n in neighbors if n.get("tanimoto") is not None]

    cos_sorted = sorted(cos, reverse=True) if cos else []
    tan_sorted = sorted(tan, reverse=True) if tan else []

    top5 = cos_sorted[:5] if cos_sorted else []
    mean_top5 = sum(top5) / len(top5) if top5 else 0.0
    cos_max = cos_sorted[0] if cos_sorted else 0.0
    tan_max = tan_sorted[0] if tan_sorted else 0.0
    n_cos_ge_tau = sum(1 for v in cos_sorted if v is not None and v >= 0.6)

    return {
        "cosine_max": float(cos_max),
        "cosine_mean_top5": float(mean_top5),
        "tanimoto_max": float(tan_max),
        "n_cos_ge_tau": int(n_cos_ge_tau),
    }

