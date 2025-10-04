"""
Non-docking score fusion (S) with z-scored components.

S = w_sim * z(sim_feature)
  + w_qsar * z(qsar_pred)      # optional
  + w_ph4 * z(ph4_best)
  + w_dock * z(docking_score)  # optional
  + w_evd * evidence_strength

Missing components are ignored in z-scores and treated as zeros in fusion.
"""

from __future__ import annotations

from typing import Dict, List
import math


def _z(vals: List[float]) -> List[float]:
    xs = [v for v in vals if v is not None]
    if not xs:
        return [0.0 for _ in vals]
    mu = sum(xs) / len(xs)
    var = sum((v - mu) ** 2 for v in xs) / max(1, len(xs) - 1)
    sd = math.sqrt(var) if var > 0 else 1.0
    return [float((0.0 if v is None else (v - mu) / sd)) for v in vals]


def fuse(records: List[Dict], weights: Dict[str, float]) -> List[Dict]:
    """Fuse non-docking features per spec into a final score.

    Expected record keys:
    - similarity: use cosine_max or cosine_mean_top5 (if present)
    - e3fp_tanimoto_max: optional E3FP similarity
    - qsar_score: optional float
    - ph4_best: optional float
    - evidence_strength: float in [0,1]
    """
    sim_vals = []
    e3fp_vals = []
    gm_vals = []
    qsar_vals = []
    ph4_vals = []
    evd_vals = []
    dock_vals = []

    for r in records:
        sim = r.get("cosine_max") or r.get("cosine_mean_top5")
        sim_vals.append(sim if sim is not None else None)
        e3fp_vals.append(r.get("e3fp_tanimoto_max"))
        gm_vals.append(r.get("gm_cosine_max"))
        qsar_vals.append(r.get("qsar_score"))
        ph4_vals.append(r.get("ph4_best"))
        evd_vals.append(r.get("evidence_strength", 0.0))
        dock_vals.append(r.get("docking_score"))

    sim_z = _z([0.0 if v is None else float(v) for v in sim_vals])
    e3fp_z = _z([0.0 if v is None else float(v) for v in e3fp_vals])
    gm_z = _z([0.0 if v is None else float(v) for v in gm_vals])
    qsar_z = _z([0.0 if v is None else float(v) for v in qsar_vals])
    ph4_z = _z([0.0 if v is None else float(v) for v in ph4_vals])
    dock_z = _z([0.0 if v is None else float(v) for v in dock_vals])

    ws = {
        "morgan_2d": float(weights.get("morgan_2d", 0.30)),
        "e3fp": float(weights.get("e3fp", 0.10)),
        "geminimol": float(weights.get("geminimol", 0.25)),
        "qsar": float(weights.get("qsar", 0.0)),
        "pharmacophore": float(weights.get("pharmacophore", 0.25)),
        "evidence": float(weights.get("evidence", 0.10)),
        "docking": float(weights.get("docking", 0.0)),
    }

    fused: List[Dict] = []
    for i, r in enumerate(records):
        score = (
            ws["morgan_2d"] * sim_z[i]
            + ws["e3fp"] * e3fp_z[i]
            + ws["geminimol"] * gm_z[i]
            + ws["qsar"] * qsar_z[i]
            + ws["pharmacophore"] * ph4_z[i]
            + ws["evidence"] * float(evd_vals[i] or 0.0)
            + ws["docking"] * dock_z[i]
        )
        out = dict(r)
        out.update({
            "fused_score": float(score),
            "morgan_2d_z": float(sim_z[i]),
            "e3fp_z": float(e3fp_z[i]),
            "geminimol_z": float(gm_z[i]),
            "qsar_z": float(qsar_z[i]),
            "ph4_z": float(ph4_z[i]),
            "docking_z": float(dock_z[i]),
        })
        fused.append(out)

    fused.sort(key=lambda x: x.get("fused_score", 0.0), reverse=True)
    for rank, r in enumerate(fused, start=1):
        r["rank"] = rank
    return fused
