"""
Results schema for ranked outputs and provenance.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional
from pydantic import BaseModel, Field


class RankedResult(BaseModel):
    compound_id: str
    compound_smiles: Optional[str] = None
    compound_inchikey: Optional[str] = None
    compound_source: Optional[str] = None
    target_uniprot: str
    target_chembl: Optional[str] = None
    target_gene: Optional[str] = None
    evidence_strength: float
    # Comparator diagnostics
    n_comparators: Optional[int] = None
    median_pchembl: Optional[float] = None
    assay_consistency: Optional[float] = None
    gate_reason: Optional[str] = None
    # Similarity features
    cosine_max: Optional[float] = None
    tanimoto_max: Optional[float] = None
    # E3FP features (nullable)
    e3fp_tanimoto_max: Optional[float] = None
    e3fp_tanimoto_mean_topk: Optional[float] = None
    # GeminiMol features (nullable)
    gm_cosine_max: Optional[float] = None
    gm_cosine_mean_topk: Optional[float] = None
    gm_profile_score: Optional[float] = None
    # Other features
    ph4_best: Optional[float] = None
    qsar_score: Optional[float] = None
    docking_score: Optional[float] = None
    fused_score: float
    rank: int
    provenance: Dict[str, Any] = Field(default_factory=dict)
