"""
Results schema for ranked outputs and provenance.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional
from pydantic import BaseModel, Field


class RankedResult(BaseModel):
    compound_id: str
    target_uniprot: str
    target_chembl: Optional[str] = None
    evidence_strength: float
    cosine_max: Optional[float] = None
    tanimoto_max: Optional[float] = None
    ph4_best: Optional[float] = None
    qsar_score: Optional[float] = None
    fused_score: float
    rank: int
    provenance: Dict[str, Any] = Field(default_factory=dict)
