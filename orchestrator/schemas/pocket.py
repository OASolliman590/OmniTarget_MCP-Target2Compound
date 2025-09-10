"""
Pocket and binding box schemas.
"""

from __future__ import annotations

from typing import Literal, Optional, Tuple, Dict
from pydantic import BaseModel, Field


class PocketInfo(BaseModel):
    """Validated pocket information used for docking box definition.

    source: 'cocrystal' when derived from a validated co-crystal ligand; 'finder'
            reserved for future pocket detection tools; 'none' when disabled.
    provenance: capture flags and metrics for traceability.
    """

    source: Literal["cocrystal", "finder", "none"] = Field(default="none")
    pdb_id: str = Field(default="")
    chain_id: Optional[str] = None
    center: Optional[Tuple[float, float, float]] = None
    size: Optional[Tuple[float, float, float]] = None
    druggability: Optional[float] = None
    provenance: Dict = Field(default_factory=dict)

