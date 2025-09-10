"""
RCSB CCD (Chemical Component Dictionary) helpers and additive filtering.

This module provides a deny-list for common crystallization additives and a
resolver stub for CCD entries. Do not fabricate CCD metadata; when offline,
the resolver should raise a clear error.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


_ADDITIVE_DENYLIST = {
    # Solvents/waters
    "HOH", "WAT", "DOD",
    # Common small additives/buffers
    "EDO", "GOL", "PEG", "MPD", "IPA", "TRS", "TRIS", "MES", "HEP", "HEZ", "HEPES",
    "SO4", "PO4", "CL", "NA", "K", "MG", "CA", "ZN", "MN", "CO", "CU",
    # Glycerol/ethylene glycols
    "GOL", "PG4", "PGE", "1PE", "2PE",
    # Detergents / lipids (common)
    "BME", "NAG", "NGZ", "BOG", "CHS",
}


def is_common_additive(het: str) -> bool:
    """Return True if CCD HET code is a common crystallization additive.

    The list is curated from frequently encountered non-biological ligands.
    It is conservative and can be extended via configuration if needed.
    """
    if not het:
        return False
    return het.strip().upper() in _ADDITIVE_DENYLIST


@dataclass
class CCDEntry:
    """Minimal CCD entry representation.

    name: Human-readable name
    type: A short category, e.g., 'ligand', 'ion', 'solvent'
    """

    name: str
    type: str


def resolve_ccd(het: str) -> CCDEntry:
    """Resolve minimal CCD information for a HET code.

    This requires either network access to the RCSB API or a local CCD dump.
    When unavailable, a RuntimeError is raised. Callers should not fabricate
    CCD details; they may fall back to `is_common_additive` for gating.
    """
    if not het:
        raise ValueError("HET code is required")
    raise RuntimeError(
        "CCD resolution requires network or local CCD dataset; unavailable in this run."
    )

