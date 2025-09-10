"""
SIFTS mapping utilities.

These functions check mappings between PDB chains and UniProt accessions, and
whether a residue contact region is covered by the UniProt-PDB mapping.

This module does not fabricate mappings. When mapping data are unavailable
(e.g., offline), functions return False or raise a clear error depending on
context. Docking is disabled by default; callers should only use these
functions when a co-crystal is present and validation is explicitly enabled.
"""

from __future__ import annotations

from typing import List

import httpx


def maps_to_uniprot_chain(pdb_id: str, chain: str, uniprot: str) -> bool:
    """Return True if PDB chain maps to the given UniProt accession (via PDBe/SIFTS).

    Uses the PDBe mappings API when network is available. If unavailable or any
    error occurs, returns False (conservative).
    """
    if not (pdb_id and chain and uniprot):
        return False
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"
    try:
        resp = httpx.get(url, timeout=10.0)
        if resp.status_code != 200:
            return False
        data = resp.json().get(pdb_id.lower(), {})
        # Data structure: { "UniProt": {"<uniprot>": [{"chain_id": "A", ...}, ...]}}
        up = data.get("UniProt", {})
        entries = up.get(uniprot, [])
        for e in entries:
            if e.get("chain_id", "").upper() == chain.upper():
                return True
        return False
    except Exception:
        return False


def contact_region_is_mapped(pdb_id: str, chain: str, uniprot: str, residues: List[int]) -> bool:
    """Return True if all residues map to UniProt via SIFTS (PDBe API).

    Conservative: any failure returns False.
    """
    if not residues:
        return False
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"
    try:
        resp = httpx.get(url, timeout=10.0)
        if resp.status_code != 200:
            return False
        data = resp.json().get(pdb_id.lower(), {})
        up = data.get("UniProt", {})
        entries = up.get(uniprot, [])
        # Each entry may contain residue ranges under 'mappings'
        ranges = []
        for e in entries:
            if e.get("chain_id", "").upper() == chain.upper():
                for m in e.get("mappings", []) or []:
                    start = m.get("start", {}).get("author_residue_number") or m.get("start", {}).get("residue_number")
                    end = m.get("end", {}).get("author_residue_number") or m.get("end", {}).get("residue_number")
                    if start is not None and end is not None:
                        ranges.append((int(start), int(end)))
        if not ranges:
            return False
        for r in residues:
            if not any(start <= r <= end for start, end in ranges):
                return False
        return True
    except Exception:
        return False
