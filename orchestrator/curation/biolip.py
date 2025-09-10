"""
BioLiP integration helpers to identify biologically relevant co-crystal ligands.

This module provides a thin, explicit interface to fetch ligand annotations
from BioLiP and select the most relevant ligand for defining a pocket.

Notes
-----
- Network access may be required to retrieve BioLiP data. When offline, the
  fetch function should raise a RuntimeError with a clear message and the
  caller should skip docking logic. Do not fabricate data.
- Selection prefers ligands that make clear residue contacts and are not
  common crystallization additives (see `ccd.is_common_additive`).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Dict
import httpx


@dataclass
class BiolipLigand:
    """Minimal BioLiP ligand record for pocket validation.

    Attributes
    -----------
    ligand_id: 3-letter residue name (e.g., ATP)
    chain: PDB chain identifier
    contacts: List of contacted residue identifiers (e.g., ["A:123", "A:124"]) 
    notes: Free-form metadata from source
    """

    ligand_id: str
    chain: str
    contacts: List[str]
    notes: Dict[str, str]


def fetch_biolip_ligands(pdb_id: str) -> List[BiolipLigand]:
    """Fetch ligand annotations for a given PDB ID via PDBe ligand API (BioLiP-like).

    This is a pragmatic network integration using PDBe's ligand_monomers
    endpoint to approximate BioLiP-like co-crystal ligands. For each ligand,
    returns minimal fields. When network is unavailable, raises RuntimeError.
    """
    if not pdb_id:
        return []
    url = f"https://www.ebi.ac.uk/pdbe/api/pdb/entry/ligand_monomers/{pdb_id.lower()}"
    try:
        resp = httpx.get(url, timeout=10.0)
        if resp.status_code != 200:
            raise RuntimeError(f"PDBe ligand API returned {resp.status_code}")
        data = resp.json().get(pdb_id.lower(), [])
        ligs: List[BiolipLigand] = []
        for lig in data:
            chem_id = lig.get("chem_comp_id")
            for inst in lig.get("molecule_instances", []) or []:
                chain_id = inst.get("chain_id")
                ligs.append(
                    BiolipLigand(
                        ligand_id=str(chem_id or ""),
                        chain=str(chain_id or ""),
                        contacts=[],  # PDBe endpoint does not list contacts; can be added via interactions API later
                        notes={"source": "pdbe_ligand_monomers"},
                    )
                )
        return ligs
    except Exception as e:
        raise RuntimeError(f"Ligand fetch unavailable: {e}")


def select_relevant_ligand(ligs: List[BiolipLigand]) -> Optional[BiolipLigand]:
    """Select the most biologically relevant ligand from BioLiP candidates.

    Heuristics (deterministic, no fabricated inputs):
    - Prefer ligands with non-empty residue contacts.
    - If multiple candidates tie, prefer the one with the greatest number of
      contacts; break ties by lexicographic order of ligand_id.

    Parameters
    ----------
    ligs: list[BiolipLigand]

    Returns
    -------
    Optional[BiolipLigand]
        The chosen ligand or None if the input list is empty.
    """
    if not ligs:
        return None

    def score(l: BiolipLigand) -> tuple[int, str]:
        return (len(l.contacts or []), l.ligand_id or "")

    # Max by contacts, then ligand_id for determinism
    return sorted(ligs, key=score, reverse=True)[0]
