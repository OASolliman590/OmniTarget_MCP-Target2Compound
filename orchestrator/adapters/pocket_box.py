"""
Pocket box helpers from validated co-crystal ligands.

Docking remains disabled by default. This adapter exposes a function to define
an axis-aligned box around a specified ligand residue in a structure file
(PDB/mmCIF), using RDKit when available.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple

from ..schemas.pocket import PocketInfo


def _try_import_rdkit():
    try:
        from rdkit import Chem  # type: ignore
        from rdkit.Chem import AllChem  # type: ignore
        return Chem, AllChem
    except Exception as e:
        raise ImportError(
            "RDKit is required for pocket box generation; not available."
        ) from e


def box_from_cocrystal(
    structure_path: Path, ligand_resname: str, margin: float = 12.0
) -> PocketInfo:
    """Compute a docking box around a co-crystal ligand by residue name.

    Parameters
    ----------
    structure_path: Path
        Path to PDB/mmCIF file.
    ligand_resname: str
        Three/four-letter residue name of the ligand to box.
    margin: float
        Padding (Å) added around the ligand's bounding box.

    Returns
    -------
    PocketInfo

    Raises
    ------
    ImportError
        If RDKit is not installed.
    RuntimeError
        If the ligand cannot be located in the structure.
    """
    Chem, _ = _try_import_rdkit()

    if not structure_path.exists():
        raise FileNotFoundError(f"Structure file not found: {structure_path}")

    mol = None
    if structure_path.suffix.lower() == ".pdb":
        mol = Chem.MolFromPDBFile(str(structure_path), sanitize=False, removeHs=False)
    else:
        # Attempt to parse via RDKit for other formats; may require additional
        # integration in future.
        mol = Chem.MolFromPDBFile(str(structure_path), sanitize=False, removeHs=False)

    if mol is None:
        raise RuntimeError("Could not parse structure for pocket definition")

    resn = ligand_resname.strip().upper()
    conf = mol.GetConformer()
    coords = []

    for atom in mol.GetAtoms():
        info = mol.GetAtomWithIdx(atom.GetIdx()).GetPDBResidueInfo()
        if info is None:
            continue
        if (info.GetResidueName() or "").strip().upper() == resn:
            pos = conf.GetAtomPosition(atom.GetIdx())
            coords.append((pos.x, pos.y, pos.z))

    if not coords:
        raise RuntimeError(f"Ligand residue {resn} not found in structure")

    xs, ys, zs = zip(*coords)
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    min_z, max_z = min(zs), max(zs)

    cx = (min_x + max_x) / 2.0
    cy = (min_y + max_y) / 2.0
    cz = (min_z + max_z) / 2.0

    size_x = (max_x - min_x) + 2 * margin
    size_y = (max_y - min_y) + 2 * margin
    size_z = (max_z - min_z) + 2 * margin

    return PocketInfo(
        source="cocrystal",
        pdb_id="",  # populate upstream
        chain_id=None,  # populate upstream
        center=(cx, cy, cz),
        size=(size_x, size_y, size_z),
        druggability=None,
        provenance={"biolip": True, "ligand_id": resn, "sifts_ok": None},
    )

