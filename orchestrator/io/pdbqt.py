"""
Utilities to convert ligands and receptors to PDBQT using Open Babel.

No placeholders: if Open Babel or inputs are unavailable, raise a clear error.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional


def _which(cmd: str) -> Optional[str]:
    return shutil.which(cmd)


def require_obabel() -> str:
    ob = _which("obabel") or _which("obabel.exe")
    if not ob:
        raise RuntimeError("Open Babel 'obabel' not found in PATH; required for PDBQT conversion")
    return ob


def ligand_smiles_to_pdbqt(smiles: str, out_path: Path) -> Path:
    """Convert a SMILES to PDBQT via RDKit SDF + Open Babel PDBQT.

    Raises RuntimeError if conversion tools are not available or fail.
    """
    try:
        from rdkit import Chem  # type: ignore
        from rdkit.Chem import AllChem  # type: ignore
    except Exception as e:
        raise RuntimeError("RDKit not available for ligand 3D generation") from e

    obabel = require_obabel()

    # Build RDKit 3D and write temp SDF
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise RuntimeError("Invalid SMILES for PDBQT conversion")
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol) != 0:
        raise RuntimeError("Failed to embed 3D coordinates for ligand")
    AllChem.MMFFOptimizeMolecule(mol)

    with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as f:
        sdf_path = Path(f.name)
    w = Chem.SDWriter(str(sdf_path))
    w.write(mol)
    w.close()

    # Convert to PDBQT with Gasteiger charges
    out_path = out_path.with_suffix(".pdbqt")
    cmd = [obabel, "-isdf", str(sdf_path), "-opdbqt", "-O", str(out_path), "--partialcharge", "gasteiger"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0 or not out_path.exists():
        raise RuntimeError(f"obabel failed for ligand: {res.stderr.strip()}")
    return out_path


def receptor_to_pdbqt(pdb_path: Path, out_path: Path) -> Path:
    """Convert receptor PDB/mmCIF to PDBQT using Open Babel.

    Raises RuntimeError if conversion tools are not available or fail.
    """
    obabel = require_obabel()
    if not pdb_path.exists():
        raise FileNotFoundError(f"Receptor file not found: {pdb_path}")
    out_path = out_path.with_suffix(".pdbqt")

    # obabel auto-detects by extension; ensure hydrogens and charges
    cmd = [obabel, str(pdb_path), "-O", str(out_path), "--partialcharge", "gasteiger"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0 or not out_path.exists():
        raise RuntimeError(f"obabel failed for receptor: {res.stderr.strip()}")
    return out_path

