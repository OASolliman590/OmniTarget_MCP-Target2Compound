"""
Molecular visualization helpers.

Provides lightweight 2D rendering grid via RDKit. 3D viewer integration
is left as an optional future extension to avoid new dependencies.
"""

from __future__ import annotations

from pathlib import Path
from typing import List


def smiles_grid(smiles_list: List[str], legends: List[str] | None = None, n_cols: int = 5):
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw
    except Exception as e:
        raise RuntimeError("RDKit is required for molecular rendering") from e

    mols = []
    for smi in smiles_list:
        m = Chem.MolFromSmiles(smi)
        if m is not None:
            mols.append(m)

    if not mols:
        raise ValueError("No valid molecules to render")

    legends = legends or ["" for _ in mols]
    img = Draw.MolsToGridImage(mols, molsPerRow=n_cols, legends=legends, subImgSize=(200, 200))
    return img

