"""
Pharmacophore feature overlap utilities using RDKit FeatureFactory.

Matches features between a query and a comparator set using ETKDG conformers.
When RDKit is unavailable or 3D embedding fails, returns an empty result
without fabricating matches.
"""

from __future__ import annotations

from typing import Dict, List, Tuple

from loguru import logger


def _try_import_rdkit():
    try:
        from rdkit import Chem  # type: ignore
        from rdkit.Chem import AllChem, ChemicalFeatures, rdMolDescriptors  # type: ignore
        from rdkit import RDConfig  # type: ignore
        import os
        factory_path = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
        return Chem, AllChem, ChemicalFeatures, factory_path
    except Exception:
        return None, None, None, None


def _embed_minimize(mol, max_confs=10, prune_rms=0.5):
    Chem, AllChem, _, _ = _try_import_rdkit()
    if Chem is None:
        return False
    try:
        params = AllChem.ETKDGv3()
        params.randomSeed = 0xC0FFEE
        params.pruneRmsThresh = float(prune_rms)
        AllChem.EmbedMultipleConfs(mol, numConfs=int(max_confs), params=params)
        for cid in range(mol.GetNumConformers()):
            try:
                AllChem.MMFFOptimizeMolecule(mol, confId=cid, maxIters=200)
            except Exception:
                pass
        return True
    except Exception as e:
        logger.warning(f"3D embedding failed: {e}")
        return False


def rdkit_feature_ph4(query_smiles: str, comparator_smiles: List[str]) -> Dict:
    """Compute simple pharmacophore feature overlap metrics between query and comparators.

    Returns a dict with: {ph4_best, ph4_mean_topk, n_matches}.
    When RDKit is unavailable, returns zeros and n_matches=0.
    """
    Chem, AllChem, ChemicalFeatures, factory_path = _try_import_rdkit()
    if Chem is None:
        logger.warning("RDKit not available; pharmacophore features disabled")
        return {"ph4_best": 0.0, "ph4_mean_topk": 0.0, "n_matches": 0}

    qmol = Chem.MolFromSmiles(query_smiles)
    if qmol is None:
        return {"ph4_best": 0.0, "ph4_mean_topk": 0.0, "n_matches": 0}

    _embed_minimize(qmol)
    try:
        factory = ChemicalFeatures.BuildFeatureFactory(factory_path)
    except Exception as e:
        logger.warning(f"Feature factory error: {e}")
        return {"ph4_best": 0.0, "ph4_mean_topk": 0.0, "n_matches": 0}

    qfeats = factory.GetFeaturesForMol(qmol)
    qcoords = [(f.GetFamily(), f.GetPos()) for f in qfeats]

    scores: List[float] = []
    for smi in comparator_smiles:
        cmol = Chem.MolFromSmiles(smi)
        if cmol is None:
            continue
        _embed_minimize(cmol)
        cfeats = factory.GetFeaturesForMol(cmol)
        ccoords = [(f.GetFamily(), f.GetPos()) for f in cfeats]
        # Simple overlap by matching families and proximity (<=2.0 Å)
        matches = 0
        for fam1, p1 in qcoords:
            for fam2, p2 in ccoords:
                if fam1 != fam2:
                    continue
                d = (p1 - p2).Length()
                if d <= 2.0:
                    matches += 1
                    break
        if (len(qcoords) + len(ccoords)) > 0:
            score = matches / max(1, min(len(qcoords), len(ccoords)))
            scores.append(float(score))

    if not scores:
        return {"ph4_best": 0.0, "ph4_mean_topk": 0.0, "n_matches": 0}

    scores.sort(reverse=True)
    topk = scores[: min(5, len(scores))]
    return {"ph4_best": float(scores[0]), "ph4_mean_topk": float(sum(topk) / len(topk)), "n_matches": len(scores)}

