"""
E3FP 3D molecular fingerprint adapter.

E3FP (Extended 3-Dimensional FingerPrint) is a 3D, RDKit-integrated fingerprint 
method inspired by ECFP. This adapter provides 3D similarity analysis for 
ligand-based target prediction.

References:
- GitHub: https://github.com/keiserlab/e3fp
- PyPI: https://pypi.org/project/e3fp/
- Docs: https://e3fp.readthedocs.io/
"""

from __future__ import annotations

import numpy as np
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from loguru import logger


def _try_import_e3fp():
    """Try to import E3FP and RDKit dependencies."""
    try:
        import e3fp
        from e3fp.conformer.util import mol_from_sdf
        from e3fp.fingerprint.fprint import Fingerprint
        from rdkit import Chem
        from rdkit.Chem import AllChem
        return e3fp, mol_from_sdf, Fingerprint, Chem, AllChem
    except ImportError as e:
        logger.debug(f"E3FP not available: {e}")
        return None, None, None, None, None


def _try_import_rdkit():
    """Try to import RDKit for fallback operations."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, DataStructs
        return Chem, AllChem, DataStructs
    except ImportError:
        return None, None, None


def compute_e3fp_for_sdf(
    sdf_path: Path,
    radius: int = 2,
    shells: int = 5,
    max_confs: int = 10,
    prune_rms: float = 0.5
) -> Dict[str, np.ndarray]:
    """
    Build (or read) 3D conformers from the SDF, then compute E3FP bitvectors.
    
    Parameters
    ----------
    sdf_path : Path
        Path to SDF file containing 3D conformers
    radius : int, default=2
        E3FP radius parameter
    shells : int, default=5
        E3FP shells parameter
    max_confs : int, default=10
        Maximum number of conformers to process
    prune_rms : float, default=0.5
        RMS threshold for conformer pruning
        
    Returns
    -------
    Dict[str, np.ndarray]
        Dictionary mapping ligand names to E3FP bitvectors.
        Returns empty dict if E3FP/RDKit not available (no fabrication).
    """
    e3fp, mol_from_sdf, Fingerprint, Chem, AllChem = _try_import_e3fp()
    
    if e3fp is None:
        logger.warning("E3FP not available; cannot compute 3D fingerprints")
        return {}
    
    if not sdf_path.exists():
        logger.warning(f"SDF file not found: {sdf_path}")
        return {}
    
    try:
        # Load molecules from SDF
        mols = mol_from_sdf(str(sdf_path))
        if not mols:
            logger.warning(f"No molecules found in SDF: {sdf_path}")
            return {}
        
        results = {}
        for i, mol in enumerate(mols[:max_confs]):
            if mol is None:
                continue
                
            # Get molecule name or use index
            name = mol.GetProp('_Name') if mol.HasProp('_Name') else f"mol_{i}"
            
            try:
                # Generate RDKit fingerprint first
                from rdkit.Chem import rdMolDescriptors
                rdkit_fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=1024)
                
                # Generate E3FP fingerprint using the correct API
                fp = Fingerprint.from_rdkit(
                    rdkit_fp, 
                    radius=radius, 
                    first=0,
                    last=shells
                )
                
                # Convert to numpy array
                bitvector = fp.to_vector()
                results[name] = bitvector
                
            except Exception as e:
                logger.warning(f"Failed to compute E3FP for {name}: {e}")
                continue
        
        logger.info(f"Computed E3FP fingerprints for {len(results)} molecules from {sdf_path}")
        return results
        
    except Exception as e:
        logger.error(f"Error processing SDF file {sdf_path}: {e}")
        return {}


def tanimoto_e3fp(query_fp: np.ndarray, ref_fp: np.ndarray) -> float:
    """
    Compute Tanimoto similarity between two E3FP bitvectors.
    
    Parameters
    ----------
    query_fp : np.ndarray
        Query E3FP fingerprint (boolean vector or sparse matrix)
    ref_fp : np.ndarray
        Reference E3FP fingerprint (boolean vector or sparse matrix)
        
    Returns
    -------
    float
        Tanimoto similarity coefficient [0, 1]
    """
    if query_fp.shape != ref_fp.shape:
        logger.warning("E3FP fingerprints have different shapes")
        return 0.0
    
    # Handle sparse matrices
    from scipy.sparse import issparse
    if issparse(query_fp) or issparse(ref_fp):
        # Convert sparse matrices to dense arrays
        query_dense = query_fp.toarray().flatten() if issparse(query_fp) else query_fp.flatten()
        ref_dense = ref_fp.toarray().flatten() if issparse(ref_fp) else ref_fp.flatten()
    else:
        query_dense = query_fp.flatten()
        ref_dense = ref_fp.flatten()
    
    # Convert to boolean if needed
    query_bits = query_dense.astype(bool)
    ref_bits = ref_dense.astype(bool)
    
    # Compute Tanimoto coefficient
    intersection = np.sum(query_bits & ref_bits)
    union = np.sum(query_bits | ref_bits)
    
    if union == 0:
        return 0.0
    
    return float(intersection / union)


def topk_similarity_matrix_e3fp(
    query_fps: Dict[str, np.ndarray],
    ref_fps: Dict[str, np.ndarray],
    k: int = 25
) -> Dict[str, Dict[str, Any]]:
    """
    Compute top-k E3FP similarities between query and reference fingerprints.
    
    Parameters
    ----------
    query_fps : Dict[str, np.ndarray]
        Query E3FP fingerprints {name: bitvector}
    ref_fps : Dict[str, np.ndarray]
        Reference E3FP fingerprints {name: bitvector}
    k : int, default=25
        Number of top similarities to return
        
    Returns
    -------
    Dict[str, Dict[str, Any]]
        For each query: {
            'e3fp_tanimoto_max': float,
            'e3fp_tanimoto_mean_topk': float,
            'neighbors': [(ref_id, score), ...]
        }
    """
    if not query_fps or not ref_fps:
        logger.warning("Empty fingerprint dictionaries provided")
        return {}
    
    results = {}
    
    for query_name, query_fp in query_fps.items():
        similarities = []
        
        for ref_name, ref_fp in ref_fps.items():
            try:
                sim = tanimoto_e3fp(query_fp, ref_fp)
                similarities.append((ref_name, sim))
            except Exception as e:
                logger.warning(f"Error computing E3FP similarity {query_name} vs {ref_name}: {e}")
                continue
        
        # Sort by similarity (descending)
        similarities.sort(key=lambda x: x[1], reverse=True)
        
        # Get top-k
        top_k = similarities[:k]
        
        if top_k:
            max_sim = top_k[0][1]
            mean_topk = sum(sim for _, sim in top_k) / len(top_k)
        else:
            max_sim = 0.0
            mean_topk = 0.0
        
        results[query_name] = {
            'e3fp_tanimoto_max': max_sim,
            'e3fp_tanimoto_mean_topk': mean_topk,
            'neighbors': top_k
        }
    
    logger.info(f"Computed E3FP similarities for {len(results)} query molecules")
    return results


def compute_e3fp_features_for_ligands(
    ligand_smiles: List[str],
    ligand_names: List[str],
    radius: int = 2,
    shells: int = 5
) -> Dict[str, np.ndarray]:
    """
    Compute E3FP fingerprints directly from SMILES strings.
    
    This is a fallback when SDF files are not available.
    
    Parameters
    ----------
    ligand_smiles : List[str]
        List of SMILES strings
    ligand_names : List[str]
        List of ligand names
    radius : int, default=2
        E3FP radius parameter
    shells : int, default=5
        E3FP shells parameter
        
    Returns
    -------
    Dict[str, np.ndarray]
        Dictionary mapping ligand names to E3FP bitvectors
    """
    e3fp, _, Fingerprint, Chem, AllChem = _try_import_e3fp()
    
    if e3fp is None:
        logger.warning("E3FP not available; cannot compute 3D fingerprints from SMILES")
        return {}
    
    results = {}
    
    for smiles, name in zip(ligand_smiles, ligand_names):
        try:
            # Convert SMILES to molecule
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                logger.warning(f"Invalid SMILES for {name}: {smiles}")
                continue

            # Generate 3D conformer with ETKDG, robust fallback
            mol = Chem.AddHs(mol)
            try:
                if hasattr(AllChem, "ETKDGv3"):
                    params = AllChem.ETKDGv3()
                else:
                    params = AllChem.ETKDG()
                params.randomSeed = 42
                params.pruneRmsThresh = 0.5
                conf_id = AllChem.EmbedMolecule(mol, params)
            except Exception:
                conf_id = AllChem.EmbedMolecule(mol, randomSeed=42)

            if conf_id < 0:
                logger.warning(f"Failed to embed 3D conformer for {name}")
                continue

            # Optimize geometry (MMFF then UFF fallback)
            try:
                AllChem.MMFFOptimizeMolecule(mol)
            except Exception:
                try:
                    AllChem.UFFOptimizeMolecule(mol)
                except Exception:
                    logger.warning(f"Forcefield optimization failed for {name}")

            # Generate RDKit fingerprint first
            from rdkit.Chem import rdMolDescriptors
            rdkit_fp = rdMolDescriptors.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=1024)

            # Generate E3FP fingerprint
            fp = Fingerprint.from_rdkit(
                rdkit_fp,
                radius=radius,
                first=0,
                last=shells
            )

            bitvector = fp.to_vector()
            results[name] = bitvector

        except Exception as e:
            logger.warning(f"Failed to compute E3FP for {name}: {e}")
            continue
    
    logger.info(f"Computed E3FP fingerprints for {len(results)} molecules from SMILES")
    return results
