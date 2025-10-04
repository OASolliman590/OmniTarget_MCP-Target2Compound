"""
GeminiMol adapter for molecular embeddings and PharmProfiler integration.

GeminiMol provides molecular embeddings and PharmProfiler for ligand-based 
target identification. This adapter integrates with the GeminiMol repository
following the official setup and usage patterns.

References:
- GitHub: https://github.com/Wang-Lin-boop/GeminiMol
- PharmProfiler for target identification and virtual screening
"""

from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
import numpy as np
from loguru import logger


def _check_geminimol_setup(repo_path: str) -> Tuple[bool, str]:
    """
    Check if GeminiMol is properly set up.
    
    Parameters
    ----------
    repo_path : str
        Path to GeminiMol repository
        
    Returns
    -------
    Tuple[bool, str]
        (is_setup, error_message)
    """
    geminimol_path = Path(repo_path)
    
    if not geminimol_path.exists():
        return False, f"GeminiMol repository not found at {repo_path}"
    
    # Check for required files (README + model dir). PharmProfiler.py may live
    # under geminimol/ in upstream; don't hard-require it here.
    if not (geminimol_path / "README.md").exists():
        return False, "Required README.md not found in GeminiMol repository"
    
    # Check for model directory
    model_path = geminimol_path / "models" / "GeminiMol"
    if not model_path.exists():
        return False, f"Model directory not found at {model_path}"
    
    return True, ""


def _check_environment_variables() -> Tuple[bool, List[str]]:
    """
    Check if required GeminiMol environment variables are set.
    
    Returns
    -------
    Tuple[bool, List[str]]
        (all_set, missing_vars)
    """
    required_vars = [
        "GeminiMol",
        "geminimol_app", 
        "geminimol_lib",
        "geminimol_data"
    ]
    
    missing = []
    for var in required_vars:
        if not os.getenv(var):
            missing.append(var)
    
    return len(missing) == 0, missing


def encode_smiles(smiles_list: List[str], repo_path: str, model_path: str) -> Optional[np.ndarray]:
    """
    Generate GeminiMol embeddings for SMILES strings.
    
    Parameters
    ----------
    smiles_list : List[str]
        List of SMILES strings to encode
    repo_path : str
        Path to GeminiMol repository
    model_path : str
        Path to GeminiMol model
        
    Returns
    -------
    Optional[np.ndarray]
        Embeddings array of shape (n_molecules, embedding_dim) or None if failed
    """
    try:
        # Check setup
        is_setup, error_msg = _check_geminimol_setup(repo_path)
        if not is_setup:
            logger.warning(f"GeminiMol setup check failed: {error_msg}")
            return None
        
        # Check environment variables
        env_ok, missing_vars = _check_environment_variables()
        if not env_ok:
            logger.warning(f"Missing GeminiMol environment variables: {missing_vars}")
            return None
        
        # Create temporary input file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.smi', delete=False) as f:
            for smiles in smiles_list:
                f.write(f"{smiles}\n")
            temp_input = f.name
        
        try:
            # Create temporary output file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.npy', delete=False) as f:
                temp_output = f.name
            
            # Run GeminiMol encoding (placeholder - actual command depends on GeminiMol API)
            # This is a placeholder implementation - the actual command would depend on
            # the specific GeminiMol interface
            cmd = [
                "python", 
                str(Path(repo_path) / "encode_smiles.py"),  # Placeholder script
                "--input", temp_input,
                "--output", temp_output,
                "--model", model_path
            ]
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300  # 5 minute timeout
            )
            
            if result.returncode != 0:
                logger.warning(f"GeminiMol encoding failed: {result.stderr}")
                return None
            
            # Load embeddings
            if Path(temp_output).exists():
                embeddings = np.load(temp_output)
                logger.info(f"Generated GeminiMol embeddings for {len(smiles_list)} molecules")
                return embeddings
            else:
                logger.warning("GeminiMol output file not found")
                return None
                
        finally:
            # Clean up temporary files
            for temp_file in [temp_input, temp_output]:
                try:
                    Path(temp_file).unlink()
                except Exception:
                    pass
                    
    except Exception as e:
        logger.warning(f"GeminiMol encoding failed: {e}")
        return None


def pharm_profile_score(
    query_smiles: str,
    ref_smiles: List[str], 
    ref_labels: List[str],
    repo_path: str,
    model_path: str,
    output_dir: Path
) -> Optional[Dict[str, Any]]:
    """
    Run PharmProfiler for target identification.
    
    Parameters
    ----------
    query_smiles : str
        Query SMILES string
    ref_smiles : List[str]
        Reference SMILES strings
    ref_labels : List[str]
        Reference target labels
    repo_path : str
        Path to GeminiMol repository
    model_path : str
        Path to GeminiMol model
    output_dir : Path
        Directory to save PharmProfiler output
        
    Returns
    -------
    Optional[Dict[str, Any]]
        PharmProfiler results or None if failed
    """
    try:
        # Check setup
        is_setup, error_msg = _check_geminimol_setup(repo_path)
        if not is_setup:
            logger.warning(f"GeminiMol setup check failed: {error_msg}")
            return None
        
        # Check environment variables
        env_ok, missing_vars = _check_environment_variables()
        if not env_ok:
            logger.warning(f"Missing GeminiMol environment variables: {missing_vars}")
            return None
        
        # Create reference data file
        ref_file = output_dir / "pharm_profiler_ref.csv"
        with open(ref_file, 'w') as f:
            f.write("smiles,target\n")
            for smiles, label in zip(ref_smiles, ref_labels):
                f.write(f"{smiles},{label}\n")
        
        # Create query file
        query_file = output_dir / "pharm_profiler_query.smi"
        with open(query_file, 'w') as f:
            f.write(f"{query_smiles}\n")
        
        # Run PharmProfiler
        # Resolve PharmProfiler.py path across layouts or env var
        env_app = os.getenv("geminimol_app")
        candidates = []
        if env_app:
            candidates.append(Path(env_app) / "PharmProfiler.py")
        candidates.extend([
            Path(repo_path) / "geminimol" / "PharmProfiler.py",
            Path(repo_path) / "PharmProfiler.py",
        ])
        pharm_profiler_script = next((p for p in candidates if p and p.exists()), None)
        if pharm_profiler_script is None:
            logger.warning("PharmProfiler.py not found in expected locations; skipping PharmProfiler")
            return None
        
        output_file = output_dir / "pharm_profiler_results.tsv"
        
        cmd = [
            "python",
            str(pharm_profiler_script),
            "--query", str(query_file),
            "--reference", str(ref_file),
            "--output", str(output_file),
            "--model", model_path
        ]
        
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=600  # 10 minute timeout
        )
        
        if result.returncode != 0:
            logger.warning(f"PharmProfiler failed: {result.stderr}")
            return None
        
        # Parse results
        if output_file.exists():
            results = _parse_pharm_profiler_output(output_file)
            logger.info(f"PharmProfiler completed for query: {query_smiles}")
            return results
        else:
            logger.warning("PharmProfiler output file not found")
            return None
            
    except Exception as e:
        logger.warning(f"PharmProfiler failed: {e}")
        return None


def _parse_pharm_profiler_output(output_file: Path) -> Dict[str, Any]:
    """
    Parse PharmProfiler TSV output file.
    
    Parameters
    ----------
    output_file : Path
        Path to PharmProfiler output file
        
    Returns
    -------
    Dict[str, Any]
        Parsed results
    """
    try:
        import pandas as pd
        
        # Read TSV file
        df = pd.read_csv(output_file, sep='\t')
        
        # Extract key metrics
        results = {
            "n_targets": len(df),
            "max_score": float(df.iloc[:, -1].max()) if len(df.columns) > 1 else 0.0,
            "mean_score": float(df.iloc[:, -1].mean()) if len(df.columns) > 1 else 0.0,
            "top_target": df.iloc[df.iloc[:, -1].idxmax(), 0] if len(df.columns) > 1 else None,
            "scores": df.iloc[:, -1].tolist() if len(df.columns) > 1 else []
        }
        
        return results
        
    except Exception as e:
        logger.warning(f"Failed to parse PharmProfiler output: {e}")
        return {"error": str(e)}


def compute_geminimol_features(
    query_smiles: List[str],
    ref_smiles: List[str],
    ref_labels: List[str],
    repo_path: str,
    model_path: str,
    output_dir: Path,
    use_pharm_profiler: bool = True
) -> Dict[str, Dict[str, Any]]:
    """
    Compute GeminiMol features for query molecules.
    
    Parameters
    ----------
    query_smiles : List[str]
        Query SMILES strings
    ref_smiles : List[str]
        Reference SMILES strings
    ref_labels : List[str]
        Reference target labels
    repo_path : str
        Path to GeminiMol repository
    model_path : str
        Path to GeminiMol model
    output_dir : Path
        Output directory for results
    use_pharm_profiler : bool
        Whether to use PharmProfiler
        
    Returns
    -------
    Dict[str, Dict[str, Any]]
        Features for each query molecule
    """
    results = {}
    
    # Generate embeddings
    logger.info("Computing GeminiMol embeddings...")
    query_embeddings = encode_smiles(query_smiles, repo_path, model_path)
    ref_embeddings = encode_smiles(ref_smiles, repo_path, model_path)
    
    if query_embeddings is None or ref_embeddings is None:
        logger.warning("Failed to generate GeminiMol embeddings")
        # Return empty results with None values
        for smiles in query_smiles:
            results[smiles] = {
                "gm_cosine_max": None,
                "gm_cosine_mean_topk": None,
                "gm_profile_score": None
            }
        return results
    
    # Compute cosine similarities
    from sklearn.metrics.pairwise import cosine_similarity
    
    similarities = cosine_similarity(query_embeddings, ref_embeddings)
    
    for i, smiles in enumerate(query_smiles):
        sim_scores = similarities[i]
        
        # Get top similarities
        top_indices = np.argsort(sim_scores)[::-1][:10]  # Top 10
        top_scores = sim_scores[top_indices]
        
        features = {
            "gm_cosine_max": float(np.max(sim_scores)),
            "gm_cosine_mean_topk": float(np.mean(top_scores)),
            "gm_profile_score": None
        }
        
        # Run PharmProfiler if enabled
        if use_pharm_profiler:
            profile_result = pharm_profile_score(
                smiles, ref_smiles, ref_labels, repo_path, model_path, output_dir
            )
            if profile_result and "max_score" in profile_result:
                features["gm_profile_score"] = profile_result["max_score"]
        
        results[smiles] = features
    
    return results


def cosine_similarity_geminimol(
    query_embeddings: np.ndarray,
    ref_embeddings: np.ndarray
) -> np.ndarray:
    """
    Compute cosine similarity between GeminiMol embeddings.
    
    Parameters
    ----------
    query_embeddings : np.ndarray
        Query embeddings
    ref_embeddings : np.ndarray
        Reference embeddings
        
    Returns
    -------
    np.ndarray
        Similarity matrix
    """
    try:
        from sklearn.metrics.pairwise import cosine_similarity
        return cosine_similarity(query_embeddings, ref_embeddings)
    except ImportError:
        # Fallback implementation
        # Normalize embeddings
        query_norm = query_embeddings / np.linalg.norm(query_embeddings, axis=1, keepdims=True)
        ref_norm = ref_embeddings / np.linalg.norm(ref_embeddings, axis=1, keepdims=True)
        
        # Compute cosine similarity
        return np.dot(query_norm, ref_norm.T)
