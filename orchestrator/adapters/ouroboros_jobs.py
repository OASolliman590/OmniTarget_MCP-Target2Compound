"""
Ouroboros adapter for Chemical modes: Check, Exploration, Migration, and Fusion.

Ouroboros provides directed chemical evolution through various modes:
- ChemicalCheck: reconstruction similarity testing
- ChemicalExploration: scaffold hopping and optimization
- ChemicalMigration: directed evolution between molecules
- ChemicalFusion: chimera generation

References:
- Official page: https://www.aideepmed.com/Ouroboros/
- Preprint: Directed Chemical Evolution via Navigating Molecular Space
"""

from __future__ import annotations

import json
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional
from loguru import logger


def _check_ouroboros_setup(repo_path: str) -> tuple[bool, str]:
    """
    Check if Ouroboros is properly set up.
    
    Parameters
    ----------
    repo_path : str
        Path to Ouroboros repository
        
    Returns
    -------
    Tuple[bool, str]
        (is_setup, error_message)
    """
    ouroboros_path = Path(repo_path)
    
    if not ouroboros_path.exists():
        return False, f"Ouroboros repository not found at {repo_path}"
    
    # Check for required files/directories. Upstream repo may place scripts under
    # a top-level 'scripts/' dir or under 'ouroboros/'. Accept either layout.
    if not (ouroboros_path / "README.md").exists():
        return False, "Required README.md not found in Ouroboros repository"

    scripts_root = None
    if (ouroboros_path / "scripts").exists():
        scripts_root = ouroboros_path / "scripts"
    elif (ouroboros_path / "ouroboros").exists():
        scripts_root = ouroboros_path / "ouroboros"

    if scripts_root is None:
        return False, "No 'scripts' or 'ouroboros' directory found in Ouroboros repository"

    # Verify at least one Chemical*.py exists
    expected = [
        "ChemicalCheck.py",
        "ChemicalExploration.py",
        "ChemicalMigration.py",
        "ChemicalFusion.py",
    ]
    if not any((scripts_root / name).exists() for name in expected):
        return False, f"No Chemical*.py scripts found under {scripts_root}"

    return True, ""


def run_chemical_check(
    start_smiles: str,
    model_name: str,
    job_name: str,
    repo_path: str,
    output_dir: Path,
    **kwargs
) -> Dict[str, Any]:
    """
    Run ChemicalCheck mode to test reconstruction similarity.
    
    Parameters
    ----------
    start_smiles : str
        Starting SMILES string
    model_name : str
        Ouroboros model name
    job_name : str
        Job identifier
    repo_path : str
        Path to Ouroboros repository
    output_dir : Path
        Output directory for results
    **kwargs
        Additional parameters
        
    Returns
    -------
    Dict[str, Any]
        Job results and metadata
    """
    try:
        # Check setup
        is_setup, error_msg = _check_ouroboros_setup(repo_path)
        if not is_setup:
            logger.warning(f"Ouroboros setup check failed: {error_msg}")
            return {"error": error_msg, "success": False}
        
        # Create output directory
        job_output_dir = output_dir / f"ouroboros_{job_name}"
        job_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create input file
        input_file = job_output_dir / "input.smi"
        with open(input_file, 'w') as f:
            f.write(f"{start_smiles}\n")
        
        # Output file
        output_file = job_output_dir / f"{job_name}_generation.csv"
        
        # Construct command (placeholder - actual command depends on Ouroboros API)
        # Resolve script path across possible layouts
        script_candidates = [
            Path(repo_path) / "scripts" / "ChemicalCheck.py",
            Path(repo_path) / "ChemicalCheck.py",
            Path(repo_path) / "ouroboros" / "ChemicalCheck.py",
        ]
        script_path = next((p for p in script_candidates if p.exists()), None)
        
        if script_path is None:
            logger.warning(f"ChemicalCheck script not found in {repo_path}")
            return {"error": "ChemicalCheck script not found", "success": False}
        
        cmd = [
            "python",
            str(script_path),
            "--input", str(input_file),
            "--output", str(output_file),
            "--model", model_name,
            "--start_smiles", start_smiles
        ]
        
        # Add additional parameters
        for key, value in kwargs.items():
            if value is not None:
                cmd.extend([f"--{key}", str(value)])
        
        logger.info(f"Running ChemicalCheck: {' '.join(cmd)}")
        
        # Run command
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=1800,  # 30 minute timeout
            cwd=repo_path
        )
        
        # Parse results
        results = {
            "job_type": "chemical_check",
            "job_name": job_name,
            "start_smiles": start_smiles,
            "model_name": model_name,
            "command": " ".join(cmd),
            "success": result.returncode == 0,
            "stdout": result.stdout,
            "stderr": result.stderr,
            "output_file": str(output_file) if output_file.exists() else None
        }
        
        if result.returncode == 0 and output_file.exists():
            # Parse CSV output
            parsed_results = _parse_generation_output(output_file)
            results.update(parsed_results)
            logger.info(f"ChemicalCheck completed successfully: {job_name}")
        else:
            logger.warning(f"ChemicalCheck failed: {result.stderr}")
            results["error"] = result.stderr
        
        return results
        
    except Exception as e:
        logger.error(f"ChemicalCheck failed: {e}")
        return {"error": str(e), "success": False}


def run_chemical_exploration(
    start_smiles: str,
    mode: str,
    model_name: str,
    job_name: str,
    repo_path: str,
    output_dir: Path,
    replica_num: int = 1,
    steps: int = 100,
    step_interval: int = 10,
    loud: bool = False,
    temperature: float = 1.0,
    learning_rate: float = 0.01,
    **kwargs
) -> Dict[str, Any]:
    """
    Run ChemicalExploration mode for scaffold hopping and optimization.
    
    Parameters
    ----------
    start_smiles : str
        Starting SMILES string
    mode : str
        Exploration mode (scaffold_hopping, directional_optimization, etc.)
    model_name : str
        Ouroboros model name
    job_name : str
        Job identifier
    repo_path : str
        Path to Ouroboros repository
    output_dir : Path
        Output directory for results
    replica_num : int
        Number of replicas
    steps : int
        Number of steps
    step_interval : int
        Step interval for output
    loud : bool
        Verbose output
    temperature : float
        Sampling temperature
    learning_rate : float
        Learning rate
    **kwargs
        Additional parameters
        
    Returns
    -------
    Dict[str, Any]
        Job results and metadata
    """
    try:
        # Check setup
        is_setup, error_msg = _check_ouroboros_setup(repo_path)
        if not is_setup:
            logger.warning(f"Ouroboros setup check failed: {error_msg}")
            return {"error": error_msg, "success": False}
        
        # Create output directory
        job_output_dir = output_dir / f"ouroboros_{job_name}"
        job_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create input file
        input_file = job_output_dir / "input.smi"
        with open(input_file, 'w') as f:
            f.write(f"{start_smiles}\n")
        
        # Output file
        output_file = job_output_dir / f"{job_name}_exploration.csv"
        
        # Construct command
        script_candidates = [
            Path(repo_path) / "scripts" / "ChemicalExploration.py",
            Path(repo_path) / "ChemicalExploration.py",
            Path(repo_path) / "ouroboros" / "ChemicalExploration.py",
        ]
        script_path = next((p for p in script_candidates if p.exists()), None)

        if script_path is None:
            logger.warning(f"ChemicalExploration script not found in {repo_path}")
            return {"error": "ChemicalExploration script not found", "success": False}
        
        cmd = [
            "python",
            str(script_path),
            "--input", str(input_file),
            "--output", str(output_file),
            "--model", model_name,
            "--mode", mode,
            "--start_smiles", start_smiles,
            "--replica_num", str(replica_num),
            "--steps", str(steps),
            "--step_interval", str(step_interval),
            "--temperature", str(temperature),
            "--learning_rate", str(learning_rate)
        ]
        
        if loud:
            cmd.append("--loud")
        
        # Add additional parameters
        for key, value in kwargs.items():
            if value is not None:
                cmd.extend([f"--{key}", str(value)])
        
        logger.info(f"Running ChemicalExploration: {' '.join(cmd)}")
        
        # Run command
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600,  # 60 minute timeout
            cwd=repo_path
        )
        
        # Parse results
        results = {
            "job_type": "chemical_exploration",
            "job_name": job_name,
            "start_smiles": start_smiles,
            "mode": mode,
            "model_name": model_name,
            "command": " ".join(cmd),
            "success": result.returncode == 0,
            "stdout": result.stdout,
            "stderr": result.stderr,
            "output_file": str(output_file) if output_file.exists() else None,
            "parameters": {
                "replica_num": replica_num,
                "steps": steps,
                "step_interval": step_interval,
                "temperature": temperature,
                "learning_rate": learning_rate
            }
        }
        
        if result.returncode == 0 and output_file.exists():
            # Parse CSV output
            parsed_results = _parse_generation_output(output_file)
            results.update(parsed_results)
            logger.info(f"ChemicalExploration completed successfully: {job_name}")
        else:
            logger.warning(f"ChemicalExploration failed: {result.stderr}")
            results["error"] = result.stderr
        
        return results
        
    except Exception as e:
        logger.error(f"ChemicalExploration failed: {e}")
        return {"error": str(e), "success": False}


def run_chemical_migration(
    start_smiles: str,
    ref_smiles: str,
    model_name: str,
    job_name: str,
    repo_path: str,
    output_dir: Path,
    **kwargs
) -> Dict[str, Any]:
    """
    Run ChemicalMigration mode for directed evolution between molecules.
    
    Parameters
    ----------
    start_smiles : str
        Starting SMILES string
    ref_smiles : str
        Reference SMILES string
    model_name : str
        Ouroboros model name
    job_name : str
        Job identifier
    repo_path : str
        Path to Ouroboros repository
    output_dir : Path
        Output directory for results
    **kwargs
        Additional parameters
        
    Returns
    -------
    Dict[str, Any]
        Job results and metadata
    """
    try:
        # Check setup
        is_setup, error_msg = _check_ouroboros_setup(repo_path)
        if not is_setup:
            logger.warning(f"Ouroboros setup check failed: {error_msg}")
            return {"error": error_msg, "success": False}
        
        # Create output directory
        job_output_dir = output_dir / f"ouroboros_{job_name}"
        job_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create input file
        input_file = job_output_dir / "input.smi"
        with open(input_file, 'w') as f:
            f.write(f"{start_smiles}\n{ref_smiles}\n")
        
        # Output file
        output_file = job_output_dir / f"{job_name}_migration.csv"
        
        # Construct command
        script_candidates = [
            Path(repo_path) / "scripts" / "ChemicalMigration.py",
            Path(repo_path) / "ChemicalMigration.py",
            Path(repo_path) / "ouroboros" / "ChemicalMigration.py",
        ]
        script_path = next((p for p in script_candidates if p.exists()), None)

        if script_path is None:
            logger.warning(f"ChemicalMigration script not found in {repo_path}")
            return {"error": "ChemicalMigration script not found", "success": False}
        
        cmd = [
            "python",
            str(script_path),
            "--input", str(input_file),
            "--output", str(output_file),
            "--model", model_name,
            "--start_smiles", start_smiles,
            "--ref_smiles", ref_smiles
        ]
        
        # Add additional parameters
        for key, value in kwargs.items():
            if value is not None:
                cmd.extend([f"--{key}", str(value)])
        
        logger.info(f"Running ChemicalMigration: {' '.join(cmd)}")
        
        # Run command
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600,  # 60 minute timeout
            cwd=repo_path
        )
        
        # Parse results
        results = {
            "job_type": "chemical_migration",
            "job_name": job_name,
            "start_smiles": start_smiles,
            "ref_smiles": ref_smiles,
            "model_name": model_name,
            "command": " ".join(cmd),
            "success": result.returncode == 0,
            "stdout": result.stdout,
            "stderr": result.stderr,
            "output_file": str(output_file) if output_file.exists() else None
        }
        
        if result.returncode == 0 and output_file.exists():
            # Parse CSV output
            parsed_results = _parse_generation_output(output_file)
            results.update(parsed_results)
            logger.info(f"ChemicalMigration completed successfully: {job_name}")
        else:
            logger.warning(f"ChemicalMigration failed: {result.stderr}")
            results["error"] = result.stderr
        
        return results
        
    except Exception as e:
        logger.error(f"ChemicalMigration failed: {e}")
        return {"error": str(e), "success": False}


def run_chemical_fusion(
    probe_dataset_csv: str,
    fusion_targets: List[str],
    model_name: str,
    job_name: str,
    repo_path: str,
    output_dir: Path,
    fusion_temperature: float = 1.0,
    **kwargs
) -> Dict[str, Any]:
    """
    Run ChemicalFusion mode for chimera generation.
    
    Parameters
    ----------
    probe_dataset_csv : str
        Path to probe dataset CSV
    fusion_targets : List[str]
        List of fusion targets (e.g., ["AURKA:PI3Kg"])
    model_name : str
        Ouroboros model name
    job_name : str
        Job identifier
    repo_path : str
        Path to Ouroboros repository
    output_dir : Path
        Output directory for results
    fusion_temperature : float
        Fusion temperature
    **kwargs
        Additional parameters
        
    Returns
    -------
    Dict[str, Any]
        Job results and metadata
    """
    try:
        # Check setup
        is_setup, error_msg = _check_ouroboros_setup(repo_path)
        if not is_setup:
            logger.warning(f"Ouroboros setup check failed: {error_msg}")
            return {"error": error_msg, "success": False}
        
        # Create output directory
        job_output_dir = output_dir / f"ouroboros_{job_name}"
        job_output_dir.mkdir(parents=True, exist_ok=True)
        
        # Copy probe dataset
        probe_file = job_output_dir / "probe_dataset.csv"
        import shutil
        shutil.copy2(probe_dataset_csv, probe_file)
        
        # Output file
        output_file = job_output_dir / f"{job_name}_fusion.csv"
        
        # Construct command
        script_candidates = [
            Path(repo_path) / "scripts" / "ChemicalFusion.py",
            Path(repo_path) / "ChemicalFusion.py",
            Path(repo_path) / "ouroboros" / "ChemicalFusion.py",
        ]
        script_path = next((p for p in script_candidates if p.exists()), None)

        if script_path is None:
            logger.warning(f"ChemicalFusion script not found in {repo_path}")
            return {"error": "ChemicalFusion script not found", "success": False}
        
        cmd = [
            "python",
            str(script_path),
            "--probe_dataset", str(probe_file),
            "--output", str(output_file),
            "--model", model_name,
            "--fusion_temperature", str(fusion_temperature)
        ]
        
        # Add fusion targets
        for target in fusion_targets:
            cmd.extend(["--fusion_target", target])
        
        # Add additional parameters
        for key, value in kwargs.items():
            if value is not None:
                cmd.extend([f"--{key}", str(value)])
        
        logger.info(f"Running ChemicalFusion: {' '.join(cmd)}")
        
        # Run command
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600,  # 60 minute timeout
            cwd=repo_path
        )
        
        # Parse results
        results = {
            "job_type": "chemical_fusion",
            "job_name": job_name,
            "probe_dataset": probe_dataset_csv,
            "fusion_targets": fusion_targets,
            "model_name": model_name,
            "command": " ".join(cmd),
            "success": result.returncode == 0,
            "stdout": result.stdout,
            "stderr": result.stderr,
            "output_file": str(output_file) if output_file.exists() else None,
            "parameters": {
                "fusion_temperature": fusion_temperature
            }
        }
        
        if result.returncode == 0 and output_file.exists():
            # Parse CSV output
            parsed_results = _parse_generation_output(output_file)
            results.update(parsed_results)
            logger.info(f"ChemicalFusion completed successfully: {job_name}")
        else:
            logger.warning(f"ChemicalFusion failed: {result.stderr}")
            results["error"] = result.stderr
        
        return results
        
    except Exception as e:
        logger.error(f"ChemicalFusion failed: {e}")
        return {"error": str(e), "success": False}


def _parse_generation_output(output_file: Path) -> Dict[str, Any]:
    """
    Parse Ouroboros generation output CSV file.
    
    Parameters
    ----------
    output_file : Path
        Path to output CSV file
        
    Returns
    -------
    Dict[str, Any]
        Parsed results
    """
    try:
        import pandas as pd
        
        if not output_file.exists():
            return {"error": "Output file not found"}
        
        # Read CSV file
        df = pd.read_csv(output_file)
        
        # Extract key metrics
        results = {
            "n_generated": len(df),
            "columns": list(df.columns)
        }
        
        # Look for similarity columns
        similarity_cols = [col for col in df.columns if 'similarity' in col.lower()]
        if similarity_cols:
            results["max_similarity"] = float(df[similarity_cols[0]].max())
            results["mean_similarity"] = float(df[similarity_cols[0]].mean())
        
        # Look for SMILES column
        smiles_cols = [col for col in df.columns if 'smiles' in col.lower()]
        if smiles_cols:
            results["generated_smiles"] = df[smiles_cols[0]].tolist()
        
        # Look for property columns
        property_cols = [col for col in df.columns if col not in similarity_cols + smiles_cols]
        if property_cols:
            results["properties"] = {col: df[col].tolist() for col in property_cols}
        
        return results
        
    except Exception as e:
        logger.warning(f"Failed to parse Ouroboros output: {e}")
        return {"error": str(e)}


def run_ouroboros_jobs(
    jobs: List[Dict[str, Any]],
    repo_path: str,
    output_dir: Path
) -> Dict[str, Any]:
    """
    Run multiple Ouroboros jobs.
    
    Parameters
    ----------
    jobs : List[Dict[str, Any]]
        List of job configurations
    repo_path : str
        Path to Ouroboros repository
    output_dir : Path
        Output directory for results
        
    Returns
    -------
    Dict[str, Any]
        Summary of all job results
    """
    results = {
        "total_jobs": len(jobs),
        "successful_jobs": 0,
        "failed_jobs": 0,
        "job_results": []
    }
    
    for job_config in jobs:
        job_type = job_config.get("type")
        job_name = job_config.get("job_name", f"job_{len(results['job_results'])}")
        
        try:
            if job_type == "check":
                result = run_chemical_check(
                    start_smiles=job_config["start_smiles"],
                    model_name=job_config["model_name"],
                    job_name=job_name,
                    repo_path=repo_path,
                    output_dir=output_dir,
                    **{k: v for k, v in job_config.items() 
                       if k not in ["type", "start_smiles", "model_name", "job_name"]}
                )
            elif job_type == "exploration":
                result = run_chemical_exploration(
                    start_smiles=job_config["start_smiles"],
                    mode=job_config["mode"],
                    model_name=job_config["model_name"],
                    job_name=job_name,
                    repo_path=repo_path,
                    output_dir=output_dir,
                    **{k: v for k, v in job_config.items() 
                       if k not in ["type", "start_smiles", "mode", "model_name", "job_name"]}
                )
            elif job_type == "migration":
                result = run_chemical_migration(
                    start_smiles=job_config["start_smiles"],
                    ref_smiles=job_config["ref_smiles"],
                    model_name=job_config["model_name"],
                    job_name=job_name,
                    repo_path=repo_path,
                    output_dir=output_dir,
                    **{k: v for k, v in job_config.items() 
                       if k not in ["type", "start_smiles", "ref_smiles", "model_name", "job_name"]}
                )
            elif job_type == "fusion":
                result = run_chemical_fusion(
                    probe_dataset_csv=job_config["probe_dataset_csv"],
                    fusion_targets=job_config["fusion_targets"],
                    model_name=job_config["model_name"],
                    job_name=job_name,
                    repo_path=repo_path,
                    output_dir=output_dir,
                    **{k: v for k, v in job_config.items() 
                       if k not in ["type", "probe_dataset_csv", "fusion_targets", "model_name", "job_name"]}
                )
            else:
                result = {
                    "error": f"Unknown job type: {job_type}",
                    "success": False
                }
            
            results["job_results"].append(result)
            
            if result.get("success", False):
                results["successful_jobs"] += 1
            else:
                results["failed_jobs"] += 1
                
        except Exception as e:
            error_result = {
                "job_name": job_name,
                "job_type": job_type,
                "error": str(e),
                "success": False
            }
            results["job_results"].append(error_result)
            results["failed_jobs"] += 1
    
    # Save summary
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        summary_file = output_dir / "ouroboros_jobs_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(results, f, indent=2)
    except Exception as e:
        logger.warning(f"Failed to save Ouroboros summary: {e}")
    
    logger.info(f"Ouroboros jobs completed: {results['successful_jobs']}/{results['total_jobs']} successful")
    
    return results
