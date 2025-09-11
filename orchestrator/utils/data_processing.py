"""
Data processing utilities for aggregating, validating, and summarizing
pipeline outputs for reporting and visualization.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd


REQUIRED_RESULT_COLUMNS = [
    "compound_id",
    "target_uniprot",
    "fused_score",
]


@dataclass
class RunData:
    results: pd.DataFrame
    manifest: Optional[dict]
    extras: Dict[str, Path]


def load_results(results_csv: Path) -> pd.DataFrame:
    df = pd.read_csv(results_csv)
    return df


def validate_results(df: pd.DataFrame) -> List[str]:
    errors: List[str] = []
    for col in REQUIRED_RESULT_COLUMNS:
        if col not in df.columns:
            errors.append(f"Missing required column: {col}")
    return errors


def summarize_results(df: pd.DataFrame) -> Dict[str, float]:
    summary: Dict[str, float] = {}
    if "fused_score" in df.columns:
        summary["n_rows"] = float(len(df))
        summary["fused_score_mean"] = float(df["fused_score"].mean())
        summary["fused_score_std"] = float(df["fused_score"].std(ddof=1) if len(df) > 1 else 0.0)
        summary["fused_score_min"] = float(df["fused_score"].min())
        summary["fused_score_max"] = float(df["fused_score"].max())
    if "cosine_max" in df.columns:
        summary["cosine_max_mean"] = float(df["cosine_max"].mean())
    if "tanimoto_max" in df.columns:
        summary["tanimoto_max_mean"] = float(df["tanimoto_max"].mean())
    return summary


def load_manifest(manifest_path: Path) -> Optional[dict]:
    try:
        import json
        with open(manifest_path, "r") as f:
            return json.load(f)
    except Exception:
        return None


def collect_run_data(results_csv: Path, manifest_path: Optional[Path] = None, **extras: Path) -> RunData:
    df = load_results(results_csv)
    manifest = load_manifest(manifest_path) if manifest_path else None
    return RunData(results=df, manifest=manifest, extras=extras)

