"""
Export utilities for saving figures and data in multiple formats.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import pandas as pd


def export_dataframe(df: pd.DataFrame, out_dir: Path, basename: str) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = {}
    csv_path = out_dir / f"{basename}.csv"
    json_path = out_dir / f"{basename}.json"
    df.to_csv(csv_path, index=False)
    df.to_json(json_path, orient="records")
    paths["csv"] = csv_path
    paths["json"] = json_path
    try:
        parquet_path = out_dir / f"{basename}.parquet"
        df.to_parquet(parquet_path, index=False)
        paths["parquet"] = parquet_path
    except Exception:
        pass
    return paths


def save_figure(fig, out_path: Path, formats: Iterable[str] = ("png", "svg"), dpi: int = 150) -> dict:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    saved = {}
    for fmt in formats:
        path = out_path.with_suffix(f".{fmt}")
        try:
            fig.savefig(path, dpi=dpi, bbox_inches="tight")
            saved[fmt] = path
        except Exception:
            continue
    return saved

