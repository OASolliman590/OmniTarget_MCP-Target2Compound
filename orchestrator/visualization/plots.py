"""
Static plot generation for pipeline results using matplotlib/seaborn.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

sns.set_context("talk")


def plot_score_distribution(df: pd.DataFrame, score_col: str = "fused_score"):
    fig, ax = plt.subplots(figsize=(8, 5))
    if score_col in df.columns:
        sns.histplot(df[score_col], kde=True, ax=ax)
        ax.set_title("Fused Score Distribution")
        ax.set_xlabel(score_col)
    else:
        ax.text(0.5, 0.5, f"Missing column: {score_col}", ha="center", va="center")
        ax.set_axis_off()
    return fig


def plot_top_compounds(df: pd.DataFrame, k: int = 20):
    fig, ax = plt.subplots(figsize=(10, max(4, k * 0.3)))
    cols = [c for c in ("compound_id", "fused_score") if c in df.columns]
    if len(cols) == 2:
        top = df.sort_values("fused_score", ascending=False).head(k)
        sns.barplot(y=top["compound_id"], x=top["fused_score"], ax=ax, orient="h")
        ax.set_title(f"Top {len(top)} Compounds by Fused Score")
        ax.set_xlabel("Fused Score")
        ax.set_ylabel("Compound")
    else:
        ax.text(0.5, 0.5, "Missing columns for top compounds", ha="center", va="center")
        ax.set_axis_off()
    return fig


def plot_similarity_vs_score(df: pd.DataFrame, sim_col: str = "cosine_max", score_col: str = "fused_score"):
    fig, ax = plt.subplots(figsize=(6, 5))
    if sim_col in df.columns and score_col in df.columns:
        sns.scatterplot(x=df[sim_col], y=df[score_col], s=30, ax=ax)
        ax.set_title("Similarity vs Fused Score")
        ax.set_xlabel(sim_col)
        ax.set_ylabel(score_col)
    else:
        ax.text(0.5, 0.5, "Missing columns for scatter", ha="center", va="center")
        ax.set_axis_off()
    return fig

