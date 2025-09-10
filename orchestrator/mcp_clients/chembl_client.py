"""
ChEMBL comparators client helpers.

These wrappers use the existing MCP `ChEMBLClient` (if configured) to fetch
active ligands per target and to summarize comparator quality characteristics.

All thresholds and gating should be driven by external configuration.
"""

from __future__ import annotations

from statistics import median
from typing import Any, Dict, List

from loguru import logger

from .chembl import ChEMBLClient as MCPChEMBLClient


def fetch_target_actives(chembl_target_id: str, min_pchembl: float) -> List[Dict]:
    """Fetch active compounds for a ChEMBL target via MCP.

    Returns records with keys: {chembl_id, smiles, pchembl_value, assay_type, organism}.
    If the MCP is offline, returns an empty list (do not fabricate data).
    """
    client = MCPChEMBLClient()
    try:
        # Underlying MCP returns activities; adapt and filter by pChEMBL
        # The existing MCP client exposes `get_bioactivity_data(target_id, compound_id=None)`
        # Here we only pass target_id
        # Note: The exact response schema depends on the MCP; adapt mapping as needed.
        loop = None
        try:
            import asyncio
            loop = asyncio.get_event_loop()
        except RuntimeError:
            pass

        async def _fetch():
            return await client.get_bioactivity_data(target_id=chembl_target_id)

        if loop and loop.is_running():
            # If already running in async context, caller should integrate async.
            logger.warning("fetch_target_actives called in running event loop; returning empty.")
            return []
        else:
            import asyncio
            activities = asyncio.run(_fetch())

        records: List[Dict[str, Any]] = []
        for a in activities or []:
            pv = a.get("pchembl_value") or a.get("pChEMBL_value")
            try:
                pv = float(pv) if pv is not None else None
            except Exception:
                pv = None
            if pv is None or pv < float(min_pchembl):
                continue
            rec = {
                "chembl_id": a.get("molecule_chembl_id") or a.get("chembl_id"),
                "smiles": a.get("canonical_smiles") or a.get("smiles"),
                "pchembl_value": pv,
                "assay_type": a.get("assay_type"),
                "organism": a.get("target_organism") or a.get("organism"),
            }
            # Only accept records with a SMILES string
            if rec["smiles"]:
                records.append(rec)
        return records
    except Exception as e:
        logger.warning(f"ChEMBL MCP unavailable or error: {e}")
        return []


def summarize_comparator_quality(records: List[Dict]) -> Dict:
    """Summarize comparator set quality into 0–1 metrics.

    Computes size, median pChEMBL, and a crude assay consistency score in [0,1].
    The latter penalizes heterogeneous assay_type distributions.
    """
    size = len(records)
    pvals = [r.get("pchembl_value") for r in records if r.get("pchembl_value") is not None]
    med = float(median(pvals)) if pvals else 0.0

    # Assay consistency: 1.0 when all assays are same type; lower otherwise.
    if size <= 1:
        assay_consistency = 1.0 if size == 1 else 0.0
    else:
        types = [r.get("assay_type") or "UNK" for r in records]
        total = float(len(types))
        max_frac = 0.0
        if total > 0:
            from collections import Counter
            c = Counter(types)
            max_frac = max(c.values()) / total
        assay_consistency = float(max_frac)

    return {
        "size": size,
        "median_pchembl": med,
        "assay_consistency_score": assay_consistency,
    }

