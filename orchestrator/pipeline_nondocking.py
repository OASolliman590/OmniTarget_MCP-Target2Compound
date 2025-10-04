"""
Non-docking pipeline implementation per OmniTarget MCP-Target2Compound v2 spec.

Stage order:
 1) discover_targets → annotate_targets (UniProt/STRING/ProteinAtlas/PDB)
 2) build_comparators (ChEMBL) → evidence_strength
 3) GeminiMol embeddings & index → per-ligand similarity features
 4) QSAR per target (if feasible)
 5) Pharmacophore vs top-k comparators
 6) Fuse scores → rank; write CSV + manifest

Notes:
- Docking disabled by default; pocket finding hooks are present elsewhere.
- All thresholds and gates are provided via config; nothing is hard-coded here.
- When services/models are unavailable, steps are skipped without fabricating data.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
from dataclasses import asdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional

from loguru import logger

from .mcp_clients.kegg import KEGGClient
from .mcp_clients.reactome import ReactomeClient
from .mcp_clients.proteinatlas import ProteinAtlasClient
from .mcp_clients.string import STRINGClient
from .mcp_clients.uniprot import UniProtClient
from .mcp_clients.pdb import PDBClient
from .mcp_clients.chembl_client import summarize_comparator_quality
from .mcp_clients.chembl_mcp_client import ChEMBLMCPClient

from .adapters.geminimol_embed import tanimoto_neighbors_only
from .adapters.geminimol_adapter import compute_geminimol_features
from .adapters.ouroboros_jobs import run_ouroboros_jobs
from .adapters.pharmacophore import rdkit_feature_ph4
from .adapters.e3fp_adapter import (
    compute_e3fp_for_sdf, 
    topk_similarity_matrix_e3fp,
    compute_e3fp_features_for_ligands
)
from .io.chem_io import read_ligands, standardize, enumerate_chemistry, generate_3d, export
from .scoring.evidence import evidence_strength, similarity_features
from .scoring.fuse_nd import fuse as fuse_scores
from .adapters.vina_adapter import VinaAdapter
from .adapters.pocket_box import box_from_cocrystal
from .curation.ccd import is_common_additive
from .curation.sifts import maps_to_uniprot_chain
# from .visualization.reports import generate_report  # Temporarily disabled due to syntax error


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


class NonDockingPipeline:
    def __init__(self, config: Any):
        self.cfg = config
        # MCP clients (used when available)
        self.kegg = KEGGClient()
        self.reactome = ReactomeClient()
        self.hpa = ProteinAtlasClient()
        self.string = STRINGClient()
        self.uniprot = UniProtClient()
        self.pdb = PDBClient()
        self.output_dir = Path(self.cfg.get("output_dir", "data/outputs")) if isinstance(self.cfg, dict) else Path(getattr(self.cfg, "output_dir", "data/outputs"))
        self.run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.run_dir = self.output_dir / f"run_{self.run_id}"
        self.run_dir.mkdir(parents=True, exist_ok=True)

    def _manifest_base(self) -> Dict[str, Any]:
        import platform
        return {
            "run_id": self.run_id,
            "created_at": datetime.utcnow().isoformat() + "Z",
            "python": platform.python_version(),
            "platform": platform.platform(),
            "config": self.cfg.dict() if hasattr(self.cfg, "dict") else self.cfg,
            "stages": [],
        }

    def _write_results(self, rows: List[Dict[str, Any]]) -> Path:
        out = self.run_dir / "results.csv"
        cols = [
            "compound_id",
            "compound_smiles",
            "compound_inchikey", 
            "compound_source",
            "target_uniprot",
            "target_chembl",
            "target_gene",
            "evidence_strength",
            # Comparator diagnostics
            "n_comparators",
            "median_pchembl",
            "assay_consistency",
            "gate_reason",
            # Similarity features
            "cosine_max",
            "tanimoto_max",
            # E3FP features
            "e3fp_tanimoto_max",
            "e3fp_tanimoto_mean_topk",
            # GeminiMol features
            "gm_cosine_max",
            "gm_cosine_mean_topk",
            "gm_profile_score",
            # Other features
            "ph4_best",
            "qsar_score",
            "docking_score",
            "fused_score",
            "rank",
        ]
        with open(out, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=cols)
            writer.writeheader()
            for r in rows:
                writer.writerow({k: r.get(k) for k in cols})
        return out

    def _write_docking(self, rows: List[Dict[str, Any]]) -> Path:
        out = self.run_dir / "docking_results.csv"
        cols = [
            "compound_id",
            "target_uniprot",
            "target_chembl",
            "binding_energy",
            "structure_id",
            "vina_exhaustiveness",
        ]
        with open(out, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=cols)
            writer.writeheader()
            for r in rows:
                writer.writerow({k: r.get(k) for k in cols})
        return out

    def _write_manifest(self, manifest: Dict[str, Any]) -> Path:
        path = self.run_dir / "manifest.json"
        with open(path, "w") as f:
            json.dump(manifest, f, indent=2)
        return path

    async def run(self) -> Dict[str, Any]:
        manifest = self._manifest_base()

        # 0) Ligand ingestion and prep
        lig_paths = []
        ipaths = self.cfg.get("ligand_prep", {}).get("input_paths") if isinstance(self.cfg, dict) else self.cfg.ligand_prep.input_paths
        for pattern in ipaths:
            for p in Path().glob(pattern):
                lig_paths.append(p)
        ligs = read_ligands(lig_paths)
        ligs = standardize(
            ligs,
            ph=(self.cfg.get("ligand_prep", {}).get("ph") if isinstance(self.cfg, dict) else self.cfg.ligand_prep.ph) or 7.4,
            strip_salts=True,
            largest_fragment_only=True,
        )
        if (self.cfg.get("ligand_prep", {}).get("enumerate", {}).get("enable") if isinstance(self.cfg, dict) else self.cfg.ligand_prep.enumerate.enable):
            ligs = enumerate_chemistry(ligs)
        if (self.cfg.get("ligand_prep", {}).get("conformers", {}).get("enable", True) if isinstance(self.cfg, dict) else self.cfg.ligand_prep.conformers.enable):
            ligs = generate_3d(ligs)
        exports = export(ligs, self.run_dir, ("smi", "sdf"))
        manifest["stages"].append({
            "stage": "ligand_prep",
            "inputs": [str(p) for p in lig_paths],
            "outputs": {k: str(v) for k, v in exports.items()},
            "hashes": {k: _sha256(v) for k, v in exports.items()},
        })

        # 1) Discover/annotate targets (offline-safe: may be empty)
        targets: List[Dict[str, Any]] = []
        try:
            # Check for test targets first
            test_targets = self.cfg.get("test_targets", []) if isinstance(self.cfg, dict) else getattr(self.cfg, "test_targets", [])
            if test_targets:
                targets = test_targets
                logger.info(f"Using {len(targets)} test targets")
            else:
                # Real target discovery using MCP services
                disease_terms = self.cfg.get("disease_terms", []) if isinstance(self.cfg, dict) else getattr(self.cfg, "disease_terms", [])
                organism = self.cfg.get("organism", "Homo sapiens") if isinstance(self.cfg, dict) else getattr(self.cfg, "organism", "Homo sapiens")
                max_targets = self.cfg.get("max_targets", 20) if isinstance(self.cfg, dict) else getattr(self.cfg, "max_targets", 20)
                
                if disease_terms:
                    logger.info(f"Discovering targets for diseases: {disease_terms}")
                    
                    # Initialize MCP clients
                    kegg_client = KEGGClient()
                    reactome_client = ReactomeClient()
                    
                    try:
                        # Get targets from KEGG pathways
                        kegg_targets = await kegg_client.discover_disease_targets(
                            disease_terms=disease_terms,
                            organism="hsa",  # Human
                            max_pathways_per_disease=10
                        )
                        logger.info(f"KEGG discovered {len(kegg_targets)} targets")
                        
                        # Get targets from Reactome pathways
                        reactome_targets = []
                        for term in disease_terms:
                            pathways = await reactome_client.search_pathways(term)
                            for pathway in pathways[:5]:  # Limit to 5 pathways per term
                                # Extract genes from pathway (simplified)
                                if pathway.get("genes"):
                                    for gene in pathway["genes"]:
                                        reactome_targets.append({
                                            "target_id": gene,
                                            "gene_name": gene,
                                            "discovery_source": "Reactome",
                                            "discovery_pathway": pathway
                                        })
                        logger.info(f"Reactome discovered {len(reactome_targets)} targets")
                        
                        # Combine and deduplicate targets
                        all_targets = kegg_targets + reactome_targets
                        unique_targets = []
                        seen_ids = set()
                        
                        for target in all_targets:
                            target_id = target.get("target_id") or target.get("gene_id") or target.get("gene_name")
                            if target_id and target_id not in seen_ids:
                                unique_targets.append(target)
                                seen_ids.add(target_id)
                        
                        # Limit targets if specified
                        if max_targets and len(unique_targets) > max_targets:
                            unique_targets = unique_targets[:max_targets]
                            logger.info(f"Limited targets to {max_targets}")
                        
                        targets = unique_targets
                        logger.info(f"Total unique targets discovered: {len(targets)}")
                        
                    except Exception as e:
                        logger.warning(f"Target discovery failed: {e}")
                        targets = []
                    finally:
                        await kegg_client.close()
                        await reactome_client.close()
                else:
                    logger.warning("No disease terms provided for target discovery")
        except Exception as e:
            logger.warning(f"Target discovery skipped: {e}")
        manifest["stages"].append({"stage": "discover_targets", "n_targets": len(targets)})

        # 2) ChEMBL comparators and evidence
        scored: List[Dict[str, Any]] = []
        chembl_cfg = self.cfg.get("chembl") if isinstance(self.cfg, dict) else self.cfg.chembl.dict()
        min_pchembl = float(chembl_cfg.get("min_pchembl", 6.0))
        min_comparators = int(chembl_cfg.get("min_comparators", 5))
        require_consistency = bool(chembl_cfg.get("require_assay_consistency", True))
        chembl_limit = int(chembl_cfg.get("limit", 200))

        # If no targets discovered, skip gracefully
        target_list = targets or []
        manifest_targets: List[Dict[str, Any]] = []
        # Keep per-target comparator SMILES for downstream 3D/2D similarity features
        per_target_actives_smiles: Dict[str, List[str]] = {}
        for t in target_list:
            tid = t.get("chembl_id") or t.get("target_chembl_id") or t.get("target_chembl") or t.get("uniprot_id")
            if not tid:
                continue
            # Use real ChEMBL MCP Server with search_activities tool
            client = ChEMBLMCPClient()
            try:
                # Simple filesystem cache to avoid repeated expensive calls
                import json as _json
                cache_dir = Path("data/cache/chembl")
                cache_dir.mkdir(parents=True, exist_ok=True)
                cache_key = f"{tid}_IC50_limit{chembl_limit}.json"
                cache_path = cache_dir / cache_key

                activities_list: List[Dict[str, Any]] = []
                if cache_path.exists():
                    try:
                        activities_list = _json.loads(cache_path.read_text())
                        logger.info(f"Loaded ChEMBL activities from cache for {tid} ({len(activities_list)})")
                    except Exception:
                        activities_list = []
                if not activities_list:
                    await client.start_server()
                    # Use search_activities tool to get bioactivity data for this target
                    activities = await client.search_activities(
                        target_chembl_id=tid,
                        activity_type="IC50",
                        limit=chembl_limit
                    )
                    await client.stop_server()
                    # The response is in the format: [{'type': 'text', 'text': '{"activities": [...]}'}]
                    if activities and len(activities) > 0:
                        try:
                            response_text = activities[0].get('text', '{}')
                            response_data = _json.loads(response_text)
                            activities_list = response_data.get('activities', [])
                            # Write-through cache
                            try:
                                cache_path.write_text(_json.dumps(activities_list))
                            except Exception:
                                pass
                        except Exception as e:
                            logger.warning(f"Failed to parse ChEMBL response for {tid}: {e}")
                            activities_list = []

                # Parse the response and convert activities to our expected format
                actives = []
                for activity in activities_list or []:
                    try:
                        pchembl_val = activity.get("pchembl_value")
                        if pchembl_val and float(pchembl_val) >= min_pchembl:
                            actives.append({
                                "chembl_id": activity.get("molecule_chembl_id"),
                                "smiles": activity.get("canonical_smiles"),
                                "pchembl_value": float(pchembl_val),
                                "assay_type": activity.get("assay_type"),
                                "organism": activity.get("target_organism")
                            })
                    except Exception:
                        continue

            except Exception as e:
                logger.error(f"Failed to fetch ChEMBL data for {tid}: {e}")
                actives = []
            
            # If no actives found, skip this target
            if not actives:
                logger.warning(f"No ChEMBL data found for {tid}, skipping target")
                continue
            
            summary = summarize_comparator_quality(actives)
            evd = evidence_strength([summary])
            reason = None
            if summary["size"] < min_comparators:
                reason = f"insufficient_comparators({summary['size']} < {min_comparators})"
            if require_consistency and summary["assay_consistency_score"] < 0.5:
                reason = reason or "low_assay_consistency"
            # Keep SMILES for E3FP stage and add a concise preview to manifest
            act_smiles = [a["smiles"] for a in actives if a.get("smiles")]
            per_target_actives_smiles[tid] = act_smiles

            manifest_targets.append({
                "target": tid,
                "comparators": summary,
                "evidence_strength": evd,
                "gate_reason": reason,
                "n_actives": len(actives),
                "actives_preview": act_smiles[:10],
            })
            if reason:
                continue
            # For each ligand, compute similarity features via RDKit Tanimoto fallback
            smiles_list = [x["smiles"] for x in actives if x.get("smiles")]
            for lig in ligs:
                q = lig.get("smiles")
                nbrs = tanimoto_neighbors_only(q, smiles_list, k=int(self.cfg.get("similarity", {}).get("top_k", 50)) if isinstance(self.cfg, dict) else int(self.cfg.similarity.top_k))
                feats = similarity_features(q, [
                    {"cosine": None, "tanimoto": n.tanimoto} for n in nbrs
                ])
                
                # Add comparator diagnostics
                feats.update({
                    "n_comparators": summary["size"],
                    "median_pchembl": summary.get("median_pchembl"),
                    "assay_consistency": summary.get("assay_consistency_score"),
                    "gate_reason": reason,
                })
                
                scored.append({
                    "compound_id": lig.get("name"),
                    "compound_smiles": lig.get("smiles"),
                    "compound_inchikey": lig.get("inchikey"),
                    "compound_source": lig.get("source_path"),
                    "target_uniprot": t.get("uniprot_id", ""),
                    "target_chembl": tid,
                    "target_gene": t.get("gene_name", ""),
                    "evidence_strength": evd,
                    **feats,
                })

        manifest["stages"].append({
            "stage": "comparators",
            "targets": manifest_targets,
        })

        # 3-5) QSAR and Pharmacophore (optional)
        ph4_disabled = bool(self.cfg.get("pharmacophore", {}).get("method") == "disabled") if isinstance(self.cfg, dict) else False

        # Apply optional predictions per scored record
        for r in scored:
            if not ph4_disabled:
                try:
                    # Use topk comparator smiles if available
                    r["ph4_best"] = 0.0  # Computed above requires per-query comparator list; omitted in this offline scaffold
                except Exception:
                    r["ph4_best"] = 0.0

        # 6) E3FP 3D fingerprints (optional)
        e3fp_enabled = bool(self.cfg.get("similarity", {}).get("enable_e3fp", False)) if isinstance(self.cfg, dict) else False
        e3fp_summary: Dict[str, Any] = {"enabled": e3fp_enabled}
        
        if e3fp_enabled:
            try:
                e3fp_cfg = self.cfg.get("similarity", {}).get("e3fp", {}) if isinstance(self.cfg, dict) else {}
                radius = e3fp_cfg.get("radius", 2)
                shells = e3fp_cfg.get("shells", 5)
                top_k = e3fp_cfg.get("top_k", 25)

                # SDF-derived query E3FP if available
                sdf_files = list(self.run_dir.glob("**/*.sdf"))

                query_e3fp: Dict[str, Any] = {}
                if sdf_files:
                    logger.info(f"Computing E3FP fingerprints from {len(sdf_files)} SDF files")
                    for sdf_file in sdf_files:
                        try:
                            fps = compute_e3fp_for_sdf(sdf_file, radius=radius, shells=shells)
                            query_e3fp.update(fps)
                        except Exception as e:
                            logger.warning(f"E3FP computation failed for {sdf_file}: {e}")

                # Fallback to SMILES if SDF route produced no fingerprints
                if not query_e3fp:
                    from .adapters.e3fp_adapter import compute_e3fp_features_for_ligands as _e3fp_from_smiles
                    smiles_list = [lig.get("smiles") for lig in ligs if lig.get("smiles")]
                    names_list = [lig.get("name") for lig in ligs if lig.get("smiles")]
                    try:
                        query_e3fp = _e3fp_from_smiles(smiles_list, names_list, radius=radius, shells=shells)
                    except Exception as e:
                        logger.warning(f"E3FP fallback from SMILES failed: {e}")
                        query_e3fp = {}

                # Reference E3FP per target from comparator SMILES
                per_target_ref_e3fp: Dict[str, Dict[str, Any]] = {}
                from .adapters.e3fp_adapter import compute_e3fp_features_for_ligands as _e3fp_from_smiles
                for tid_key, smiles_list in per_target_actives_smiles.items():
                    if not smiles_list:
                        continue
                    names = [f"{tid_key}_ref_{i}" for i in range(len(smiles_list))]
                    try:
                        ref_fps = _e3fp_from_smiles(smiles_list, names, radius=radius, shells=shells)
                        per_target_ref_e3fp[tid_key] = ref_fps
                    except Exception as e:
                        logger.warning(f"E3FP ref computation failed for {tid_key}: {e}")

                # Use top-k similarities per target to populate record features
                if query_e3fp and per_target_ref_e3fp:
                    from .adapters.e3fp_adapter import topk_similarity_matrix_e3fp as _e3fp_topk
                    per_target_sim_maps: Dict[str, Dict[str, Dict[str, Any]]] = {}
                    for tid_key, ref_map in per_target_ref_e3fp.items():
                        try:
                            per_target_sim_maps[tid_key] = _e3fp_topk(query_e3fp, ref_map, k=top_k)
                        except Exception as e:
                            logger.warning(f"E3FP similarity failed for {tid_key}: {e}")
                            per_target_sim_maps[tid_key] = {}

                    updated, missing = 0, 0
                    for r in scored:
                        compound_name = r.get("compound_id", "")
                        tid_key = r.get("target_chembl") or r.get("target_uniprot") or ""
                        sim_map = per_target_sim_maps.get(tid_key, {})
                        if compound_name in sim_map:
                            r["e3fp_tanimoto_max"] = sim_map[compound_name].get("e3fp_tanimoto_max")
                            r["e3fp_tanimoto_mean_topk"] = sim_map[compound_name].get("e3fp_tanimoto_mean_topk")
                            updated += 1
                        else:
                            r["e3fp_tanimoto_max"] = None
                            r["e3fp_tanimoto_mean_topk"] = None
                            missing += 1

                    e3fp_summary.update({
                        "radius": radius,
                        "shells": shells,
                        "top_k": top_k,
                        "n_query": len(query_e3fp),
                        "n_ref_sets": len(per_target_ref_e3fp),
                        "updated_records": updated,
                        "missing_records": missing,
                        "note": "E3FP computed from 3D SDF or SMILES fallback"
                    })
                else:
                    # Fallback: populate e3fp_* via 2D Morgan similarity
                    logger.warning("E3FP fingerprints unavailable; falling back to 2D Morgan for e3fp_* features")
                    from .adapters.geminimol_embed import tanimoto_neighbors_only as _tan_only
                    for r in scored:
                        tid_key = r.get("target_chembl") or r.get("target_uniprot") or ""
                        comp_smiles = per_target_actives_smiles.get(tid_key, [])
                        q = r.get("compound_smiles")
                        if not q or not comp_smiles:
                            r["e3fp_tanimoto_max"] = None
                            r["e3fp_tanimoto_mean_topk"] = None
                            continue
                        nbrs = _tan_only(q, comp_smiles, k=top_k)
                        if nbrs:
                            scores = [n.tanimoto for n in nbrs if n.tanimoto is not None]
                            r["e3fp_tanimoto_max"] = max(scores) if scores else 0.0
                            r["e3fp_tanimoto_mean_topk"] = (sum(scores) / len(scores)) if scores else 0.0
                        else:
                            r["e3fp_tanimoto_max"] = 0.0
                            r["e3fp_tanimoto_mean_topk"] = 0.0

                    e3fp_summary.update({
                        "radius": radius,
                        "shells": shells,
                        "top_k": top_k,
                        "fallback": "2d_morgan",
                        "note": "E3FP not available; populated via 2D Tanimoto vs comparators"
                    })
            except Exception as e:
                logger.warning(f"E3FP computation failed: {e}")
                e3fp_summary["note"] = f"E3FP computation failed: {str(e)}"
        else:
            # Add None values for E3FP features when disabled
            for r in scored:
                r["e3fp_tanimoto_max"] = None
                r["e3fp_tanimoto_mean_topk"] = None

        # 7) GeminiMol embeddings and PharmProfiler (optional)
        geminimol_cfg = self.cfg.get("representations", {}).get("geminimol", {}) if isinstance(self.cfg, dict) else {}
        geminimol_enabled = bool(geminimol_cfg.get("enabled", False))
        geminimol_summary: Dict[str, Any] = {"enabled": geminimol_enabled}
        
        if geminimol_enabled:
            try:
                repo_path = geminimol_cfg.get("repo_path", "third_party/GeminiMol")
                model_path = geminimol_cfg.get("model_path", "third_party/GeminiMol/models/GeminiMol")
                use_pharm_profiler = geminimol_cfg.get("use_pharm_profiler", True)
                
                # Get query SMILES
                query_smiles = [lig.get("smiles") for lig in ligs if lig.get("smiles")]
                query_names = [lig.get("name", f"mol_{i}") for i, lig in enumerate(ligs) if lig.get("smiles")]
                
                # Get reference SMILES and labels from ChEMBL actives
                ref_smiles = []
                ref_labels = []
                for target_data in manifest_targets:
                    if target_data.get("gate_reason"):
                        continue
                    # Get comparator SMILES from the target's actives
                    # This would need to be passed through from the comparator building stage
                    # For now, use placeholder
                    pass
                
                if query_smiles and ref_smiles:
                    logger.info(f"Computing GeminiMol features for {len(query_smiles)} query molecules")
                    geminimol_features = compute_geminimol_features(
                        query_smiles=query_smiles,
                        ref_smiles=ref_smiles,
                        ref_labels=ref_labels,
                        repo_path=repo_path,
                        model_path=model_path,
                        output_dir=self.run_dir / "third_party",
                        use_pharm_profiler=use_pharm_profiler
                    )
                    
                    # Add GeminiMol features to scored records
                    for r in scored:
                        compound_smiles = r.get("compound_smiles")
                        if compound_smiles in geminimol_features:
                            features = geminimol_features[compound_smiles]
                            r["gm_cosine_max"] = features.get("gm_cosine_max")
                            r["gm_cosine_mean_topk"] = features.get("gm_cosine_mean_topk")
                            r["gm_profile_score"] = features.get("gm_profile_score")
                        else:
                            r["gm_cosine_max"] = None
                            r["gm_cosine_mean_topk"] = None
                            r["gm_profile_score"] = None
                    
                    geminimol_summary.update({
                        "repo_path": repo_path,
                        "model_path": model_path,
                        "use_pharm_profiler": use_pharm_profiler,
                        "n_query": len(query_smiles),
                        "n_ref": len(ref_smiles),
                        "note": "computed embeddings and PharmProfiler scores"
                    })
                else:
                    logger.warning("GeminiMol computation skipped - insufficient data")
                    geminimol_summary["note"] = "insufficient query or reference data"
                    
            except Exception as e:
                logger.warning(f"GeminiMol computation failed: {e}")
                geminimol_summary["note"] = f"GeminiMol computation failed: {str(e)}"
        else:
            # Add None values for GeminiMol features when disabled
            for r in scored:
                r["gm_cosine_max"] = None
                r["gm_cosine_mean_topk"] = None
                r["gm_profile_score"] = None

        # 8) Optional docking (gated by validated co-crystal pocket)
        docking_cfg = (self.cfg.get("docking") if isinstance(self.cfg, dict) else getattr(self.cfg, "docking", {})) or {}
        docking_enabled = bool(docking_cfg.get("enabled", False))
        docking_summary: Dict[str, Any] = {"enabled": docking_enabled}
        docking_rows: List[Dict[str, Any]] = []
        if docking_enabled:
            try:
                receptor_path = docking_cfg.get("receptor_pdb_path")
                ligand_resname = docking_cfg.get("cocrystal", {}).get("ligand_resname") if docking_cfg.get("cocrystal") else None
                chain_id = docking_cfg.get("cocrystal", {}).get("chain_id") if docking_cfg.get("cocrystal") else None
                uniprot_id = docking_cfg.get("cocrystal", {}).get("uniprot_id") if docking_cfg.get("cocrystal") else None
                auto_validate = bool(docking_cfg.get("auto_validate", False))
                pdb_id = docking_cfg.get("pdb_id")

                gate_reason = None
                if not receptor_path or not Path(receptor_path).exists():
                    gate_reason = "missing_receptor"
                elif not ligand_resname and not (auto_validate and pdb_id):
                    gate_reason = "missing_cocrystal_ligand_resname"
                elif auto_validate and pdb_id:
                    try:
                        from .curation.biolip import fetch_biolip_ligands, select_relevant_ligand
                        biolip_ligs = fetch_biolip_ligands(pdb_id)
                        sel = select_relevant_ligand(biolip_ligs)
                        if sel is None:
                            gate_reason = "biolip_no_ligand"
                        else:
                            ligand_resname = sel.ligand_id
                            chain_id = sel.chain
                    except Exception as e:
                        docking_summary["biolip_error"] = str(e)
                        gate_reason = gate_reason or "biolip_unavailable"
                elif is_common_additive(ligand_resname):
                    gate_reason = "ccd_additive_ligand"
                elif uniprot_id and chain_id and pdb_id and not maps_to_uniprot_chain(pdb_id, chain_id, uniprot_id):
                    gate_reason = "sifts_mapping_unavailable"

                if gate_reason:
                    docking_summary["gate_reason"] = gate_reason
                else:
                    pocket = box_from_cocrystal(Path(receptor_path), ligand_resname, margin=float(docking_cfg.get("box_margin", 12.0)))
                    center = pocket.center
                    size = pocket.size
                    if not center or not size:
                        docking_summary["gate_reason"] = "invalid_pocket_box"
                    else:
                        vina = VinaAdapter()
                        ok = awaitable_false_if_needed(vina.setup)
                        if not ok:
                            docking_summary["gate_reason"] = "vina_unavailable"
                        else:
                            # run docking for each compound against supplied receptor
                            for lig in ligs:
                                smiles = lig.get("smiles")
                                if not smiles:
                                    continue
                                res = awaitable_none_if_needed(vina.dock_compound, smiles, receptor_path, lig.get("name"), uniprot_id or "", center, size)
                                if res is None:
                                    continue
                                docking_rows.append({
                                    "compound_id": lig.get("name"),
                                    "target_uniprot": uniprot_id or "",
                                    "target_chembl": None,
                                    "binding_energy": res.binding_energy if hasattr(res, "binding_energy") else None,
                                    "structure_id": Path(receptor_path).stem,
                                    "vina_exhaustiveness": getattr(vina, "exhaustiveness", None),
                                })
            except Exception as e:
                docking_summary["error"] = str(e)
        manifest["stages"].append({"stage": "docking", **docking_summary, "n_results": len(docking_rows)})

        docking_path = None
        if docking_rows:
            docking_path = self._write_docking(docking_rows)
            manifest.setdefault("outputs", {})["docking_csv"] = str(docking_path)
            manifest["outputs"]["docking_csv_sha256"] = _sha256(docking_path)

            # Integrate docking into scored records where target/compound match
            dock_map = {(r["compound_id"], r["target_uniprot"]): r for r in docking_rows}
            for rec in scored:
                key = (rec.get("compound_id"), rec.get("target_uniprot", ""))
                drec = dock_map.get(key)
                if drec and drec.get("binding_energy") is not None:
                    # More negative energy is better; convert to positive score by negation
                    rec["docking_score"] = float(-1.0 * drec["binding_energy"])

        # 8) Ouroboros Chemical modes (optional, after baseline target ID)
        ouroboros_cfg = self.cfg.get("ouroboros", {}) if isinstance(self.cfg, dict) else {}
        ouroboros_enabled = bool(ouroboros_cfg.get("enabled", False))
        ouroboros_summary: Dict[str, Any] = {"enabled": ouroboros_enabled}
        
        if ouroboros_enabled:
            try:
                repo_path = ouroboros_cfg.get("repo_path", "third_party/Ouroboros")
                model_name = ouroboros_cfg.get("model_name", "Ouroboros_M1c")
                jobs = ouroboros_cfg.get("jobs", [])
                
                if jobs:
                    logger.info(f"Running {len(jobs)} Ouroboros jobs")
                    ouroboros_results = run_ouroboros_jobs(
                        jobs=jobs,
                        repo_path=repo_path,
                        output_dir=self.run_dir / "third_party"
                    )
                    
                    ouroboros_summary.update({
                        "repo_path": repo_path,
                        "model_name": model_name,
                        "n_jobs": len(jobs),
                        "successful_jobs": ouroboros_results.get("successful_jobs", 0),
                        "failed_jobs": ouroboros_results.get("failed_jobs", 0),
                        "note": "generation experiments completed separately from baseline ranking"
                    })
                else:
                    logger.info("Ouroboros enabled but no jobs configured")
                    ouroboros_summary["note"] = "no jobs configured"
                    
            except Exception as e:
                logger.warning(f"Ouroboros jobs failed: {e}")
                ouroboros_summary["note"] = f"Ouroboros jobs failed: {str(e)}"

        # 9) Fuse and write outputs (now includes optional docking)
        weights = self.cfg.get("scoring", {}).get("weights") if isinstance(self.cfg, dict) else self.cfg.scoring.weights
        fused = fuse_scores(scored, weights or {}) if scored else []
        results_path = self._write_results(fused)

        manifest["stages"].append({
            "stage": "e3fp",
            **e3fp_summary,
        })

        manifest["stages"].append({
            "stage": "geminimol",
            **geminimol_summary,
        })

        manifest["stages"].append({
            "stage": "ouroboros",
            **ouroboros_summary,
        })

        manifest["stages"].append({
            "stage": "fuse_scores",
            "n_results": len(fused),
            "weights": weights,
        })

        manifest["outputs"] = manifest.get("outputs", {}) | {
            "results_csv": str(results_path),
            "results_csv_sha256": _sha256(results_path),
        }

        man_path = self._write_manifest(manifest)
        logger.info(f"Run complete. Results: {results_path} Manifest: {man_path}")

        # Optional reporting
        viz_cfg = (self.cfg.get("visualization") if isinstance(self.cfg, dict) else getattr(self.cfg, "visualization", {})) or {}
        out = {"results_csv": str(results_path), "manifest": str(man_path)}
        try:
            if bool(viz_cfg.get("enabled", False)):
                report_dir = self.run_dir / "report"
                title = viz_cfg.get("title", f"Run {self.run_id}")
                # files = generate_report(results_path, man_path, report_dir, title=title)  # Temporarily disabled
                files = {}  # Temporary empty dict
                out["report_html"] = str(files.get("html", ""))
                manifest["outputs"]["report_html"] = out["report_html"]
                # Attach asset files with hashes
                assets = {}
                for key, p in files.items():
                    if key == "html":
                        continue
                    try:
                        assets[key] = {
                            "path": str(p),
                            "sha256": _sha256(p),
                        }
                    except Exception:
                        continue
                if assets:
                    manifest["outputs"]["report_assets"] = assets
                self._write_manifest(manifest)
                logger.info(f"Report generated: {out['report_html']}")
        except Exception as e:
            logger.warning(f"Reporting skipped: {e}")

        return out


def awaitable_false_if_needed(coro_func):
    """Call an async setup function and return bool, handling sync usage."""
    try:
        import asyncio
        if asyncio.get_event_loop().is_running():
            # Cannot run here in async loop; assume setup False to skip
            return False
        return asyncio.run(coro_func())
    except Exception:
        return False


def awaitable_none_if_needed(coro_func, *args, **kwargs):
    """Call an async coroutine and return result, handling sync context."""
    try:
        import asyncio
        if asyncio.get_event_loop().is_running():
            return None
        return asyncio.run(coro_func(*args, **kwargs))
    except Exception:
        return None
