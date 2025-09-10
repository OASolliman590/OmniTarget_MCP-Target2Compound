"""
Chemical IO utilities for multi-format ligand ingestion and export.

Supports input formats: .smi, .sdf, .mol, .mol2 (best-effort).
Outputs: merged .smi and .sdf (if RDKit available). No synthetic data.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple, TypedDict, Literal

from loguru import logger


def _try_import_rdkit():
    try:
        from rdkit import Chem  # type: ignore
        from rdkit.Chem import AllChem, Descriptors  # type: ignore
        from rdkit.Chem.MolStandardize import charge  # type: ignore
        from rdkit.Chem.rdmolfiles import SDWriter  # type: ignore
        from rdkit.Chem.rdmolfiles import MolToMolBlock  # type: ignore
        return Chem, AllChem, charge, Descriptors, SDWriter, MolToMolBlock
    except Exception:
        return None, None, None, None, None, None


class LigandRecord(TypedDict, total=False):
    source_path: str
    name: str
    smiles: str
    inchikey: str
    mol: object
    ph_used: float
    enumeration_id: Optional[str]
    notes: Dict[str, str]


def read_ligands(paths: Iterable[Path]) -> List[LigandRecord]:
    recs: List[LigandRecord] = []
    Chem, *_ = _try_import_rdkit()
    for p in paths:
        if not p.exists():
            logger.warning(f"Ligand file not found: {p}")
            continue
        suf = p.suffix.lower()
        if suf == ".smi":
            with open(p) as f:
                for i, line in enumerate(f):
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split()
                    smi = parts[0]
                    name = parts[1] if len(parts) > 1 else f"{p.stem}_{i+1}"
                    recs.append({
                        "source_path": str(p),
                        "name": name,
                        "smiles": smi,
                        "inchikey": "",
                        "mol": None,
                        "ph_used": 7.4,
                        "enumeration_id": None,
                        "notes": {},
                    })
        elif suf in {".sdf", ".mol", ".mol2"} and Chem is not None:
            try:
                suppl = Chem.SDMolSupplier(str(p), removeHs=False, sanitize=True) if suf == ".sdf" else None
                if suf in {".mol", ".mol2"}:
                    mol = Chem.MolFromMolFile(str(p), sanitize=True) if suf == ".mol" else Chem.MolFromMol2File(str(p), sanitize=True)
                    suppl = [mol] if mol is not None else []
                for idx, mol in enumerate(suppl or []):
                    if mol is None:
                        continue
                    smi = Chem.MolToSmiles(mol)
                    name = mol.GetProp('_Name') if mol.HasProp('_Name') else f"{p.stem}_{idx+1}"
                    recs.append({
                        "source_path": str(p),
                        "name": name,
                        "smiles": smi,
                        "inchikey": Chem.InchiToInchiKey(Chem.MolToInchi(mol)) if Chem is not None else "",
                        "mol": mol,
                        "ph_used": 7.4,
                        "enumeration_id": None,
                        "notes": {},
                    })
            except Exception as e:
                logger.warning(f"Failed reading {p}: {e}")
        else:
            logger.warning(f"Unsupported ligand format (or RDKit missing) for {p}")
    return recs


def standardize(
    recs: List[LigandRecord],
    ph: float = 7.4,
    strip_salts: bool = True,
    largest_fragment_only: bool = True,
    kekulize: bool = True,
) -> List[LigandRecord]:
    Chem, AllChem, charge, Descriptors, *_ = _try_import_rdkit()
    out: List[LigandRecord] = []
    for r in recs:
        r = dict(r)
        r["ph_used"] = ph
        if Chem is None:
            out.append(r)
            continue
        mol = r.get("mol") or Chem.MolFromSmiles(r.get("smiles", ""))
        if mol is None:
            continue
        try:
            if strip_salts or largest_fragment_only:
                frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
                if frags:
                    frags = sorted(frags, key=lambda m: m.GetNumHeavyAtoms(), reverse=True)
                    mol = frags[0]
            if kekulize:
                Chem.Kekulize(mol, clearAromaticFlags=True)
            # Neutralize/protonate approximately
            if charge is not None:
                uncharger = charge.Uncharger()
                mol = uncharger.uncharge(mol)
            r["mol"] = mol
            r["smiles"] = Chem.MolToSmiles(mol)
            try:
                r["inchikey"] = Chem.InchiToInchiKey(Chem.MolToInchi(mol))
            except Exception:
                r["inchikey"] = ""
            out.append(r)
        except Exception as e:
            logger.warning(f"Standardization failed for {r.get('name')}: {e}")
    return out


def enumerate_chemistry(
    recs: List[LigandRecord],
    max_tautomers: int = 8,
    max_stereoisomers: int = 4,
    allow_chiral_perms: bool = False,
) -> List[LigandRecord]:
    # Keep as a no-op scaffold (do not fabricate enumerations in this environment)
    return recs


def generate_3d(
    recs: List[LigandRecord],
    max_confs: int = 10,
    prune_rms_thresh: float = 0.5,
    ff: Literal["MMFF94", "UFF"] = "MMFF94",
    max_iter: int = 200,
) -> List[LigandRecord]:
    Chem, AllChem, *_ = _try_import_rdkit()
    if Chem is None:
        return recs
    for r in recs:
        mol = r.get("mol")
        if mol is None:
            m = Chem.MolFromSmiles(r.get("smiles", ""))
            if m is None:
                continue
            mol = Chem.AddHs(m)
        try:
            params = AllChem.ETKDGv3()
            params.randomSeed = 0xC0FFEE
            params.pruneRmsThresh = float(prune_rms_thresh)
            AllChem.EmbedMultipleConfs(mol, numConfs=int(max_confs), params=params)
            for cid in range(mol.GetNumConformers()):
                try:
                    if ff == "MMFF94":
                        AllChem.MMFFOptimizeMolecule(mol, confId=cid, maxIters=int(max_iter))
                    else:
                        AllChem.UFFOptimizeMolecule(mol, confId=cid, maxIters=int(max_iter))
                except Exception:
                    pass
            r["mol"] = mol
        except Exception as e:
            logger.warning(f"3D generation failed for {r.get('name')}: {e}")
    return recs


def export(
    recs: List[LigandRecord],
    out_dir: Path,
    formats: Tuple[Literal["smi", "sdf", "mol", "pdbqt"], ...] = ("smi", "sdf"),
) -> Dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    Chem, AllChem, charge, Descriptors, SDWriter, MolToMolBlock = _try_import_rdkit()
    outputs: Dict[str, Path] = {}

    # SMI export
    if "smi" in formats:
        smi_path = out_dir / "ligands.smi"
        with open(smi_path, "w") as f:
            for r in recs:
                smi = r.get("smiles")
                if not smi:
                    continue
                name = r.get("name") or r.get("inchikey") or "ligand"
                f.write(f"{smi}\t{name}\n")
        outputs["smi"] = smi_path

    # SDF export (only if RDKit available)
    if "sdf" in formats and Chem is not None and SDWriter is not None:
        sdf_path = out_dir / "ligands.sdf"
        writer = SDWriter(str(sdf_path))
        for r in recs:
            mol = r.get("mol")
            if mol is None:
                continue
            try:
                mol.SetProp("_Name", r.get("name", "ligand"))
            except Exception:
                pass
            writer.write(mol)
        writer.close()
        outputs["sdf"] = sdf_path

    # PDBQT export intentionally omitted (docking disabled by default)
    return outputs

