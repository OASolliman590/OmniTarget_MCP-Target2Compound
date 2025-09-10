"""
AutoDock Vina adapter for molecular docking.

This adapter does not fabricate outputs. It requires a real Vina binary and
valid PDBQT files. If prerequisites are missing, methods return None or raise
with clear messages so callers can skip docking gracefully.
"""

import os
import tempfile
import subprocess
import asyncio
from typing import List, Dict, Any, Optional, Tuple
from pathlib import Path

from loguru import logger

from ..settings import settings
from ..schemas.scoring import DockingScore
from ..io.pdbqt import ligand_smiles_to_pdbqt, receptor_to_pdbqt, require_obabel


class VinaAdapter:
    """Adapter for AutoDock Vina molecular docking."""
    
    def __init__(self):
        self.vina_binary = settings.models.vina_binary_path
        self.exhaustiveness = settings.models.vina_exhaustiveness
        self.num_modes = settings.models.vina_num_modes
        self.temp_dir = Path(settings.files.temp_dir)
        
        logger.info("Initialized Vina adapter")
    
    async def setup(self) -> bool:
        """Set up Vina environment."""
        try:
            # Check if Vina binary exists (either as absolute path or in PATH)
            vina_path = Path(self.vina_binary)
            if not vina_path.exists():
                # Try to find it in PATH
                import shutil
                vina_in_path = shutil.which(self.vina_binary)
                if vina_in_path:
                    self.vina_binary = vina_in_path
                    logger.info(f"Found Vina binary in PATH: {vina_in_path}")
                else:
                    logger.error(f"Vina binary not found: {self.vina_binary}")
                    return False
            
            # Create temp directory
            self.temp_dir.mkdir(exist_ok=True, parents=True)
            
            logger.info("Vina setup completed")
            return True
        except Exception as e:
            logger.error(f"Error setting up Vina: {e}")
            return False
    
    async def dock_compound(
        self,
        smiles: str,
        receptor_pdb_path: str,
        compound_id: str = None,
        target_id: str = None,
        center: Tuple[float, float, float] = None,
        size: Tuple[float, float, float] = (20, 20, 20)
    ) -> Optional[DockingScore]:
        """Dock compound to receptor structure.
        
        Args:
            smiles: Compound SMILES
            receptor_pdb_path: Path to receptor PDB file
            compound_id: Compound identifier
            target_id: Target identifier
            center: Docking box center coordinates
            size: Docking box size
            
        Returns:
            DockingScore object or None if docking failed
        """
        try:
            # Generate 3D coordinates for ligand
            ligand_pdbqt = await self._prepare_ligand(smiles)
            if not ligand_pdbqt:
                return None
            
            # Prepare receptor
            receptor_pdbqt = await self._prepare_receptor(receptor_pdb_path)
            if not receptor_pdbqt:
                return None
            
            # Determine docking box
            if center is None:
                center = await self._calculate_binding_site_center(receptor_pdb_path)
            
            # Run Vina docking
            result = await self._run_vina_docking(
                ligand_pdbqt, receptor_pdbqt, center, size
            )
            
            if result:
                score = DockingScore(
                    compound_id=compound_id or f"compound_{hash(smiles) % 10000}",
                    target_id=target_id or f"target_{hash(receptor_pdb_path) % 10000}",
                    binding_energy=result["binding_energy"],
                    structure_id=Path(receptor_pdb_path).stem,
                    vina_score=result["binding_energy"],
                    vina_exhaustiveness=self.exhaustiveness
                )
                return score
            
            return None
            
        except Exception as e:
            logger.error(f"Error docking compound: {e}")
            return None
    
    async def _prepare_ligand(self, smiles: str) -> Optional[str]:
        """Prepare ligand PDBQT from SMILES using RDKit + Open Babel."""
        try:
            require_obabel()
            with tempfile.NamedTemporaryFile(suffix=".pdbqt", delete=False) as f:
                out = Path(f.name)
            pdbqt = ligand_smiles_to_pdbqt(smiles, out)
            return str(pdbqt)
        except Exception as e:
            logger.error(f"Error preparing ligand PDBQT: {e}")
            return None
    
    async def _prepare_receptor(self, pdb_path: str) -> Optional[str]:
        """Prepare receptor PDBQT from PDB using Open Babel."""
        try:
            require_obabel()
            pdbqt = receptor_to_pdbqt(Path(pdb_path), Path(pdb_path).with_suffix(".pdbqt"))
            return str(pdbqt)
        except Exception as e:
            logger.error(f"Error preparing receptor PDBQT: {e}")
            return None
    
    async def _calculate_binding_site_center(self, pdb_path: str) -> Tuple[float, float, float]:
        """Binding site center is provided by pocket validation; no fallback."""
        raise RuntimeError("Docking box center must be provided from validated pocket")
    
    async def _run_vina_docking(
        self,
        ligand_pdbqt: str,
        receptor_pdbqt: str,
        center: Tuple[float, float, float],
        size: Tuple[float, float, float]
    ) -> Optional[Dict[str, Any]]:
        """Run Vina docking calculation."""
        try:
            with tempfile.NamedTemporaryFile(suffix=".pdbqt", delete=False) as f:
                output_path = f.name
            
            # Build Vina command
            cmd = [
                self.vina_binary,
                "--receptor", receptor_pdbqt,
                "--ligand", ligand_pdbqt,
                "--out", output_path,
                "--center_x", str(center[0]),
                "--center_y", str(center[1]),
                "--center_z", str(center[2]),
                "--size_x", str(size[0]),
                "--size_y", str(size[1]),
                "--size_z", str(size[2]),
                "--exhaustiveness", str(self.exhaustiveness),
                "--num_modes", str(self.num_modes)
            ]
            
            # Run Vina
            process = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE
            )
            
            stdout, stderr = await process.communicate()
            
            if process.returncode == 0:
                out_txt = stdout.decode()
                be = _parse_vina_best_energy(out_txt)
                if be is None:
                    raise RuntimeError("Could not parse Vina energy from output")
                return {
                    "binding_energy": be,
                    "output_path": output_path,
                    "stdout": out_txt,
                    "stderr": stderr.decode(),
                }
            else:
                logger.error(f"Vina failed: {stderr.decode()}")
                return None
                
        except Exception as e:
            logger.error(f"Error running Vina: {e}")
            return None


def _parse_vina_best_energy(stdout_text: str) -> Optional[float]:
    """Parse the best binding energy from Vina stdout text.

    Looks for the docking results table and returns the first energy value.
    Returns None if not found.
    """
    lines = stdout_text.splitlines()
    # Heuristic: lines with mode index at col 1 and energy in col 2
    for line in lines:
        s = line.strip()
        if not s or not s[0].isdigit():
            continue
        parts = s.split()
        if len(parts) >= 2:
            try:
                return float(parts[1])
            except Exception:
                continue
    return None
