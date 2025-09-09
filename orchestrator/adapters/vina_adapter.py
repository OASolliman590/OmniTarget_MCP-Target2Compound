"""
AutoDock Vina adapter for molecular docking.
"""

import os
import tempfile
import subprocess
import asyncio
from typing import List, Dict, Any, Optional, Tuple
from pathlib import Path

from loguru import logger
from rdkit import Chem
from rdkit.Chem import AllChem

from ..settings import settings
from ..schemas.scoring import DockingScore


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
        """Prepare ligand PDBQT from SMILES."""
        try:
            # Convert SMILES to 3D structure
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            
            mol = Chem.AddHs(mol)
            
            # Generate 3D coordinates
            if AllChem.EmbedMolecule(mol) != 0:
                logger.warning(f"Failed to generate 3D coordinates for {smiles}")
                return None
            
            AllChem.MMFFOptimizeMolecule(mol)
            
            # Save as SDF temporarily
            with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as f:
                sdf_path = f.name
            
            writer = Chem.SDWriter(sdf_path)
            writer.write(mol)
            writer.close()
            
            # Convert to PDBQT (would need openbabel or similar)
            # For now, return placeholder path
            pdbqt_path = sdf_path.replace(".sdf", ".pdbqt")
            
            # Placeholder conversion (would use obabel in real implementation)
            with open(pdbqt_path, "w") as f:
                f.write("# Placeholder PDBQT file\n")
            
            return pdbqt_path
            
        except Exception as e:
            logger.error(f"Error preparing ligand: {e}")
            return None
    
    async def _prepare_receptor(self, pdb_path: str) -> Optional[str]:
        """Prepare receptor PDBQT from PDB."""
        try:
            # Convert PDB to PDBQT (would need AutoDockTools or similar)
            pdbqt_path = pdb_path.replace(".pdb", ".pdbqt")
            
            # Placeholder conversion
            with open(pdbqt_path, "w") as f:
                f.write("# Placeholder receptor PDBQT file\n")
            
            return pdbqt_path
            
        except Exception as e:
            logger.error(f"Error preparing receptor: {e}")
            return None
    
    async def _calculate_binding_site_center(self, pdb_path: str) -> Tuple[float, float, float]:
        """Calculate binding site center from PDB."""
        # Placeholder implementation - would analyze PDB structure
        return (0.0, 0.0, 0.0)
    
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
                # Parse results (simplified)
                binding_energy = -7.5  # Placeholder
                
                return {
                    "binding_energy": binding_energy,
                    "output_path": output_path,
                    "stdout": stdout.decode(),
                    "stderr": stderr.decode()
                }
            else:
                logger.error(f"Vina failed: {stderr.decode()}")
                return None
                
        except Exception as e:
            logger.error(f"Error running Vina: {e}")
            return None
