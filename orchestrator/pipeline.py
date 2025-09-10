"""
Main pipeline orchestration for drug discovery.
"""

import asyncio
import hashlib
import json
import time
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Any, Optional

import pandas as pd
from loguru import logger

from .settings import settings
from .schemas.config import RunConfig
from .schemas.pipeline import PipelineRun, PipelineStatus, PipelineStage, StageResult, PipelineResult
from .schemas.targets import Target
from .schemas.compounds import Compound
from .schemas.scoring import RankedHit, DockingScore, EvidenceRecord

from .mcp_clients import (
    KEGGClient, ReactomeClient, ProteinAtlasClient, 
    STRINGClient, UniProtClient, PDBClient, ChEMBLClient
)
from .adapters import GeminiMolAdapter, VinaAdapter, OuroborosAdapter
from .scoring import ScoreNormalizer, ScoreFusion


class DrugDiscoveryPipeline:
    """Main drug discovery pipeline orchestrator."""
    
    def __init__(self):
        """Initialize pipeline with all components."""
        # MCP clients
        self.kegg_client = KEGGClient()
        self.reactome_client = ReactomeClient()
        self.proteinatlas_client = ProteinAtlasClient()
        self.string_client = STRINGClient()
        self.uniprot_client = UniProtClient()
        self.pdb_client = PDBClient()
        self.chembl_client = ChEMBLClient()
        
        # ML adapters
        self.geminimol_adapter = GeminiMolAdapter()
        self.vina_adapter = VinaAdapter()
        self.ouroboros_adapter = OuroborosAdapter()
        
        # Scoring
        self.score_normalizer = ScoreNormalizer()
        
        logger.info("Initialized drug discovery pipeline")
    
    async def setup(self) -> bool:
        """Set up all pipeline components."""
        try:
            logger.info("Setting up pipeline components...")
            
            # Setup adapters
            adapters = [
                self.geminimol_adapter,
                self.vina_adapter,
                self.ouroboros_adapter
            ]
            
            setup_results = await asyncio.gather(
                *[adapter.setup() for adapter in adapters],
                return_exceptions=True
            )
            
            for i, result in enumerate(setup_results):
                if isinstance(result, Exception):
                    logger.error(f"Failed to setup adapter {i}: {result}")
                elif not result:
                    logger.warning(f"Adapter {i} setup returned False")
            
            logger.info("Pipeline setup completed")
            return True
            
        except Exception as e:
            logger.error(f"Error setting up pipeline: {e}")
            return False
    
    async def run(self, config: RunConfig) -> PipelineResult:
        """Execute the complete drug discovery pipeline.
        
        Args:
            config: Pipeline configuration
            
        Returns:
            Pipeline execution results
        """
        # Create run record
        run_id = self._generate_run_id(config)
        run = PipelineRun(
            run_id=run_id,
            run_name=config.run_name or run_id,
            config_hash=self._hash_config(config),
            config_data=config.dict(),
            status=PipelineStatus.RUNNING,
            started_at=datetime.now(),
            current_stage=PipelineStage.INITIALIZATION
        )
        
        try:
            logger.info(f"Starting pipeline run: {run_id}")
            
            # Stage 1: Target Discovery
            targets = await self._discover_targets(config, run)
            
            # Stage 2: Target Characterization
            targets = await self._characterize_targets(targets, config, run)
            
            # Stage 3: Compound Loading
            compounds = await self._load_compounds(config, run)
            
            # Stage 4: Feature Generation
            compounds = await self._generate_features(compounds, run)
            
            # Stage 5: Molecular Docking
            docking_scores = await self._dock_compounds(compounds, targets, run)
            
            # Stage 6: Evidence Validation
            evidence_records = await self._validate_evidence(compounds, targets, run)
            
            # Stage 7: Score Integration (sequence-based predictor removed)
            ranked_hits = await self._integrate_scores(docking_scores, evidence_records, config, run)
            
            # Stage 8: Result Generation
            result = await self._generate_results(ranked_hits, config, run)
            
            # Update run status
            run.status = PipelineStatus.COMPLETED
            run.completed_at = datetime.now()
            run.duration = (run.completed_at - run.started_at).total_seconds()
            
            logger.info(f"Pipeline run completed: {run_id}")
            return result
            
        except Exception as e:
            logger.error(f"Pipeline run failed: {e}")
            run.status = PipelineStatus.FAILED
            run.error_message = str(e)
            run.completed_at = datetime.now()
            raise
    
    async def _discover_targets(self, config: RunConfig, run: PipelineRun) -> List[Target]:
        """Stage 1: Discover targets from disease terms."""
        stage = PipelineStage.TARGET_DISCOVERY
        stage_result = StageResult(stage=stage, status=PipelineStatus.RUNNING, started_at=datetime.now())
        
        try:
            logger.info(f"Discovering targets for diseases: {config.disease_terms}")
            
            targets = []
            
            # Get pathways from KEGG
            kegg_pathways = await self.kegg_client.search_disease_pathways(config.disease_terms)
            
            # Get enriched pathways from Reactome
            reactome_pathways = []
            for term in config.disease_terms:
                pathways = await self.reactome_client.search_pathways(term)
                reactome_pathways.extend(pathways)
            
            # Extract genes from pathways
            all_genes = set()
            for pathway in kegg_pathways + reactome_pathways:
                if pathway.get("genes"):
                    all_genes.update(pathway["genes"])
            
            # Create target objects
            for gene in all_genes:
                target = Target(
                    target_id=gene,
                    gene_name=gene,
                    discovery_source=["KEGG", "Reactome"],
                    discovery_pathways=kegg_pathways + reactome_pathways
                )
                targets.append(target)
            
            # Limit targets if specified
            if config.max_targets:
                targets = targets[:config.max_targets]
            
            stage_result.status = PipelineStatus.COMPLETED
            stage_result.completed_at = datetime.now()
            stage_result.items_processed = len(targets)
            stage_result.items_successful = len(targets)
            
            logger.info(f"Discovered {len(targets)} targets")
            return targets
            
        except Exception as e:
            stage_result.status = PipelineStatus.FAILED
            stage_result.error_message = str(e)
            logger.error(f"Target discovery failed: {e}")
            raise
        finally:
            stage_result.duration = (datetime.now() - stage_result.started_at).total_seconds()
            run.add_stage_result(stage_result)
    
    async def _characterize_targets(self, targets: List[Target], config: RunConfig, run: PipelineRun) -> List[Target]:
        """Stage 2: Characterize targets with additional data."""
        stage = PipelineStage.TARGET_CHARACTERIZATION
        stage_result = StageResult(stage=stage, status=PipelineStatus.RUNNING, started_at=datetime.now())
        
        try:
            logger.info(f"Characterizing {len(targets)} targets")
            
            characterized_targets = []
            
            for target in targets:
                try:
                    # Get protein sequence from UniProt
                    protein_info = await self.uniprot_client.search_proteins(target.gene_name)
                    if protein_info:
                        target.sequence_data = protein_info[0]
                        target.uniprot_id = protein_info[0].get("accession")
                    
                    # Get expression data from Protein Atlas
                    expression_data = await self.proteinatlas_client.get_expression_data(target.gene_name)
                    target.expression_data = expression_data
                    
                    # Get PPI data from STRING
                    if target.uniprot_id:
                        ppi_data = await self.string_client.get_protein_interactions([target.uniprot_id])
                        target.ppi_data = ppi_data
                    
                    # Get structure data from PDB
                    if target.uniprot_id:
                        structures = await self.pdb_client.search_structures(target.uniprot_id)
                        target.structure_data = structures
                    
                    # Update quality flags
                    target.update_flags()
                    characterized_targets.append(target)
                    
                except Exception as e:
                    logger.warning(f"Failed to characterize target {target.target_id}: {e}")
                    characterized_targets.append(target)  # Keep even if characterization failed
            
            stage_result.status = PipelineStatus.COMPLETED
            stage_result.completed_at = datetime.now()
            stage_result.items_processed = len(targets)
            stage_result.items_successful = len(characterized_targets)
            
            logger.info(f"Characterized {len(characterized_targets)} targets")
            return characterized_targets
            
        except Exception as e:
            stage_result.status = PipelineStatus.FAILED
            stage_result.error_message = str(e)
            logger.error(f"Target characterization failed: {e}")
            raise
        finally:
            stage_result.duration = (datetime.now() - stage_result.started_at).total_seconds()
            run.add_stage_result(stage_result)
    
    async def _load_compounds(self, config: RunConfig, run: PipelineRun) -> List[Compound]:
        """Stage 3: Load and process compounds."""
        stage = PipelineStage.COMPOUND_LOADING
        stage_result = StageResult(stage=stage, status=PipelineStatus.RUNNING, started_at=datetime.now())
        
        try:
            logger.info(f"Loading compounds from {len(config.compounds.input_paths)} files")
            
            compounds = []
            compound_id_counter = 1
            
            for file_path in config.compounds.input_paths:
                try:
                    file_path_obj = Path(file_path)
                    
                    if file_path_obj.suffix.lower() == ".smi":
                        # Read SMILES file
                        with open(file_path, "r") as f:
                            for line_num, line in enumerate(f, 1):
                                line = line.strip()
                                if line and not line.startswith("#"):
                                    parts = line.split()
                                    smiles = parts[0]
                                    name = parts[1] if len(parts) > 1 else f"compound_{compound_id_counter}"
                                    
                                    compound = Compound(
                                        compound_id=f"comp_{compound_id_counter:06d}",
                                        smiles=smiles,
                                        name=name,
                                        source_file=str(file_path),
                                        source_line=line_num
                                    )
                                    compounds.append(compound)
                                    compound_id_counter += 1
                    
                    # Add support for SDF/MOL files here if needed
                    
                except Exception as e:
                    logger.error(f"Error reading file {file_path}: {e}")
                    continue
            
            # Apply filters
            if config.compounds.filter_duplicates:
                compounds = self._deduplicate_compounds(compounds)
            
            if config.compounds.max_compounds:
                compounds = compounds[:config.compounds.max_compounds]
            
            stage_result.status = PipelineStatus.COMPLETED
            stage_result.completed_at = datetime.now()
            stage_result.items_processed = len(compounds)
            stage_result.items_successful = len(compounds)
            
            logger.info(f"Loaded {len(compounds)} compounds")
            return compounds
            
        except Exception as e:
            stage_result.status = PipelineStatus.FAILED
            stage_result.error_message = str(e)
            logger.error(f"Compound loading failed: {e}")
            raise
        finally:
            stage_result.duration = (datetime.now() - stage_result.started_at).total_seconds()
            run.add_stage_result(stage_result)
    
    def _deduplicate_compounds(self, compounds: List[Compound]) -> List[Compound]:
        """Remove duplicate compounds based on canonical SMILES."""
        seen_smiles = set()
        unique_compounds = []
        
        for compound in compounds:
            if compound.smiles not in seen_smiles:
                seen_smiles.add(compound.smiles)
                unique_compounds.append(compound)
        
        logger.info(f"Removed {len(compounds) - len(unique_compounds)} duplicate compounds")
        return unique_compounds
    
    async def _generate_features(self, compounds: List[Compound], run: PipelineRun) -> List[Compound]:
        """Generate molecular features for compounds."""
        stage_result = StageResult(
            stage=PipelineStage.FEATURE_GENERATION,
            status=PipelineStatus.RUNNING,
            started_at=datetime.now()
        )
        
        try:
            logger.info(f"Generating features for {len(compounds)} compounds")
            
            # Generate GeminiMol embeddings for each compound
            for compound in compounds:
                try:
                    embedding = await self.geminimol_adapter.compute_embedding(compound.smiles)
                    compound.features = {
                        'geminimol_embedding': embedding,
                        'molecular_weight': self._calculate_molecular_weight(compound.smiles),
                        'logp': self._calculate_logp(compound.smiles)
                    }
                except Exception as e:
                    logger.warning(f"Failed to generate features for {compound.compound_id}: {e}")
                    compound.features = {}
            
            stage_result.status = PipelineStatus.COMPLETED
            logger.info(f"Generated features for {len(compounds)} compounds")
            return compounds
            
        except Exception as e:
            stage_result.status = PipelineStatus.FAILED
            stage_result.error_message = str(e)
            logger.error(f"Feature generation failed: {e}")
            raise
        finally:
            stage_result.duration = (datetime.now() - stage_result.started_at).total_seconds()
            run.add_stage_result(stage_result)
    
    def _calculate_molecular_weight(self, smiles: str) -> float:
        """Calculate molecular weight from SMILES."""
        try:
            from rdkit import Chem
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                return Chem.rdMolDescriptors.CalcExactMolWt(mol)
        except:
            pass
        return 0.0
    
    def _calculate_logp(self, smiles: str) -> float:
        """Calculate LogP from SMILES."""
        try:
            from rdkit import Chem
            from rdkit.Chem import Crippen
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                return Crippen.MolLogP(mol)
        except:
            pass
        return 0.0
    
    # Sequence-based predictor scoring removed
    
    async def _dock_compounds(self, compounds: List[Compound], targets: List[Target], run: PipelineRun) -> List[DockingScore]:
        """Dock compounds against targets using Vina."""
        stage_result = StageResult(
            stage=PipelineStage.DOCKING,
            status=PipelineStatus.RUNNING,
            started_at=datetime.now()
        )
        
        try:
            logger.info(f"Docking {len(compounds)} compounds against {len(targets)} targets")
            
            scores = []
            for compound in compounds:
                for target in targets:
                    try:
                        # Use placeholder PDB path for now
                        score = await self.vina_adapter.dock_compound(
                            compound.smiles, 
                            "placeholder.pdb",
                            compound_id=compound.compound_id,
                            target_id=target.target_id
                        )
                        if score:
                            scores.append(score)
                    except Exception as e:
                        logger.warning(f"Docking failed for {compound.compound_id}-{target.target_id}: {e}")
            
            stage_result.status = PipelineStatus.COMPLETED
            logger.info(f"Generated {len(scores)} docking scores")
            return scores
            
        except Exception as e:
            stage_result.status = PipelineStatus.FAILED
            stage_result.error_message = str(e)
            logger.error(f"Molecular docking failed: {e}")
            raise
        finally:
            stage_result.duration = (datetime.now() - stage_result.started_at).total_seconds()
            run.add_stage_result(stage_result)
    
    async def _validate_evidence(self, compounds: List[Compound], targets: List[Target], run: PipelineRun) -> List[EvidenceRecord]:
        """Validate compound-target interactions with evidence databases."""
        stage_result = StageResult(
            stage=PipelineStage.EVIDENCE_VALIDATION,
            status=PipelineStatus.RUNNING,
            started_at=datetime.now()
        )
        
        try:
            logger.info(f"Validating evidence for {len(compounds)} compounds against {len(targets)} targets")
            
            # Placeholder implementation - would query ChEMBL, PDB, etc.
            evidence_records = []
            for compound in compounds:
                for target in targets:
                    evidence = EvidenceRecord(
                        compound_id=compound.compound_id,
                        target_id=target.target_id,
                        evidence_type="placeholder",
                        confidence=0.5,
                        source="placeholder"
                    )
                    evidence_records.append(evidence)
            
            stage_result.status = PipelineStatus.COMPLETED
            logger.info(f"Generated {len(evidence_records)} evidence records")
            return evidence_records
            
        except Exception as e:
            stage_result.status = PipelineStatus.FAILED
            stage_result.error_message = str(e)
            logger.error(f"Evidence validation failed: {e}")
            raise
        finally:
            stage_result.duration = (datetime.now() - stage_result.started_at).total_seconds()
            run.add_stage_result(stage_result)
    
    async def _integrate_scores(self, docking_scores: List[DockingScore], 
                              evidence_records: List[EvidenceRecord], config: RunConfig, run: PipelineRun) -> List[RankedHit]:
        """Integrate all scores and rank hits."""
        stage_result = StageResult(
            stage=PipelineStage.SCORE_INTEGRATION,
            status=PipelineStatus.RUNNING,
            started_at=datetime.now()
        )
        
        try:
            logger.info("Score integration skipped (sequence-based predictor removed in legacy pipeline)")
            stage_result.status = PipelineStatus.COMPLETED
            return []
            
        except Exception as e:
            stage_result.status = PipelineStatus.FAILED
            stage_result.error_message = str(e)
            logger.error(f"Score integration failed: {e}")
            raise
        finally:
            stage_result.duration = (datetime.now() - stage_result.started_at).total_seconds()
            run.add_stage_result(stage_result)
    
    async def _generate_results(self, ranked_hits: List[RankedHit], config: RunConfig, run: PipelineRun) -> PipelineResult:
        """Generate final pipeline results."""
        stage_result = StageResult(
            stage=PipelineStage.RESULT_GENERATION,
            status=PipelineStatus.RUNNING,
            started_at=datetime.now()
        )
        
        try:
            logger.info("Generating final results")
            
            # Create pipeline result
            result = PipelineResult(
                run_id=run.run_id,
                status=run.status,
                compounds=[],  # Would be populated from compounds
                targets=[],    # Would be populated from targets
                results=ranked_hits,
                total_compounds=0,  # Would be populated from actual compounds
                total_targets=0,    # Would be populated from actual targets
                total_pairs=len(ranked_hits),
                top_hits_count=len(ranked_hits),
                total_time=run.duration or 0.0,
                created_at=datetime.now()
            )
            
            stage_result.status = PipelineStatus.COMPLETED
            logger.info("Results generated successfully")
            return result
            
        except Exception as e:
            stage_result.status = PipelineStatus.FAILED
            stage_result.error_message = str(e)
            logger.error(f"Result generation failed: {e}")
            raise
        finally:
            stage_result.duration = (datetime.now() - stage_result.started_at).total_seconds()
            run.add_stage_result(stage_result)
    
    def _generate_run_id(self, config: RunConfig) -> str:
        """Generate unique run identifier."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        config_hash = self._hash_config(config)[:8]
        return f"run_{timestamp}_{config_hash}"
    
    def _hash_config(self, config: RunConfig) -> str:
        """Generate hash of configuration."""
        config_str = json.dumps(config.dict(), sort_keys=True)
        return hashlib.md5(config_str.encode()).hexdigest()
