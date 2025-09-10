"""
Configuration schemas for pipeline runs.
"""

from typing import List, Optional, Dict, Any
from pydantic import BaseModel, Field, validator
from pathlib import Path


class PPIConfig(BaseModel):
    """Protein-protein interaction configuration."""
    
    min_confidence: float = Field(
        default=0.7, 
        ge=0.0, 
        le=1.0,
        description="Minimum STRING confidence score (0.0-1.0)"
    )
    max_partners: int = Field(
        default=200,
        ge=1,
        description="Maximum number of protein partners to consider"
    )
    include_experimental: bool = Field(
        default=True,
        description="Include experimentally determined interactions"
    )
    include_predicted: bool = Field(
        default=True,
        description="Include computationally predicted interactions"
    )


class ExpressionConfig(BaseModel):
    """Gene/protein expression configuration."""
    
    tissue_whitelist: List[str] = Field(
        default_factory=list,
        description="Only include targets expressed in these tissues (empty = no filter)"
    )
    tissue_blacklist: List[str] = Field(
        default_factory=list,
        description="Exclude targets expressed in these tissues"
    )
    reliability_min: str = Field(
        default="Supported",
        description="Minimum reliability level from Human Protein Atlas"
    )
    expression_level_min: Optional[float] = Field(
        default=None,
        ge=0.0,
        description="Minimum expression level (if available)"
    )

    @validator('reliability_min')
    def validate_reliability(cls, v):
        valid_levels = ["Approved", "Supported", "Uncertain", "Not supported"]
        if v not in valid_levels:
            raise ValueError(f"reliability_min must be one of {valid_levels}")
        return v


class StructureConfig(BaseModel):
    """Protein structure configuration."""
    
    prefer_experimental: bool = Field(
        default=True,
        description="Prefer experimental structures from PDB over predicted"
    )
    allow_alphafold: bool = Field(
        default=True,
        description="Allow AlphaFold predicted structures when PDB unavailable"
    )
    min_resolution: Optional[float] = Field(
        default=None,
        ge=0.0,
        description="Minimum resolution for experimental structures (Angstroms)"
    )
    max_resolution: Optional[float] = Field(
        default=3.0,
        ge=0.0,
        description="Maximum resolution for experimental structures (Angstroms)"
    )
    require_structure: bool = Field(
        default=False,
        description="Skip targets without available structures"
    )


class CompoundConfig(BaseModel):
    """Compound input configuration."""
    
    input_paths: List[str] = Field(
        description="Paths to compound files (.smi, .sdf, .mol)"
    )
    max_compounds: Optional[int] = Field(
        default=None,
        ge=1,
        description="Maximum number of compounds to process"
    )
    filter_duplicates: bool = Field(
        default=True,
        description="Remove duplicate compounds based on canonical SMILES"
    )
    min_molecular_weight: Optional[float] = Field(
        default=None,
        ge=0.0,
        description="Minimum molecular weight filter"
    )
    max_molecular_weight: Optional[float] = Field(
        default=None,
        ge=0.0,
        description="Maximum molecular weight filter"
    )

    @validator('input_paths')
    def validate_paths(cls, v):
        for path in v:
            if not Path(path).exists():
                raise ValueError(f"Compound file does not exist: {path}")
            if not Path(path).suffix.lower() in ['.smi', '.sdf', '.mol']:
                raise ValueError(f"Unsupported file format: {path}")
        return v


class ScoringConfig(BaseModel):
    """Scoring and ranking configuration."""
    
    weights: Dict[str, float] = Field(
        default={"similarity": 0.5, "pharmacophore": 0.2, "docking": 0.1, "evidence": 0.2},
        description="Weights for different scoring components"
    )
    normalize_scores: bool = Field(
        default=True,
        description="Z-score normalize individual scores before combining"
    )
    top_n_results: int = Field(
        default=100,
        ge=1,
        description="Number of top results to return"
    )
    min_combined_score: Optional[float] = Field(
        default=None,
        description="Minimum combined score threshold"
    )

    @validator('weights')
    def validate_weights(cls, v):
        total = sum(v.values())
        if abs(total - 1.0) > 0.01:
            raise ValueError(f"Weights must sum to 1.0, got {total}")
        
        # Basic sanity: allow flexible keys; ensure non-negative
        
        for key, weight in v.items():
            if weight < 0:
                raise ValueError(f"Weight for {key} must be non-negative")
        
        return v


class RunConfig(BaseModel):
    """Main configuration for a pipeline run."""
    
    # Basic parameters
    disease_terms: List[str] = Field(
        description="Disease terms for target discovery"
    )
    organism: str = Field(
        default="Homo sapiens",
        description="Target organism"
    )
    
    # Component configurations
    ppi: PPIConfig = Field(default_factory=PPIConfig)
    expression: ExpressionConfig = Field(default_factory=ExpressionConfig)
    structures: StructureConfig = Field(default_factory=StructureConfig)
    compounds: CompoundConfig = Field(description="Compound input configuration")
    scoring: ScoringConfig = Field(default_factory=ScoringConfig)
    
    # Output settings
    output_dir: str = Field(
        default="./data/outputs",
        description="Directory for output files"
    )
    run_name: Optional[str] = Field(
        default=None,
        description="Custom name for this run (auto-generated if not provided)"
    )
    
    # Advanced settings
    max_targets: Optional[int] = Field(
        default=None,
        ge=1,
        description="Maximum number of targets to process"
    )
    parallel_jobs: int = Field(
        default=4,
        ge=1,
        description="Number of parallel jobs for CPU-intensive tasks"
    )
    cache_results: bool = Field(
        default=True,
        description="Cache intermediate results"
    )
    
    # Debugging and development
    debug_mode: bool = Field(
        default=False,
        description="Enable debug logging and outputs"
    )
    dry_run: bool = Field(
        default=False,
        description="Validate configuration without running pipeline"
    )

    @validator('disease_terms')
    def validate_disease_terms(cls, v):
        if not v:
            raise ValueError("At least one disease term must be provided")
        return v

    @validator('organism')
    def validate_organism(cls, v):
        # Common organisms - could be expanded
        valid_organisms = [
            "Homo sapiens", "Mus musculus", "Rattus norvegicus",
            "Drosophila melanogaster", "Caenorhabditis elegans",
            "Saccharomyces cerevisiae", "Escherichia coli"
        ]
        if v not in valid_organisms:
            # Allow other organisms but warn
            pass
        return v

    class Config:
        extra = "forbid"  # Don't allow extra fields
        validate_assignment = True
