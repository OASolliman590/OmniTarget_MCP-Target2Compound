"""
Scoring and ranking schemas.
"""

from typing import List, Optional, Dict, Any
from pydantic import BaseModel, Field, validator
from datetime import datetime
from enum import Enum


class ScoreType(str, Enum):
    """Types of scoring methods."""
    DOCKING = "docking"
    EVIDENCE = "evidence"
    COMBINED = "combined"


"""Sequence-based predictor removed."""


class DockingScore(BaseModel):
    """Molecular docking score and pose information."""
    
    compound_id: str = Field(description="Compound identifier")
    target_id: str = Field(description="Target identifier")
    binding_energy: float = Field(description="Binding energy (kcal/mol)")
    
    # Docking details
    structure_id: str = Field(description="Structure used for docking")
    pose_rank: int = Field(default=1, description="Pose ranking (1 = best)")
    rmsd: Optional[float] = Field(default=None, description="RMSD from reference pose")
    
    # Vina-specific scores
    vina_score: float = Field(description="AutoDock Vina score")
    efficiency_index: Optional[float] = Field(default=None, description="Ligand efficiency index")
    
    # Pose information
    pose_file: Optional[str] = Field(default=None, description="Path to pose file")
    binding_site: Optional[Dict[str, Any]] = Field(default=None, description="Binding site information")
    
    # Processing metadata
    vina_exhaustiveness: Optional[int] = Field(default=None, description="Vina exhaustiveness parameter")
    processing_time: Optional[float] = Field(default=None, description="Docking time in seconds")
    created_at: datetime = Field(default_factory=datetime.now)
    
    class Config:
        extra = "allow"


class EvidenceRecord(BaseModel):
    """Bioactivity evidence from ChEMBL or literature."""
    
    compound_id: str = Field(description="Compound identifier")
    target_id: str = Field(description="Target identifier")
    
    # Activity data
    activity_type: str = Field(description="Type of activity (IC50, Ki, etc.)")
    activity_value: Optional[float] = Field(default=None, description="Activity value")
    activity_units: Optional[str] = Field(default=None, description="Activity units")
    pchembl_value: Optional[float] = Field(default=None, description="pChEMBL value")
    
    # Assay information
    assay_id: Optional[str] = Field(default=None, description="Assay identifier")
    assay_type: Optional[str] = Field(default=None, description="Assay type")
    organism: Optional[str] = Field(default=None, description="Test organism")
    
    # Source information
    source: str = Field(description="Data source (ChEMBL, literature, etc.)")
    reference: Optional[str] = Field(default=None, description="Literature reference")
    confidence: Optional[str] = Field(default=None, description="Confidence level")
    
    # Quality flags
    is_active: Optional[bool] = Field(default=None, description="Compound is active")
    is_high_confidence: bool = Field(default=False, description="High confidence evidence")
    
    created_at: datetime = Field(default_factory=datetime.now)
    
    class Config:
        extra = "allow"


class IntegratedScore(BaseModel):
    """Integrated score combining multiple scoring methods."""
    
    compound_id: str = Field(description="Compound identifier")
    target_id: str = Field(description="Target identifier")
    
    # Individual scores
    docking_score: Optional[float] = Field(default=None, description="Normalized docking score")
    evidence_score: Optional[float] = Field(default=None, description="Evidence-based score")
    
    # Combined score
    combined_score: float = Field(description="Weighted combined score")
    rank: Optional[int] = Field(default=None, description="Overall ranking")
    
    # Scoring metadata
    weights_used: Dict[str, float] = Field(description="Weights used for combination")
    normalization_method: str = Field(default="zscore", description="Normalization method")
    
    # Confidence and quality
    confidence: Optional[float] = Field(default=None, description="Overall confidence")
    num_evidence_sources: int = Field(default=0, description="Number of evidence sources")
    
    created_at: datetime = Field(default_factory=datetime.now)
    
    @validator('weights_used')
    def validate_weights(cls, v):
        """Ensure weights sum to approximately 1.0."""
        total = sum(v.values())
        if abs(total - 1.0) > 0.01:
            raise ValueError(f"Weights must sum to 1.0, got {total}")
        return v
    
    class Config:
        extra = "allow"


class RankedHit(BaseModel):
    """A ranked compound-target pair with all associated data."""
    
    # Identifiers
    compound_id: str = Field(description="Compound identifier")
    target_id: str = Field(description="Target identifier")
    rank: int = Field(description="Overall ranking")
    
    # Scores
    integrated_score: IntegratedScore = Field(description="Integrated scoring data")
    docking_scores: List[DockingScore] = Field(default_factory=list, description="Docking results")
    evidence_records: List[EvidenceRecord] = Field(default_factory=list, description="Bioactivity evidence")
    
    # Summary information
    compound_name: Optional[str] = Field(default=None, description="Compound name")
    compound_smiles: str = Field(description="Compound SMILES")
    target_name: Optional[str] = Field(default=None, description="Target name")
    target_gene: Optional[str] = Field(default=None, description="Target gene symbol")
    
    # Quality indicators
    has_experimental_evidence: bool = Field(default=False, description="Has experimental bioactivity data")
    has_structure: bool = Field(default=False, description="Target has 3D structure")
    confidence_tier: str = Field(default="medium", description="Confidence tier (high/medium/low)")
    
    # Metadata
    created_at: datetime = Field(default_factory=datetime.now)
    
    @validator('confidence_tier')
    def validate_confidence_tier(cls, v):
        """Validate confidence tier values."""
        valid_tiers = ["high", "medium", "low"]
        if v not in valid_tiers:
            raise ValueError(f"confidence_tier must be one of {valid_tiers}")
        return v
    
    class Config:
        extra = "allow"
    
    def update_quality_flags(self) -> None:
        """Update quality flags based on available data."""
        self.has_experimental_evidence = len(self.evidence_records) > 0
        
        # Determine confidence tier based on available data
        score = 0
        # Sequence-based predictor removed
        if self.docking_scores:
            score += 2  # Docking weighted higher
        if self.evidence_records:
            score += 3  # Evidence weighted highest
        
        if score >= 5:
            self.confidence_tier = "high"
        elif score >= 3:
            self.confidence_tier = "medium"
        else:
            self.confidence_tier = "low"


class ScoringReport(BaseModel):
    """Summary report of scoring results."""
    
    run_id: str = Field(description="Pipeline run identifier")
    total_pairs: int = Field(description="Total compound-target pairs evaluated")
    
    # Score statistics
    docking_count: int = Field(default=0, description="Pairs with docking scores")
    evidence_count: int = Field(default=0, description="Pairs with experimental evidence")
    
    # Result summary
    top_hits: List[RankedHit] = Field(description="Top-ranked hits")
    score_distribution: Dict[str, Any] = Field(default_factory=dict, description="Score distribution statistics")
    
    # Processing metadata
    weights_used: Dict[str, float] = Field(description="Scoring weights")
    normalization_method: str = Field(description="Score normalization method")
    processing_time: float = Field(description="Total scoring time in seconds")
    
    created_at: datetime = Field(default_factory=datetime.now)
    
    class Config:
        extra = "allow"
