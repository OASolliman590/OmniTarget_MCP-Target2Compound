"""
Pydantic schemas for the drug discovery pipeline.
"""

from .config import RunConfig, ScoringConfig, PPIConfig, ExpressionConfig, StructureConfig, CompoundConfig
from .targets import Target, PathwayHit, ExpressionRecord, PPIRecord, SequenceRecord, StructureRecord
from .compounds import Compound, CompoundFeatures
from .scoring import DeepDTAScore, DockingScore, EvidenceRecord, IntegratedScore, RankedHit
from .pipeline import PipelineRun, PipelineStatus, PipelineResult
from .api import (
    RunRequest, 
    RunResponse, 
    StatusResponse, 
    ResultsResponse, 
    HealthResponse, 
    ErrorResponse
)

__all__ = [
    # Config
    "RunConfig",
    "ScoringConfig", 
    "PPIConfig",
    "ExpressionConfig",
    "StructureConfig",
    "CompoundConfig",
    
    # Targets
    "Target",
    "PathwayHit", 
    "ExpressionRecord",
    "PPIRecord",
    "SequenceRecord",
    "StructureRecord",
    
    # Compounds
    "Compound",
    "CompoundFeatures",
    
    # Scoring
    "DeepDTAScore",
    "DockingScore", 
    "EvidenceRecord",
    "IntegratedScore",
    "RankedHit",
    
    # Pipeline
    "PipelineRun",
    "PipelineStatus",
    "PipelineResult",
    
    # API
    "RunRequest",
    "RunResponse",
    "StatusResponse", 
    "ResultsResponse",
    "HealthResponse",
    "ErrorResponse",
]
