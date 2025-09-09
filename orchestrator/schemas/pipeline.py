"""
Pipeline execution and status schemas.
"""

from typing import List, Optional, Dict, Any
from pydantic import BaseModel, Field
from datetime import datetime
from enum import Enum


class PipelineStatus(str, Enum):
    """Pipeline execution status."""
    PENDING = "pending"
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"
    CANCELLED = "cancelled"


class PipelineStage(str, Enum):
    """Pipeline execution stages."""
    INITIALIZATION = "initialization"
    TARGET_DISCOVERY = "target_discovery"
    TARGET_CHARACTERIZATION = "target_characterization"
    COMPOUND_LOADING = "compound_loading"
    FEATURE_GENERATION = "feature_generation"
    DEEPDTA_SCORING = "deepdta_scoring"
    DOCKING = "docking"
    EVIDENCE_VALIDATION = "evidence_validation"
    SCORE_INTEGRATION = "score_integration"
    RESULT_GENERATION = "result_generation"
    FINALIZATION = "finalization"


class StageResult(BaseModel):
    """Result of a pipeline stage execution."""
    
    stage: PipelineStage = Field(description="Stage identifier")
    status: PipelineStatus = Field(description="Stage execution status")
    
    # Timing
    started_at: Optional[datetime] = Field(default=None, description="Stage start time")
    completed_at: Optional[datetime] = Field(default=None, description="Stage completion time")
    duration: Optional[float] = Field(default=None, description="Duration in seconds")
    
    # Results
    items_processed: int = Field(default=0, description="Number of items processed")
    items_successful: int = Field(default=0, description="Number of successful items")
    items_failed: int = Field(default=0, description="Number of failed items")
    
    # Error handling
    error_message: Optional[str] = Field(default=None, description="Error message if failed")
    warnings: List[str] = Field(default_factory=list, description="Stage warnings")
    
    # Stage-specific data
    metadata: Dict[str, Any] = Field(default_factory=dict, description="Stage-specific metadata")
    
    class Config:
        extra = "allow"


class PipelineRun(BaseModel):
    """A complete pipeline execution record."""
    
    # Identifiers
    run_id: str = Field(description="Unique run identifier")
    run_name: Optional[str] = Field(default=None, description="User-provided run name")
    
    # Configuration
    config_hash: str = Field(description="Hash of configuration used")
    config_data: Dict[str, Any] = Field(description="Full configuration data")
    
    # Status and timing
    status: PipelineStatus = Field(default=PipelineStatus.PENDING)
    created_at: datetime = Field(default_factory=datetime.now)
    started_at: Optional[datetime] = Field(default=None)
    completed_at: Optional[datetime] = Field(default=None)
    duration: Optional[float] = Field(default=None, description="Total duration in seconds")
    
    # Progress tracking
    current_stage: Optional[PipelineStage] = Field(default=None)
    progress_percent: float = Field(default=0.0, ge=0.0, le=100.0)
    stage_results: List[StageResult] = Field(default_factory=list)
    
    # Resource usage
    cpu_time: Optional[float] = Field(default=None, description="CPU time used")
    memory_peak: Optional[float] = Field(default=None, description="Peak memory usage (MB)")
    disk_usage: Optional[float] = Field(default=None, description="Disk space used (MB)")
    
    # Error handling
    error_message: Optional[str] = Field(default=None, description="Global error message")
    warnings: List[str] = Field(default_factory=list, description="Global warnings")
    
    # Output information
    output_dir: Optional[str] = Field(default=None, description="Output directory path")
    result_files: List[str] = Field(default_factory=list, description="Generated result files")
    
    class Config:
        extra = "allow"
    
    def update_progress(self) -> None:
        """Update progress based on completed stages."""
        total_stages = len(PipelineStage)
        completed_stages = sum(1 for stage_result in self.stage_results 
                             if stage_result.status == PipelineStatus.COMPLETED)
        self.progress_percent = (completed_stages / total_stages) * 100
    
    def get_stage_result(self, stage: PipelineStage) -> Optional[StageResult]:
        """Get result for a specific stage."""
        for stage_result in self.stage_results:
            if stage_result.stage == stage:
                return stage_result
        return None
    
    def add_stage_result(self, stage_result: StageResult) -> None:
        """Add or update a stage result."""
        # Remove existing result for this stage
        self.stage_results = [sr for sr in self.stage_results if sr.stage != stage_result.stage]
        # Add new result
        self.stage_results.append(stage_result)
        self.update_progress()


class PipelineResult(BaseModel):
    """Final results of a pipeline execution."""
    
    run_id: str = Field(description="Pipeline run identifier")
    
    # Summary statistics
    total_targets: int = Field(description="Total targets discovered")
    total_compounds: int = Field(description="Total compounds processed")
    total_pairs: int = Field(description="Total compound-target pairs evaluated")
    
    # Results by method
    deepdta_predictions: int = Field(default=0, description="DeepDTA predictions made")
    docking_results: int = Field(default=0, description="Docking calculations completed")
    evidence_matches: int = Field(default=0, description="Evidence records found")
    
    # Top results
    top_hits_count: int = Field(description="Number of top hits returned")
    best_score: Optional[float] = Field(default=None, description="Best combined score")
    score_range: Optional[List[float]] = Field(default=None, description="Score range [min, max]")
    
    # Quality metrics
    high_confidence_hits: int = Field(default=0, description="High confidence hits")
    medium_confidence_hits: int = Field(default=0, description="Medium confidence hits")
    low_confidence_hits: int = Field(default=0, description="Low confidence hits")
    
    # File outputs
    results_file: Optional[str] = Field(default=None, description="Main results CSV file")
    manifest_file: Optional[str] = Field(default=None, description="Provenance manifest")
    plots_dir: Optional[str] = Field(default=None, description="Plots directory")
    
    # Timing and resources
    total_time: float = Field(description="Total execution time in seconds")
    target_discovery_time: Optional[float] = Field(default=None, description="Target discovery time")
    scoring_time: Optional[float] = Field(default=None, description="Scoring time")
    
    # Provenance
    environment_info: Dict[str, Any] = Field(default_factory=dict, description="Environment information")
    dependencies: Dict[str, str] = Field(default_factory=dict, description="Dependency versions")
    
    created_at: datetime = Field(default_factory=datetime.now)
    
    class Config:
        extra = "allow"
