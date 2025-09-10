"""
FastAPI request and response schemas.
"""

from typing import List, Optional, Dict, Any
from pydantic import BaseModel, Field
from datetime import datetime

from .config import RunConfig
from .pipeline import PipelineStatus, PipelineStage, PipelineRun, PipelineResult
from .scoring import RankedHit


class RunRequest(BaseModel):
    """Request to start a new pipeline run."""
    
    config: RunConfig = Field(description="Pipeline configuration")
    run_name: Optional[str] = Field(default=None, description="Optional run name")
    priority: int = Field(default=1, ge=1, le=10, description="Run priority (1-10)")
    notify_email: Optional[str] = Field(default=None, description="Email for notifications")
    
    class Config:
        schema_extra = {
            "example": {
                "config": {
                    "disease_terms": ["lung cancer"],
                    "organism": "Homo sapiens",
                    "compounds": {
                        "input_paths": ["data/compounds/test_set.smi"]
                    },
                    "scoring": {
                        "weights": {
                            "similarity": 0.5,
                            "pharmacophore": 0.2,
                            "docking": 0.1,
                            "evidence": 0.2
                        }
                    }
                },
                "run_name": "lung_cancer_pilot_study",
                "priority": 5
            }
        }


class RunResponse(BaseModel):
    """Response when starting a pipeline run."""
    
    run_id: str = Field(description="Unique run identifier")
    status: PipelineStatus = Field(description="Initial run status")
    message: str = Field(description="Status message")
    estimated_time: Optional[int] = Field(default=None, description="Estimated completion time (minutes)")
    
    class Config:
        schema_extra = {
            "example": {
                "run_id": "run_20240115_143022_abc123",
                "status": "pending",
                "message": "Pipeline run queued successfully",
                "estimated_time": 45
            }
        }


class StatusResponse(BaseModel):
    """Response for pipeline status queries."""
    
    run_id: str = Field(description="Run identifier")
    status: PipelineStatus = Field(description="Current status")
    current_stage: Optional[PipelineStage] = Field(default=None, description="Current stage")
    progress_percent: float = Field(description="Progress percentage")
    
    # Timing information
    created_at: datetime = Field(description="Run creation time")
    started_at: Optional[datetime] = Field(default=None, description="Run start time")
    estimated_completion: Optional[datetime] = Field(default=None, description="Estimated completion")
    
    # Error information
    error_message: Optional[str] = Field(default=None, description="Error message if failed")
    warnings: List[str] = Field(default_factory=list, description="Warnings")
    
    # Resource usage
    cpu_time: Optional[float] = Field(default=None, description="CPU time used (seconds)")
    memory_usage: Optional[float] = Field(default=None, description="Current memory usage (MB)")
    
    class Config:
        schema_extra = {
            "example": {
                "run_id": "run_20240115_143022_abc123",
                "status": "running",
                "current_stage": "feature_generation",
                "progress_percent": 65.5,
                "created_at": "2024-01-15T14:30:22Z",
                "started_at": "2024-01-15T14:30:25Z",
                "estimated_completion": "2024-01-15T15:15:00Z"
            }
        }


class ResultsResponse(BaseModel):
    """Response containing pipeline results."""
    
    run_id: str = Field(description="Run identifier")
    status: PipelineStatus = Field(description="Final status")
    result_summary: PipelineResult = Field(description="Result summary")
    top_hits: List[RankedHit] = Field(description="Top-ranked hits")
    
    # File information
    download_urls: Dict[str, str] = Field(default_factory=dict, description="Download URLs for result files")
    
    class Config:
        schema_extra = {
            "example": {
                "run_id": "run_20240115_143022_abc123",
                "status": "completed",
                "result_summary": {
                    "run_id": "run_20240115_143022_abc123",
                    "total_targets": 25,
                    "total_compounds": 100,
                    "total_pairs": 2500,
                    "top_hits_count": 10,
                    "best_score": 0.95
                },
                "download_urls": {
                    "results_csv": "/api/v1/runs/run_20240115_143022_abc123/files/results.csv",
                    "manifest": "/api/v1/runs/run_20240115_143022_abc123/files/manifest.json"
                }
            }
        }


class RunListResponse(BaseModel):
    """Response for listing pipeline runs."""
    
    runs: List[PipelineRun] = Field(description="List of pipeline runs")
    total_count: int = Field(description="Total number of runs")
    page: int = Field(description="Current page number")
    page_size: int = Field(description="Number of runs per page")
    
    class Config:
        schema_extra = {
            "example": {
                "runs": [
                    {
                        "run_id": "run_20240115_143022_abc123",
                        "status": "completed",
                        "created_at": "2024-01-15T14:30:22Z",
                        "progress_percent": 100.0
                    }
                ],
                "total_count": 1,
                "page": 1,
                "page_size": 10
            }
        }


class HealthResponse(BaseModel):
    """Health check response."""
    
    status: str = Field(description="Service status")
    timestamp: datetime = Field(default_factory=datetime.now, description="Check timestamp")
    version: str = Field(description="Application version")
    
    # Service availability
    services: Dict[str, bool] = Field(description="MCP service availability")
    database: bool = Field(description="Database connectivity")
    storage: bool = Field(description="Storage accessibility")
    
    # Resource status
    cpu_usage: Optional[float] = Field(default=None, description="CPU usage percentage")
    memory_usage: Optional[float] = Field(default=None, description="Memory usage percentage")
    disk_usage: Optional[float] = Field(default=None, description="Disk usage percentage")
    
    class Config:
        schema_extra = {
            "example": {
                "status": "healthy",
                "timestamp": "2024-01-15T14:30:22Z",
                "version": "0.2.1",
                "services": {
                    "kegg-mcp": True,
                    "reactome-mcp": True,
                    "chembl-mcp": True
                },
                "database": True,
                "storage": True,
                "cpu_usage": 25.5,
                "memory_usage": 45.2,
                "disk_usage": 12.8
            }
        }


class ErrorResponse(BaseModel):
    """Error response schema."""
    
    error: str = Field(description="Error type")
    message: str = Field(description="Error message")
    details: Optional[Dict[str, Any]] = Field(default=None, description="Additional error details")
    timestamp: datetime = Field(default_factory=datetime.now, description="Error timestamp")
    request_id: Optional[str] = Field(default=None, description="Request identifier for tracking")
    
    class Config:
        schema_extra = {
            "example": {
                "error": "ValidationError",
                "message": "Invalid configuration provided",
                "details": {
                    "field": "compounds.input_paths",
                    "issue": "File does not exist"
                },
                "timestamp": "2024-01-15T14:30:22Z",
                "request_id": "req_abc123"
            }
        }


class ValidateConfigRequest(BaseModel):
    """Request to validate a configuration."""
    
    config: RunConfig = Field(description="Configuration to validate")
    check_files: bool = Field(default=True, description="Whether to check file existence")
    
    class Config:
        schema_extra = {
            "example": {
                "config": {
                    "disease_terms": ["lung cancer"],
                    "compounds": {
                        "input_paths": ["data/compounds/test_set.smi"]
                    }
                },
                "check_files": True
            }
        }


class ValidateConfigResponse(BaseModel):
    """Response for configuration validation."""
    
    is_valid: bool = Field(description="Whether configuration is valid")
    errors: List[str] = Field(default_factory=list, description="Validation errors")
    warnings: List[str] = Field(default_factory=list, description="Validation warnings")
    
    # Estimated resource requirements
    estimated_time: Optional[int] = Field(default=None, description="Estimated runtime (minutes)")
    estimated_memory: Optional[float] = Field(default=None, description="Estimated memory usage (GB)")
    estimated_disk: Optional[float] = Field(default=None, description="Estimated disk usage (GB)")
    
    class Config:
        schema_extra = {
            "example": {
                "is_valid": True,
                "errors": [],
                "warnings": ["Large number of compounds may increase runtime"],
                "estimated_time": 45,
                "estimated_memory": 8.5,
                "estimated_disk": 2.1
            }
        }
