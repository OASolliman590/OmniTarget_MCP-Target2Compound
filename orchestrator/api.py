"""
FastAPI application for the drug discovery pipeline.
"""

import asyncio
import uuid
from datetime import datetime
from typing import Dict, List, Optional

from fastapi import FastAPI, HTTPException, BackgroundTasks, status
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse
import uvicorn

from .settings import settings
from .schemas.api import (
    RunRequest, RunResponse, StatusResponse, ResultsResponse,
    HealthResponse, ErrorResponse, ValidateConfigRequest, ValidateConfigResponse
)
from .schemas.config import RunConfig
from .schemas.pipeline import PipelineRun, PipelineStatus
from .pipeline import DrugDiscoveryPipeline

# Create FastAPI app
app = FastAPI(
    title=settings.api.title,
    description=settings.api.description,
    version=settings.api.version,
    docs_url="/docs",
    redoc_url="/redoc"
)

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.api.allowed_hosts,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Global pipeline instance and run storage
pipeline = DrugDiscoveryPipeline()
active_runs: Dict[str, PipelineRun] = {}
run_results: Dict[str, dict] = {}


@app.on_event("startup")
async def startup_event():
    """Initialize the pipeline on startup."""
    await pipeline.setup()


@app.get("/health", response_model=HealthResponse)
async def health_check():
    """Health check endpoint."""
    
    # Check MCP services
    services = {}
    try:
        services["kegg-mcp"] = await pipeline.kegg_client.test_connection()
        services["reactome-mcp"] = await pipeline.reactome_client.test_connection()
        services["proteinatlas-mcp"] = await pipeline.proteinatlas_client.test_connection()
        services["string-mcp"] = await pipeline.string_client.test_connection()
        services["uniprot-mcp"] = await pipeline.uniprot_client.test_connection()
        services["pdb-mcp"] = await pipeline.pdb_client.test_connection()
        services["chembl-mcp"] = await pipeline.chembl_client.test_connection()
    except Exception:
        # If any service check fails, mark all as unknown
        services = {name: False for name in services}
    
    # Determine overall status
    all_services_up = all(services.values())
    status_text = "healthy" if all_services_up else "degraded"
    
    return HealthResponse(
        status=status_text,
        version=settings.api.version,
        services=services,
        database=True,  # Placeholder
        storage=True,   # Placeholder
    )


@app.post("/runs", response_model=RunResponse)
async def create_run(
    request: RunRequest,
    background_tasks: BackgroundTasks
):
    """Start a new pipeline run."""
    
    try:
        # Generate run ID
        run_id = f"run_{datetime.now().strftime('%Y%m%d_%H%M%S')}_{str(uuid.uuid4())[:8]}"
        
        # Create run record
        run = PipelineRun(
            run_id=run_id,
            run_name=request.run_name or run_id,
            config_hash="",  # Will be set by pipeline
            config_data=request.config.dict(),
            status=PipelineStatus.PENDING
        )
        
        active_runs[run_id] = run
        
        # Start pipeline in background
        background_tasks.add_task(execute_pipeline, run_id, request.config)
        
        return RunResponse(
            run_id=run_id,
            status=PipelineStatus.PENDING,
            message="Pipeline run queued successfully",
            estimated_time=45  # Placeholder estimate
        )
        
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Failed to create run: {str(e)}"
        )


@app.get("/runs", response_model=List[PipelineRun])
async def list_runs(
    status_filter: Optional[PipelineStatus] = None,
    limit: int = 10,
    offset: int = 0
):
    """List pipeline runs."""
    
    runs = list(active_runs.values())
    
    if status_filter:
        runs = [run for run in runs if run.status == status_filter]
    
    # Sort by creation time (newest first)
    runs.sort(key=lambda x: x.created_at, reverse=True)
    
    # Apply pagination
    return runs[offset:offset + limit]


@app.get("/runs/{run_id}", response_model=StatusResponse)
async def get_run_status(run_id: str):
    """Get status of a specific run."""
    
    if run_id not in active_runs:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Run {run_id} not found"
        )
    
    run = active_runs[run_id]
    
    return StatusResponse(
        run_id=run_id,
        status=run.status,
        current_stage=run.current_stage,
        progress_percent=run.progress_percent,
        created_at=run.created_at,
        started_at=run.started_at,
        error_message=run.error_message
    )


@app.get("/runs/{run_id}/results", response_model=ResultsResponse)
async def get_run_results(run_id: str):
    """Get results of a completed run."""
    
    if run_id not in active_runs:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Run {run_id} not found"
        )
    
    run = active_runs[run_id]
    
    if run.status != PipelineStatus.COMPLETED:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Run {run_id} is not completed"
        )
    
    if run_id not in run_results:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Results for run {run_id} not found"
        )
    
    result_data = run_results[run_id]
    
    return ResultsResponse(
        run_id=run_id,
        status=run.status,
        result_summary=result_data["summary"],
        top_hits=result_data["top_hits"],
        download_urls={
            "results_csv": f"/runs/{run_id}/files/results.csv",
            "manifest": f"/runs/{run_id}/files/manifest.json"
        }
    )


@app.get("/runs/{run_id}/files/{filename}")
async def download_file(run_id: str, filename: str):
    """Download result files."""
    
    if run_id not in active_runs:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Run {run_id} not found"
        )
    
    # Construct file path (simplified)
    file_path = f"{settings.files.output_dir}/{run_id}/{filename}"
    
    try:
        return FileResponse(file_path)
    except FileNotFoundError:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"File {filename} not found"
        )


@app.post("/validate", response_model=ValidateConfigResponse)
async def validate_config(request: ValidateConfigRequest):
    """Validate a pipeline configuration."""
    
    try:
        config = request.config
        errors = []
        warnings = []
        
        # Basic validation
        if not config.disease_terms:
            errors.append("At least one disease term is required")
        
        if not config.compounds.input_paths:
            errors.append("At least one compound input path is required")
        
        # Check file existence if requested
        if request.check_files:
            from pathlib import Path
            for file_path in config.compounds.input_paths:
                if not Path(file_path).exists():
                    errors.append(f"Compound file not found: {file_path}")
        
        # Check scoring weights
        total_weight = sum(config.scoring.weights.values())
        if abs(total_weight - 1.0) > 0.01:
            errors.append(f"Scoring weights must sum to 1.0, got {total_weight}")
        
        # Generate warnings
        if len(config.compounds.input_paths) > 5:
            warnings.append("Large number of compound files may increase runtime")
        
        is_valid = len(errors) == 0
        
        return ValidateConfigResponse(
            is_valid=is_valid,
            errors=errors,
            warnings=warnings,
            estimated_time=45,  # Placeholder
            estimated_memory=8.5,  # Placeholder
            estimated_disk=2.1   # Placeholder
        )
        
    except Exception as e:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Validation error: {str(e)}"
        )


@app.delete("/runs/{run_id}")
async def cancel_run(run_id: str):
    """Cancel a running pipeline."""
    
    if run_id not in active_runs:
        raise HTTPException(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Run {run_id} not found"
        )
    
    run = active_runs[run_id]
    
    if run.status in [PipelineStatus.COMPLETED, PipelineStatus.FAILED]:
        raise HTTPException(
            status_code=status.HTTP_400_BAD_REQUEST,
            detail=f"Cannot cancel {run.status.value} run"
        )
    
    # Mark as cancelled
    run.status = PipelineStatus.CANCELLED
    run.completed_at = datetime.now()
    
    return {"message": f"Run {run_id} cancelled"}


async def execute_pipeline(run_id: str, config: RunConfig):
    """Execute pipeline in background."""
    
    try:
        run = active_runs[run_id]
        run.status = PipelineStatus.RUNNING
        run.started_at = datetime.now()
        
        # Execute pipeline
        result = await pipeline.run(config)
        
        # Store results
        run_results[run_id] = {
            "summary": result,
            "top_hits": []  # Would be populated from actual results
        }
        
        run.status = PipelineStatus.COMPLETED
        run.completed_at = datetime.now()
        
    except Exception as e:
        run = active_runs[run_id]
        run.status = PipelineStatus.FAILED
        run.error_message = str(e)
        run.completed_at = datetime.now()


if __name__ == "__main__":
    uvicorn.run(
        "orchestrator.api:app",
        host=settings.api.host,
        port=settings.api.port,
        reload=settings.api.debug
    )
