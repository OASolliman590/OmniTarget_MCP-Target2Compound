"""
Command-line interface for the drug discovery pipeline.
"""

import asyncio
import sys
from pathlib import Path
from typing import Optional

import typer
import yaml
from loguru import logger
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn
from rich.table import Table

from .settings import settings
from .schemas.config import RunConfig
# Import legacy pipeline lazily to avoid unnecessary overhead
from .pipeline_nondocking import NonDockingPipeline

# Create Typer app
app = typer.Typer(
    name="mcp-orchestrator",
    help="MCP Drug Discovery Pipeline Orchestrator",
    add_completion=False
)

console = Console()


@app.command()
def run(
    config_path: str = typer.Argument(..., help="Path to configuration YAML file (non-docking by default)"),
    legacy: bool = typer.Option(False, "--legacy", help="Use the legacy docking-enabled pipeline"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Validate configuration without running (legacy only)"),
    debug: bool = typer.Option(False, "--debug", help="Enable debug logging"),
    output_dir: Optional[str] = typer.Option(None, "--output-dir", help="Override output directory (legacy only)")
):
    """Run the pipeline. Defaults to non-docking; use --legacy for the old flow."""

    # Configure logging
    log_level = "DEBUG" if debug else settings.logging.log_level
    logger.remove()
    logger.add(sys.stderr, level=log_level, format=settings.logging.log_format)

    if legacy:
        logger.info("'run' executing legacy pipeline. For non-docking use 'run_nd'.")
        # Legacy path mirrors previous behavior
        try:
            from .pipeline import DrugDiscoveryPipeline
            config = load_config(config_path)
            if output_dir:
                config.output_dir = output_dir
            if dry_run:
                config.dry_run = True
            asyncio.run(run_pipeline(config))
        except Exception as e:
            logger.error(f"Legacy pipeline failed: {e}")
            raise typer.Exit(1)
        return

    # Non-legacy path: non-docking default
    logger.info("'run' now executes the non-docking pipeline. Use '--legacy' or 'run_legacy' for the old behavior.")
    try:
        with open(config_path) as f:
            cfg = yaml.safe_load(f)
        pipeline = NonDockingPipeline(cfg)
        result = pipeline.run()
        console.print(f"[green]✓[/green] Non-docking run completed. Results: {result['results_csv']}")
    except Exception as e:
        logger.error(f"Non-docking pipeline failed: {e}")
        raise typer.Exit(1)


@app.command()
def run_nd(
    config_path: str = typer.Argument(..., help="Path to non-docking YAML config"),
    enable_qsar: bool = typer.Option(False, "--enable-qsar", help="Enable QSAR stage (if feasible)"),
    disable_ph4: bool = typer.Option(False, "--disable-ph4", help="Disable pharmacophore features"),
    debug: bool = typer.Option(False, "--debug", help="Enable debug logging"),
):
    """Run the non-docking pipeline (comparators + similarity + pharmacophore)."""

    # Configure logging
    log_level = "DEBUG" if debug else settings.logging.log_level
    logger.remove()
    logger.add(sys.stderr, level=log_level, format=settings.logging.log_format)

    try:
        with open(config_path) as f:
            cfg = yaml.safe_load(f)

        # CLI overrides
        if enable_qsar:
            cfg.setdefault("scoring", {}).setdefault("weights", {}).setdefault("qsar", 0.15)
        if disable_ph4:
            cfg.setdefault("pharmacophore", {})["method"] = "disabled"

        pipeline = NonDockingPipeline(cfg)
        result = pipeline.run()
        console.print(f"[green]✓[/green] Non-docking run completed. Results: {result['results_csv']}")
    except Exception as e:
        logger.error(f"Non-docking pipeline failed: {e}")
        raise typer.Exit(1)


@app.command()
def validate(
    config_path: str = typer.Argument(..., help="Path to configuration YAML file"),
    check_files: bool = typer.Option(True, "--check-files/--no-check-files", help="Check file existence")
):
    """Validate pipeline configuration."""
    
    try:
        config = load_config(config_path)
        
        console.print("[green]✓[/green] Configuration loaded successfully")
        
        # Validate configuration
        errors = validate_config(config, check_files)
        
        if errors:
            console.print("\n[red]Configuration Errors:[/red]")
            for error in errors:
                console.print(f"  • {error}")
            raise typer.Exit(1)
        else:
            console.print("[green]✓[/green] Configuration is valid")
            
    except Exception as e:
        console.print(f"[red]Error:[/red] {e}")
        raise typer.Exit(1)


@app.command()
def status():
    """Check status of MCP servers and dependencies."""
    
    async def check_status():
        from .pipeline import DrugDiscoveryPipeline
        pipeline = DrugDiscoveryPipeline()
        
        # Create status table
        table = Table(title="MCP Services Status")
        table.add_column("Service", style="cyan")
        table.add_column("Status", style="bold")
        table.add_column("URL")
        
        # Check MCP clients
        clients = [
            ("KEGG", pipeline.kegg_client),
            ("Reactome", pipeline.reactome_client),
            ("ProteinAtlas", pipeline.proteinatlas_client),
            ("STRING", pipeline.string_client),
            ("UniProt", pipeline.uniprot_client),
            ("PDB", pipeline.pdb_client),
            ("ChEMBL", pipeline.chembl_client),
        ]
        
        for name, client in clients:
            try:
                is_healthy = await client.test_connection()
                status = "[green]✓ Online[/green]" if is_healthy else "[red]✗ Offline[/red]"
                table.add_row(name, status, client.base_url)
            except Exception as e:
                table.add_row(name, f"[red]✗ Error[/red]", str(e))
        
        console.print(table)
    
    asyncio.run(check_status())


@app.command()
def run_legacy(
    config_path: str = typer.Argument(..., help="Path to configuration YAML file (legacy pipeline)"),
    dry_run: bool = typer.Option(False, "--dry-run", help="Validate configuration without running"),
    debug: bool = typer.Option(False, "--debug", help="Enable debug logging"),
    output_dir: Optional[str] = typer.Option(None, "--output-dir", help="Override output directory")
):
    """Run the legacy docking-enabled pipeline explicitly."""

    # Configure logging
    log_level = "DEBUG" if debug else settings.logging.log_level
    logger.remove()
    logger.add(sys.stderr, level=log_level, format=settings.logging.log_format)

    try:
        from .pipeline import DrugDiscoveryPipeline
        config = load_config(config_path)
        if output_dir:
            config.output_dir = output_dir
        if dry_run:
            config.dry_run = True
        asyncio.run(run_pipeline(config))
    except Exception as e:
        logger.error(f"Legacy pipeline failed: {e}")
        raise typer.Exit(1)


@app.command()
def example_config(
    output_path: str = typer.Option("config.yaml", help="Output path for example configuration")
):
    """Generate example configuration file."""
    
    example_config = {
        "disease_terms": ["lung cancer"],
        "organism": "Homo sapiens",
        "ppi": {
            "min_confidence": 0.7,
            "max_partners": 200
        },
        "expression": {
            "tissue_whitelist": ["lung"],
            "reliability_min": "Supported"
        },
        "structures": {
            "prefer_experimental": True,
            "allow_alphafold": True
        },
        "compounds": {
            "input_paths": ["data/compounds/example.smi"]
        },
        "scoring": {
            "weights": {
                "similarity": 0.5,
                "pharmacophore": 0.2,
                "docking": 0.1,
                "evidence": 0.2
            }
        },
        "output_dir": "data/outputs"
    }
    
    with open(output_path, "w") as f:
        yaml.dump(example_config, f, default_flow_style=False, sort_keys=False)
    
    console.print(f"[green]✓[/green] Example configuration written to {output_path}")


@app.command()
def example_nd_config(
    output_path: str = typer.Option("run.nd.yaml", help="Output path for non-docking example configuration")
):
    """Generate non-docking example configuration file."""

    from pathlib import Path
    src = Path("configs/run.nd.example.yaml")
    if not src.exists():
        console.print("[red]Error:[/red] run.nd.example.yaml not found in configs/")
        raise typer.Exit(1)
    content = src.read_text()
    Path(output_path).write_text(content)
    console.print(f"[green]✓[/green] Non-docking example configuration written to {output_path}")


def load_config(config_path: str) -> RunConfig:
    """Load and parse configuration file."""
    
    config_file = Path(config_path)
    
    if not config_file.exists():
        raise typer.BadParameter(f"Configuration file not found: {config_path}")
    
    try:
        with open(config_file) as f:
            config_data = yaml.safe_load(f)
        
        config = RunConfig(**config_data)
        return config
        
    except yaml.YAMLError as e:
        raise typer.BadParameter(f"Invalid YAML in configuration file: {e}")
    except Exception as e:
        raise typer.BadParameter(f"Invalid configuration: {e}")


def validate_config(config: RunConfig, check_files: bool = True) -> list:
    """Validate configuration and return list of errors."""
    
    errors = []
    
    # Check compound files
    if check_files:
        for file_path in config.compounds.input_paths:
            if not Path(file_path).exists():
                errors.append(f"Compound file not found: {file_path}")
    
    # Check output directory is writable
    output_dir = Path(config.output_dir)
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
    except Exception as e:
        errors.append(f"Cannot create output directory {output_dir}: {e}")
    
    # Validate scoring weights
    total_weight = sum(config.scoring.weights.values())
    if abs(total_weight - 1.0) > 0.01:
        errors.append(f"Scoring weights must sum to 1.0, got {total_weight}")
    
    return errors


async def run_pipeline(config: RunConfig):
    """Execute the pipeline with progress tracking."""
    
    pipeline = DrugDiscoveryPipeline()
    
    # Setup pipeline
    console.print("[blue]Setting up pipeline...[/blue]")
    setup_success = await pipeline.setup()
    
    if not setup_success:
        console.print("[red]✗[/red] Pipeline setup failed")
        return
    
    console.print("[green]✓[/green] Pipeline setup completed")
    
    if config.dry_run:
        console.print("[yellow]Dry run mode - validation completed[/yellow]")
        return
    
    # Run pipeline with progress tracking
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        console=console,
    ) as progress:
        
        task = progress.add_task("Running pipeline...", total=None)
        
        try:
            result = await pipeline.run(config)
            
            progress.update(task, completed=True, description="[green]Pipeline completed![/green]")
            
            # Display results summary
            display_results_summary(result)
            
        except Exception as e:
            progress.update(task, description=f"[red]Pipeline failed: {e}[/red]")
            raise


def display_results_summary(result):
    """Display a summary of pipeline results."""
    
    table = Table(title="Pipeline Results Summary")
    table.add_column("Metric", style="cyan")
    table.add_column("Value", style="bold")
    
    table.add_row("Total Targets", str(result.total_targets))
    table.add_row("Total Compounds", str(result.total_compounds))
    table.add_row("Total Pairs Evaluated", str(result.total_pairs))
    # Sequence-based predictor removed; no predictions row
    table.add_row("Docking Results", str(result.docking_results))
    table.add_row("Evidence Matches", str(result.evidence_matches))
    table.add_row("Top Hits", str(result.top_hits_count))
    table.add_row("Best Score", f"{result.best_score:.4f}" if result.best_score else "N/A")
    table.add_row("Total Time", f"{result.total_time:.1f}s")
    
    console.print(table)
    
    if result.results_file:
        console.print(f"\n[green]Results saved to:[/green] {result.results_file}")
    if result.manifest_file:
        console.print(f"[green]Manifest saved to:[/green] {result.manifest_file}")


if __name__ == "__main__":
    app()
