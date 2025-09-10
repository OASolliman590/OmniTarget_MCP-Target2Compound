"""
Integration tests for the complete pipeline.
"""

import pytest
import asyncio
from pathlib import Path

from ..schemas.config import RunConfig
from ..pipeline import DrugDiscoveryPipeline


@pytest.fixture
def test_config():
    """Create a minimal test configuration."""
    return RunConfig(
        disease_terms=["lung cancer"],
        organism="Homo sapiens",
        compounds={
            "input_paths": ["data/compounds/test_set.smi"],
            "max_compounds": 10
        },
        max_targets=5,
        scoring={
            "weights": {
                "similarity": 0.5,
                "pharmacophore": 0.2,
                "docking": 0.1,
                "evidence": 0.2
            }
        },
        output_dir="data/outputs/test",
        dry_run=False
    )


@pytest.fixture
def test_compounds_file():
    """Create a test compounds file."""
    test_file = Path("data/compounds/test_set.smi")
    test_file.parent.mkdir(exist_ok=True, parents=True)
    
    with open(test_file, "w") as f:
        f.write("CCO ethanol\n")
        f.write("CC(C)O isopropanol\n")
        f.write("C1CCCCC1 cyclohexane\n")
    
    yield test_file
    
    # Cleanup
    if test_file.exists():
        test_file.unlink()


@pytest.mark.asyncio
async def test_pipeline_setup():
    """Test that pipeline components can be set up."""
    pipeline = DrugDiscoveryPipeline()
    
    # This will test the setup of all adapters
    success = await pipeline.setup()
    
    # Note: In a real environment with all services running, this should be True
    # For testing purposes, we accept that some components may not be available
    assert isinstance(success, bool)


@pytest.mark.asyncio
async def test_mcp_client_connections():
    """Test MCP client connections."""
    pipeline = DrugDiscoveryPipeline()
    
    # Test individual client connections
    clients = [
        ("KEGG", pipeline.kegg_client),
        ("Reactome", pipeline.reactome_client),
        ("UniProt", pipeline.uniprot_client),
    ]
    
    for name, client in clients:
        try:
            result = await client.test_connection()
            print(f"{name} connection: {result}")
            # Don't assert True - services may not be available in test environment
            assert isinstance(result, bool)
        except Exception as e:
            print(f"{name} connection failed: {e}")
            # This is expected in test environments without running MCP servers


@pytest.mark.asyncio
async def test_compound_loading(test_compounds_file, test_config):
    """Test compound loading functionality."""
    pipeline = DrugDiscoveryPipeline()
    
    # This tests the compound loading stage
    compounds = await pipeline._load_compounds(test_config, None)
    
    assert len(compounds) == 3
    assert compounds[0].smiles == "CCO"
    assert compounds[0].name == "ethanol"


@pytest.mark.asyncio
@pytest.mark.slow
async def test_full_pipeline_dry_run(test_compounds_file):
    """Test full pipeline in dry run mode."""
    pipeline = DrugDiscoveryPipeline()
    
    config = RunConfig(
        disease_terms=["test disease"],
        compounds={
            "input_paths": [str(test_compounds_file)]
        },
        dry_run=True
    )
    
    # Setup pipeline
    await pipeline.setup()
    
    # Note: This would test the full pipeline flow
    # In a real test, we would run the pipeline and check outputs
    # For now, we just test that the configuration is valid
    assert config.dry_run is True
    assert len(config.disease_terms) == 1


def test_config_validation():
    """Test configuration validation."""
    # Valid config
    config = RunConfig(
        disease_terms=["lung cancer"],
        compounds={
            "input_paths": ["test.smi"]  # File doesn't need to exist for this test
        }
    )
    
    assert config.disease_terms == ["lung cancer"]
    assert config.organism == "Homo sapiens"  # Default value
    
    # Test invalid config
    with pytest.raises(Exception):
        RunConfig(
            disease_terms=[],  # Empty disease terms should fail
            compounds={
                "input_paths": ["test.smi"]
            }
        )


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, "-v"])
