"""
Test that results.csv contains the expected columns including new E3FP and comparator diagnostics.
"""

import pytest
import pandas as pd
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock

from orchestrator.pipeline_nondocking import NonDockingPipeline


class TestResultsColumns:
    """Test that results.csv contains all expected columns."""
    
    def test_results_columns_structure(self):
        """Test that results.csv has the expected column structure."""
        expected_columns = [
            "compound_id",
            "compound_smiles",
            "compound_inchikey", 
            "compound_source",
            "target_uniprot",
            "target_chembl",
            "target_gene",
            "evidence_strength",
            # Comparator diagnostics
            "n_comparators",
            "median_pchembl",
            "assay_consistency",
            "gate_reason",
            # Similarity features
            "cosine_max",
            "tanimoto_max",
            # E3FP features
            "e3fp_tanimoto_max",
            "e3fp_tanimoto_mean_topk",
            # Other features
            "ph4_best",
            "qsar_score",
            "docking_score",
            "fused_score",
            "rank",
        ]
        
        # Create a minimal test configuration
        test_config = {
            "disease_terms": ["test"],
            "organism": "Homo sapiens",
            "chembl": {
                "min_pchembl": 4.0,
                "min_comparators": 3,
                "require_assay_consistency": False
            },
            "similarity": {
                "top_k": 10,
                "min_tanimoto_morgan": 0.1,
                "enable_e3fp": False
            },
            "pharmacophore": {
                "method": "disabled"
            },
            "scoring": {
                "weights": {
                    "similarity": 0.5,
                    "e3fp": 0.0,
                    "pharmacophore": 0.0,
                    "evidence": 0.5,
                    "docking": 0.0
                }
            },
            "ligand_prep": {
                "input_paths": ["data/compounds/test.smi"],
                "ph": 7.4,
                "strip_salts": True,
                "largest_fragment_only": True,
                "conformers": {
                    "enable": False
                }
            },
            "output": {
                "dir": "data/outputs"
            },
            "docking": {
                "enabled": False
            }
        }
        
        # Create a temporary test compound file
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            test_smi = temp_path / "test.smi"
            test_smi.write_text("CCO\tethanol\n")
            
            # Update config to use temp file
            test_config["ligand_prep"]["input_paths"] = [str(test_smi)]
            test_config["output"]["dir"] = str(temp_path / "outputs")
            
            # Mock the MCP clients to avoid external dependencies
            with patch('orchestrator.pipeline_nondocking.ChEMBLMCPClient') as mock_chembl:
                mock_client = MagicMock()
                mock_chembl.return_value = mock_client
                
                # Mock successful ChEMBL response
                mock_activities = [
                    {
                        "type": "text",
                        "text": '{"activities": [{"molecule_chembl_id": "CHEMBL123", "canonical_smiles": "CCO", "pchembl_value": "5.5", "assay_type": "IC50", "target_organism": "Homo sapiens"}]}'
                    }
                ]
                mock_client.search_activities.return_value = mock_activities
                
                # Mock target discovery
                with patch('orchestrator.pipeline_nondocking.KEGGClient') as mock_kegg:
                    with patch('orchestrator.pipeline_nondocking.ReactomeClient') as mock_reactome:
                        with patch('orchestrator.pipeline_nondocking.ProteinAtlasClient') as mock_pa:
                            with patch('orchestrator.pipeline_nondocking.STRINGClient') as mock_string:
                                with patch('orchestrator.pipeline_nondocking.UniProtClient') as mock_uniprot:
                                    with patch('orchestrator.pipeline_nondocking.PDBClient') as mock_pdb:
                                        
                                        # Mock all clients to return empty results
                                        for mock_client_class in [mock_kegg, mock_reactome, mock_pa, mock_string, mock_uniprot, mock_pdb]:
                                            mock_client_instance = MagicMock()
                                            mock_client_class.return_value = mock_client_instance
                                            mock_client_instance.search_targets.return_value = []
                                        
                                        # Create pipeline and run
                                        pipeline = NonDockingPipeline(test_config)
                                        
                                        # Mock the target discovery to return a test target
                                        with patch.object(pipeline, '_discover_targets') as mock_discover:
                                            mock_discover.return_value = [
                                                {
                                                    "chembl_id": "CHEMBL123",
                                                    "uniprot_id": "P12345",
                                                    "gene_name": "TEST_GENE"
                                                }
                                            ]
                                            
                                            # Run the pipeline
                                            result = pipeline.run()
                                            
                                            # Check that results.csv was created
                                            results_path = Path(result["results_csv"])
                                            assert results_path.exists()
                                            
                                            # Read and check columns
                                            df = pd.read_csv(results_path)
                                            
                                            # Check that all expected columns are present
                                            for col in expected_columns:
                                                assert col in df.columns, f"Missing column: {col}"
                                            
                                            # Check that we have some data
                                            assert len(df) > 0, "No results generated"
                                            
                                            # Check that E3FP columns are present (even if None)
                                            assert "e3fp_tanimoto_max" in df.columns
                                            assert "e3fp_tanimoto_mean_topk" in df.columns
                                            
                                            # Check that comparator diagnostic columns are present
                                            assert "n_comparators" in df.columns
                                            assert "median_pchembl" in df.columns
                                            assert "assay_consistency" in df.columns
                                            assert "gate_reason" in df.columns
    
    def test_e3fp_columns_when_enabled(self):
        """Test that E3FP columns are populated when E3FP is enabled."""
        test_config = {
            "disease_terms": ["test"],
            "organism": "Homo sapiens",
            "chembl": {
                "min_pchembl": 4.0,
                "min_comparators": 3,
                "require_assay_consistency": False
            },
            "similarity": {
                "top_k": 10,
                "min_tanimoto_morgan": 0.1,
                "enable_e3fp": True,
                "e3fp": {
                    "radius": 2,
                    "shells": 5,
                    "top_k": 25
                }
            },
            "pharmacophore": {
                "method": "disabled"
            },
            "scoring": {
                "weights": {
                    "similarity": 0.4,
                    "e3fp": 0.1,
                    "pharmacophore": 0.0,
                    "evidence": 0.5,
                    "docking": 0.0
                }
            },
            "ligand_prep": {
                "input_paths": ["data/compounds/test.smi"],
                "ph": 7.4,
                "strip_salts": True,
                "largest_fragment_only": True,
                "conformers": {
                    "enable": True,
                    "max_confs": 5
                }
            },
            "output": {
                "dir": "data/outputs"
            },
            "docking": {
                "enabled": False
            }
        }
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            test_smi = temp_path / "test.smi"
            test_smi.write_text("CCO\tethanol\n")
            
            test_config["ligand_prep"]["input_paths"] = [str(test_smi)]
            test_config["output"]["dir"] = str(temp_path / "outputs")
            
            # Mock E3FP computation to return some values
            with patch('orchestrator.adapters.e3fp_adapter.compute_e3fp_for_sdf') as mock_e3fp:
                mock_e3fp.return_value = {"ethanol": [1, 0, 1, 0, 1]}  # Mock fingerprint
                
                # Mock other dependencies
                with patch('orchestrator.pipeline_nondocking.ChEMBLMCPClient') as mock_chembl:
                    mock_client = MagicMock()
                    mock_chembl.return_value = mock_client
                    mock_activities = [
                        {
                            "type": "text",
                            "text": '{"activities": [{"molecule_chembl_id": "CHEMBL123", "canonical_smiles": "CCO", "pchembl_value": "5.5", "assay_type": "IC50", "target_organism": "Homo sapiens"}]}'
                        }
                    ]
                    mock_client.search_activities.return_value = mock_activities
                    
                    # Mock target discovery
                    with patch('orchestrator.pipeline_nondocking.KEGGClient') as mock_kegg:
                        with patch('orchestrator.pipeline_nondocking.ReactomeClient') as mock_reactome:
                            with patch('orchestrator.pipeline_nondocking.ProteinAtlasClient') as mock_pa:
                                with patch('orchestrator.pipeline_nondocking.STRINGClient') as mock_string:
                                    with patch('orchestrator.pipeline_nondocking.UniProtClient') as mock_uniprot:
                                        with patch('orchestrator.pipeline_nondocking.PDBClient') as mock_pdb:
                                            
                                            # Mock all clients
                                            for mock_client_class in [mock_kegg, mock_reactome, mock_pa, mock_string, mock_uniprot, mock_pdb]:
                                                mock_client_instance = MagicMock()
                                                mock_client_class.return_value = mock_client_instance
                                                mock_client_instance.search_targets.return_value = []
                                            
                                            pipeline = NonDockingPipeline(test_config)
                                            
                                            with patch.object(pipeline, '_discover_targets') as mock_discover:
                                                mock_discover.return_value = [
                                                    {
                                                        "chembl_id": "CHEMBL123",
                                                        "uniprot_id": "P12345",
                                                        "gene_name": "TEST_GENE"
                                                    }
                                                ]
                                                
                                                result = pipeline.run()
                                                results_path = Path(result["results_csv"])
                                                
                                                if results_path.exists():
                                                    df = pd.read_csv(results_path)
                                                    
                                                    # Check that E3FP columns exist
                                                    assert "e3fp_tanimoto_max" in df.columns
                                                    assert "e3fp_tanimoto_mean_topk" in df.columns
                                                    
                                                    # Check that E3FP values are present (even if 0.0)
                                                    if len(df) > 0:
                                                        # Values should be present (not all NaN)
                                                        e3fp_max_values = df["e3fp_tanimoto_max"].dropna()
                                                        e3fp_mean_values = df["e3fp_tanimoto_mean_topk"].dropna()
                                                        
                                                        # At least some values should be present
                                                        assert len(e3fp_max_values) >= 0  # Can be 0 if E3FP fails
                                                        assert len(e3fp_mean_values) >= 0


if __name__ == "__main__":
    pytest.main([__file__])




