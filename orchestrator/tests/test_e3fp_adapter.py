"""
Test E3FP adapter functionality.
"""

import pytest
import numpy as np
from pathlib import Path
from unittest.mock import patch, MagicMock

from orchestrator.adapters.e3fp_adapter import (
    compute_e3fp_for_sdf,
    tanimoto_e3fp,
    topk_similarity_matrix_e3fp,
    compute_e3fp_features_for_ligands,
    _try_import_e3fp
)


class TestE3FPAdapter:
    """Test E3FP adapter functions."""
    
    def test_try_import_e3fp_no_deps(self):
        """Test that _try_import_e3fp returns None when E3FP is not available."""
        with patch('orchestrator.adapters.e3fp_adapter.e3fp', None):
            with patch('orchestrator.adapters.e3fp_adapter.mol_from_sdf', None):
                e3fp, mol_from_sdf, Fingerprint, Chem, AllChem = _try_import_e3fp()
                assert e3fp is None
                assert mol_from_sdf is None
                assert Fingerprint is None
                assert Chem is None
                assert AllChem is None
    
    def test_compute_e3fp_for_sdf_no_deps(self):
        """Test that compute_e3fp_for_sdf returns empty dict when E3FP is not available."""
        with patch('orchestrator.adapters.e3fp_adapter._try_import_e3fp') as mock_import:
            mock_import.return_value = (None, None, None, None, None)
            
            result = compute_e3fp_for_sdf(Path("test.sdf"))
            assert result == {}
    
    def test_compute_e3fp_for_sdf_missing_file(self):
        """Test that compute_e3fp_for_sdf returns empty dict when SDF file doesn't exist."""
        with patch('orchestrator.adapters.e3fp_adapter._try_import_e3fp') as mock_import:
            mock_import.return_value = (MagicMock(), MagicMock(), MagicMock(), MagicMock(), MagicMock())
            
            result = compute_e3fp_for_sdf(Path("nonexistent.sdf"))
            assert result == {}
    
    def test_tanimoto_e3fp(self):
        """Test Tanimoto similarity computation between E3FP bitvectors."""
        # Test identical fingerprints
        fp1 = np.array([1, 0, 1, 0, 1])
        fp2 = np.array([1, 0, 1, 0, 1])
        assert tanimoto_e3fp(fp1, fp2) == 1.0
        
        # Test completely different fingerprints
        fp1 = np.array([1, 0, 1, 0, 1])
        fp2 = np.array([0, 1, 0, 1, 0])
        assert tanimoto_e3fp(fp1, fp2) == 0.0
        
        # Test partial similarity
        fp1 = np.array([1, 0, 1, 0, 1])
        fp2 = np.array([1, 0, 0, 0, 1])
        expected = 2 / 3  # 2 bits in common, 3 bits in union
        assert abs(tanimoto_e3fp(fp1, fp2) - expected) < 1e-6
    
    def test_tanimoto_e3fp_different_shapes(self):
        """Test that Tanimoto returns 0.0 for different shaped fingerprints."""
        fp1 = np.array([1, 0, 1])
        fp2 = np.array([1, 0, 1, 0])
        assert tanimoto_e3fp(fp1, fp2) == 0.0
    
    def test_topk_similarity_matrix_e3fp_empty(self):
        """Test that topk_similarity_matrix_e3fp returns empty dict for empty inputs."""
        result = topk_similarity_matrix_e3fp({}, {}, k=5)
        assert result == {}
        
        result = topk_similarity_matrix_e3fp({"query1": np.array([1, 0, 1])}, {}, k=5)
        assert result == {}
    
    def test_topk_similarity_matrix_e3fp(self):
        """Test top-k similarity matrix computation."""
        query_fps = {
            "query1": np.array([1, 0, 1, 0, 1]),
            "query2": np.array([0, 1, 0, 1, 0])
        }
        ref_fps = {
            "ref1": np.array([1, 0, 1, 0, 1]),  # Identical to query1
            "ref2": np.array([0, 1, 0, 1, 0]),  # Identical to query2
            "ref3": np.array([1, 1, 0, 0, 0])   # Different
        }
        
        result = topk_similarity_matrix_e3fp(query_fps, ref_fps, k=2)
        
        assert "query1" in result
        assert "query2" in result
        
        # query1 should have ref1 as top match (similarity = 1.0)
        assert result["query1"]["e3fp_tanimoto_max"] == 1.0
        assert result["query1"]["neighbors"][0][0] == "ref1"
        assert result["query1"]["neighbors"][0][1] == 1.0
        
        # query2 should have ref2 as top match (similarity = 1.0)
        assert result["query2"]["e3fp_tanimoto_max"] == 1.0
        assert result["query2"]["neighbors"][0][0] == "ref2"
        assert result["query2"]["neighbors"][0][1] == 1.0
    
    def test_compute_e3fp_features_for_ligands_no_deps(self):
        """Test that compute_e3fp_features_for_ligands returns empty dict when E3FP is not available."""
        with patch('orchestrator.adapters.e3fp_adapter._try_import_e3fp') as mock_import:
            mock_import.return_value = (None, None, None, None, None)
            
            result = compute_e3fp_features_for_ligands(
                ["CCO"], ["ethanol"], radius=2, shells=5
            )
            assert result == {}
    
    @pytest.mark.skipif(
        _try_import_e3fp()[0] is None,
        reason="E3FP not installed"
    )
    def test_e3fp_integration_real(self):
        """Test E3FP integration with real library (if available)."""
        # This test only runs if E3FP is actually installed
        e3fp, mol_from_sdf, Fingerprint, Chem, AllChem = _try_import_e3fp()
        
        if e3fp is not None:
            # Test with a simple SMILES
            smiles_list = ["CCO"]  # ethanol
            names = ["ethanol"]
            
            result = compute_e3fp_features_for_ligands(
                smiles_list, names, radius=2, shells=5
            )
            
            # Should have computed fingerprints
            assert len(result) > 0
            assert "ethanol" in result
            
            # Fingerprint should be a numpy array
            fp = result["ethanol"]
            assert isinstance(fp, np.ndarray)
            assert fp.dtype == bool or fp.dtype == np.uint8 or fp.dtype == np.int8
            
            # Test similarity computation
            if len(result) >= 2:
                fps = list(result.values())
                sim = tanimoto_e3fp(fps[0], fps[1])
                assert 0.0 <= sim <= 1.0


if __name__ == "__main__":
    pytest.main([__file__])




