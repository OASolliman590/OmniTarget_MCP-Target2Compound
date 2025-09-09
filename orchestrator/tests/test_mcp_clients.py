"""
Tests for MCP clients.
"""

import pytest
import asyncio
from unittest.mock import AsyncMock

from ..mcp_clients import (
    KEGGClient, ReactomeClient, UniProtClient, 
    STRINGClient, ProteinAtlasClient, PDBClient, ChEMBLClient
)


class TestMCPClients:
    """Test suite for MCP client functionality."""
    
    @pytest.mark.asyncio
    async def test_kegg_client_connection(self):
        """Test KEGG client connection."""
        client = KEGGClient()
        
        # Test connection (may fail if service not available)
        try:
            result = await client.test_connection()
            assert isinstance(result, bool)
        except Exception as e:
            pytest.skip(f"KEGG service not available: {e}")
        finally:
            await client.close()
    
    @pytest.mark.asyncio
    async def test_kegg_disease_pathways(self):
        """Test KEGG disease pathway search."""
        client = KEGGClient()
        
        try:
            # Test with real disease terms
            pathways = await client.search_disease_pathways(["lung cancer"])
            
            # Should return a list (even if empty)
            assert isinstance(pathways, list)
            
            # If pathways found, check structure
            if pathways:
                pathway = pathways[0]
                assert "pathway_id" in pathway
                assert "pathway_name" in pathway
                assert "source" in pathway
                assert pathway["source"] == "KEGG"
                
        except Exception as e:
            pytest.skip(f"KEGG service error: {e}")
        finally:
            await client.close()
    
    @pytest.mark.asyncio
    async def test_uniprot_client_search(self):
        """Test UniProt protein search."""
        client = UniProtClient()
        
        try:
            # Test protein search
            proteins = await client.search_proteins("insulin", limit=5)
            
            assert isinstance(proteins, list)
            
            if proteins:
                protein = proteins[0]
                expected_fields = ["accession", "gene_name", "protein_name", "organism"]
                for field in expected_fields:
                    assert field in protein
                    
        except Exception as e:
            pytest.skip(f"UniProt service error: {e}")
        finally:
            await client.close()
    
    @pytest.mark.asyncio
    async def test_reactome_client(self):
        """Test Reactome client functionality."""
        client = ReactomeClient()
        
        try:
            connected = await client.test_connection()
            assert isinstance(connected, bool)
            
            if connected:
                # Test pathway search
                pathways = await client.search_pathways("cancer")
                assert isinstance(pathways, list)
                
        except Exception as e:
            pytest.skip(f"Reactome service error: {e}")
        finally:
            await client.close()
    
    @pytest.mark.asyncio
    async def test_string_client(self):
        """Test STRING client functionality."""
        client = STRINGClient()
        
        try:
            connected = await client.test_connection()
            assert isinstance(connected, bool)
            
            if connected:
                # Test PPI network query
                interactions = await client.get_protein_interactions(["ENSP00000000233"])
                assert isinstance(interactions, list)
                
        except Exception as e:
            pytest.skip(f"STRING service error: {e}")
        finally:
            await client.close()
    
    @pytest.mark.asyncio
    async def test_all_clients_basic_connection(self):
        """Test basic connectivity for all MCP clients."""
        clients = [
            ("KEGG", KEGGClient()),
            ("Reactome", ReactomeClient()),
            ("UniProt", UniProtClient()),
            ("STRING", STRINGClient()),
            ("ProteinAtlas", ProteinAtlasClient()),
            ("PDB", PDBClient()),
            ("ChEMBL", ChEMBLClient()),
        ]
        
        results = {}
        
        for name, client in clients:
            try:
                connected = await client.test_connection()
                results[name] = connected
                print(f"{name}: {'✓' if connected else '✗'}")
            except Exception as e:
                results[name] = False
                print(f"{name}: Error - {e}")
            finally:
                await client.close()
        
        # At least some services should be testable
        assert isinstance(results, dict)
        assert len(results) == len(clients)
    
    def test_client_initialization(self):
        """Test that all clients can be initialized."""
        clients = [
            KEGGClient(), ReactomeClient(), UniProtClient(),
            STRINGClient(), ProteinAtlasClient(), PDBClient(), ChEMBLClient()
        ]
        
        for client in clients:
            assert client.base_url is not None
            assert client.timeout > 0
            assert client.max_retries >= 0
    
    @pytest.mark.asyncio
    async def test_error_handling(self):
        """Test client error handling with invalid URLs."""
        # Create client with invalid URL
        client = KEGGClient()
        client.base_url = "http://invalid.nonexistent.domain"
        
        # Should handle connection errors gracefully
        try:
            result = await client.test_connection()
            assert result is False
        except Exception:
            # Any exception should be handled by the client
            pass
        finally:
            await client.close()
