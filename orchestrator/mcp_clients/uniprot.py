"""
UniProt MCP client for protein sequence and annotation data.
"""

from typing import List, Dict, Any, Optional
from loguru import logger

from .base import BaseMCPClient
from ..settings import settings


class UniProtClient(BaseMCPClient):
    """Client for UniProt MCP server."""
    
    def __init__(self):
        super().__init__(
            base_url=settings.mcp.uniprot_url,
            timeout=settings.mcp.request_timeout,
            max_retries=settings.mcp.max_retries,
            rate_limit=settings.mcp.rate_limit
        )
    
    async def test_connection(self) -> bool:
        """Test UniProt MCP connection."""
        try:
            response = await self.get("status")
            return "uniprot" in str(response).lower()
        except Exception as e:
            logger.warning(f"UniProt connection test failed: {e}")
            return False
    
    async def get_protein_info(
        self,
        accession: str
    ) -> Optional[Dict[str, Any]]:
        """Get protein information by UniProt accession.
        
        Args:
            accession: UniProt accession ID
            
        Returns:
            Protein information dictionary or None
        """
        try:
            response = await self.get(f"proteins/{accession}")
            
            if "protein" in response:
                protein = response["protein"]
                return {
                    "accession": protein.get("accession"),
                    "gene_name": protein.get("geneName"),
                    "protein_name": protein.get("proteinName"),
                    "organism": protein.get("organism"),
                    "sequence": protein.get("sequence"),
                    "length": protein.get("length"),
                    "mass": protein.get("mass"),
                    "function": protein.get("function"),
                    "go_terms": protein.get("goTerms", []),
                    "active_sites": protein.get("activeSites", []),
                    "binding_sites": protein.get("bindingSites", []),
                    "domains": protein.get("domains", []),
                    "source": "UniProt"
                }
            
            return None
            
        except Exception as e:
            logger.error(f"Error getting protein info for {accession}: {e}")
            return None
    
    async def search_proteins(
        self,
        query: str,
        organism: str = "human",
        limit: int = 100
    ) -> List[Dict[str, Any]]:
        """Search proteins by query.
        
        Args:
            query: Search query
            organism: Organism filter
            limit: Maximum results
            
        Returns:
            List of protein records
        """
        try:
            response = await self.get(
                "search",
                params={
                    "query": query,
                    "organism": organism,
                    "limit": limit
                }
            )
            
            proteins = []
            if "results" in response:
                for protein in response["results"]:
                    protein_data = {
                        "accession": protein.get("accession"),
                        "gene_name": protein.get("geneName"),
                        "protein_name": protein.get("proteinName"),
                        "organism": protein.get("organism"),
                        "score": protein.get("score"),
                        "source": "UniProt"
                    }
                    proteins.append(protein_data)
            
            return proteins
            
        except Exception as e:
            logger.error(f"Error searching proteins: {e}")
            return []
