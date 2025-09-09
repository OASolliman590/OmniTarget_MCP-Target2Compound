"""
Reactome MCP client for pathway enrichment analysis.
"""

from typing import List, Dict, Any, Optional
from loguru import logger

from .base import BaseMCPClient
from ..settings import settings


class ReactomeClient(BaseMCPClient):
    """Client for Reactome MCP server."""
    
    def __init__(self):
        super().__init__(
            base_url=settings.mcp.reactome_url,
            timeout=settings.mcp.request_timeout,
            max_retries=settings.mcp.max_retries,
            rate_limit=settings.mcp.rate_limit
        )
    
    async def test_connection(self) -> bool:
        """Test Reactome MCP connection."""
        try:
            response = await self.get("species")
            return isinstance(response, (dict, list))
        except Exception as e:
            logger.warning(f"Reactome connection test failed: {e}")
            return False
    
    async def pathway_enrichment(
        self,
        gene_list: List[str],
        species: str = "Homo sapiens",
        p_value_cutoff: float = 0.05
    ) -> List[Dict[str, Any]]:
        """Perform pathway enrichment analysis.
        
        Args:
            gene_list: List of gene symbols
            species: Species name
            p_value_cutoff: P-value threshold
            
        Returns:
            List of enriched pathway records
        """
        try:
            logger.info(f"Running Reactome enrichment for {len(gene_list)} genes")
            
            response = await self.post(
                "analysis/identifiers",
                data={
                    "identifiers": gene_list,
                    "species": species,
                    "includeDisease": True,
                    "projectToHuman": True
                }
            )
            
            pathways = []
            if "pathways" in response:
                for pathway in response["pathways"]:
                    if pathway.get("pValue", 1.0) <= p_value_cutoff:
                        pathway_data = {
                            "pathway_id": pathway.get("stId"),
                            "pathway_name": pathway.get("name"),
                            "p_value": pathway.get("pValue"),
                            "fdr": pathway.get("fdr"),
                            "found_entities": pathway.get("found"),
                            "total_entities": pathway.get("total"),
                            "species": species,
                            "source": "Reactome"
                        }
                        pathways.append(pathway_data)
            
            logger.info(f"Found {len(pathways)} enriched pathways")
            return pathways
            
        except Exception as e:
            logger.error(f"Error in Reactome enrichment: {e}")
            return []
    
    async def search_pathways(
        self,
        query: str,
        species: str = "Homo sapiens"
    ) -> List[Dict[str, Any]]:
        """Search pathways by term.
        
        Args:
            query: Search query
            species: Species name
            
        Returns:
            List of pathway records
        """
        try:
            response = await self.get(
                "search/query",
                params={
                    "query": query,
                    "species": species,
                    "types": "Pathway"
                }
            )
            
            pathways = []
            if "results" in response:
                for result in response["results"]:
                    pathway_data = {
                        "pathway_id": result.get("stId"),
                        "pathway_name": result.get("name"),
                        "description": result.get("summation"),
                        "species": species,
                        "source": "Reactome"
                    }
                    pathways.append(pathway_data)
            
            return pathways
            
        except Exception as e:
            logger.error(f"Error searching Reactome pathways: {e}")
            return []
