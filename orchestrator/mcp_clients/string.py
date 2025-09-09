"""
STRING MCP client for protein-protein interaction data.
"""

from typing import List, Dict, Any, Optional
from loguru import logger
from .base import BaseMCPClient
from ..settings import settings


class STRINGClient(BaseMCPClient):
    """Client for STRING MCP server."""
    
    def __init__(self):
        super().__init__(base_url=settings.mcp.string_url)
    
    async def test_connection(self) -> bool:
        try:
            response = await self.get("version")
            return "string" in str(response).lower()
        except:
            return False
    
    async def get_protein_interactions(
        self,
        proteins: List[str],
        species: int = 9606,
        required_score: int = 700
    ) -> List[Dict[str, Any]]:
        """Get protein-protein interactions."""
        try:
            response = await self.post(
                "network",
                data={
                    "identifiers": proteins,
                    "species": species,
                    "required_score": required_score
                }
            )
            
            interactions = []
            if "interactions" in response:
                for interaction in response["interactions"]:
                    interactions.append({
                        "protein_a": interaction.get("stringId_A"),
                        "protein_b": interaction.get("stringId_B"),
                        "combined_score": interaction.get("score"),
                        "source": "STRING"
                    })
            
            return interactions
        except Exception as e:
            logger.error(f"Error getting STRING interactions: {e}")
            return []
