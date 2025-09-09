"""Protein Atlas MCP client."""
from typing import List, Dict, Any
from loguru import logger
from .base import BaseMCPClient
from ..settings import settings

class ProteinAtlasClient(BaseMCPClient):
    def __init__(self):
        super().__init__(base_url=settings.mcp.proteinatlas_url)
    
    async def test_connection(self) -> bool:
        try:
            response = await self.get("genes")
            return isinstance(response, (dict, list))
        except:
            return False
    
    async def get_expression_data(self, gene_name: str) -> List[Dict[str, Any]]:
        try:
            response = await self.get(f"genes/{gene_name}/expression")
            return response.get("expression", [])
        except Exception as e:
            logger.error(f"Error getting expression data: {e}")
            return []
