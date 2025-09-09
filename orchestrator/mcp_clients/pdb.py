"""PDB MCP client."""
from typing import List, Dict, Any, Optional
from loguru import logger
from .base import BaseMCPClient
from ..settings import settings

class PDBClient(BaseMCPClient):
    def __init__(self):
        super().__init__(base_url=settings.mcp.pdb_url)
    
    async def test_connection(self) -> bool:
        try:
            response = await self.get("status")
            return "pdb" in str(response).lower()
        except:
            return False
    
    async def search_structures(self, uniprot_id: str) -> List[Dict[str, Any]]:
        try:
            response = await self.get(f"search/uniprot/{uniprot_id}")
            return response.get("structures", [])
        except Exception as e:
            logger.error(f"Error searching PDB structures: {e}")
            return []
