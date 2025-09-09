"""ChEMBL MCP client."""
from typing import List, Dict, Any
from loguru import logger
from .base import BaseMCPClient
from ..settings import settings

class ChEMBLClient(BaseMCPClient):
    def __init__(self):
        super().__init__(base_url=settings.mcp.chembl_url)
    
    async def test_connection(self) -> bool:
        try:
            response = await self.get("status")
            return "chembl" in str(response).lower()
        except:
            return False
    
    async def get_bioactivity_data(
        self, 
        target_id: str, 
        compound_id: str = None
    ) -> List[Dict[str, Any]]:
        try:
            params = {"target": target_id}
            if compound_id:
                params["compound"] = compound_id
            
            response = await self.get("activities", params=params)
            return response.get("activities", [])
        except Exception as e:
            logger.error(f"Error getting ChEMBL bioactivity data: {e}")
            return []
