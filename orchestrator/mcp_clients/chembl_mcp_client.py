"""
Real ChEMBL MCP Server client for advanced drug discovery features.

This client integrates with the comprehensive ChEMBL MCP Server that provides
22 specialized tools for drug discovery, chemical informatics, and bioactivity analysis.
"""

from __future__ import annotations

import asyncio
import json
import subprocess
from typing import Any, Dict, List, Optional

from loguru import logger


class ChEMBLMCPClient:
    """Client for the real ChEMBL MCP Server with 22 specialized tools."""
    
    def __init__(self, server_path: str = "services/chembl-mcp/build/index.js"):
        self.server_path = server_path
        self.process: Optional[subprocess.Popen] = None
    
    async def start_server(self) -> None:
        """Start the ChEMBL MCP Server process."""
        try:
            self.process = subprocess.Popen(
                ["node", self.server_path],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            logger.info("ChEMBL MCP Server started")
        except Exception as e:
            logger.error(f"Failed to start ChEMBL MCP Server: {e}")
            raise
    
    async def stop_server(self) -> None:
        """Stop the ChEMBL MCP Server process."""
        if self.process:
            self.process.terminate()
            self.process.wait()
            self.process = None
            logger.info("ChEMBL MCP Server stopped")
    
    async def _send_request(self, request: Dict[str, Any]) -> Dict[str, Any]:
        """Send a request to the ChEMBL MCP Server."""
        if not self.process:
            await self.start_server()
        
        try:
            # Send request
            request_json = json.dumps(request) + "\n"
            self.process.stdin.write(request_json)
            self.process.stdin.flush()
            
            # Read response
            response_line = self.process.stdout.readline()
            if response_line:
                return json.loads(response_line.strip())
            else:
                raise Exception("No response from ChEMBL MCP Server")
                
        except Exception as e:
            logger.error(f"Error communicating with ChEMBL MCP Server: {e}")
            raise
    
    async def search_activities(
        self,
        target_chembl_id: Optional[str] = None,
        activity_type: Optional[str] = None,
        limit: int = 25
    ) -> List[Dict[str, Any]]:
        """Search bioactivity measurements and assay results."""
        request = {
            "jsonrpc": "2.0",
            "id": 1,
            "method": "tools/call",
            "params": {
                "name": "search_activities",
                "arguments": {
                    "limit": limit
                }
            }
        }
        
        if target_chembl_id:
            request["params"]["arguments"]["target_chembl_id"] = target_chembl_id
        if activity_type:
            request["params"]["arguments"]["activity_type"] = activity_type
        
        response = await self._send_request(request)
        
        if "result" in response and "content" in response["result"]:
            return response["result"]["content"]
        else:
            logger.warning(f"No activities found for target {target_chembl_id}")
            return []
    
    async def get_target_compounds(
        self,
        target_chembl_id: str,
        activity_type: Optional[str] = None,
        limit: int = 25
    ) -> List[Dict[str, Any]]:
        """Get compounds tested against a specific target."""
        request = {
            "jsonrpc": "2.0",
            "id": 1,
            "method": "tools/call",
            "params": {
                "name": "get_target_compounds",
                "arguments": {
                    "target_chembl_id": target_chembl_id,
                    "limit": limit
                }
            }
        }
        
        if activity_type:
            request["params"]["arguments"]["activity_type"] = activity_type
        
        response = await self._send_request(request)
        
        if "result" in response and "content" in response["result"]:
            return response["result"]["content"]
        else:
            logger.warning(f"No compounds found for target {target_chembl_id}")
            return []
    
    async def search_similar_compounds(
        self,
        smiles: str,
        similarity_threshold: float = 0.7,
        limit: int = 25
    ) -> List[Dict[str, Any]]:
        """Find chemically similar compounds using Tanimoto similarity."""
        request = {
            "jsonrpc": "2.0",
            "id": 1,
            "method": "tools/call",
            "params": {
                "name": "search_similar_compounds",
                "arguments": {
                    "smiles": smiles,
                    "similarity_threshold": similarity_threshold,
                    "limit": limit
                }
            }
        }
        
        response = await self._send_request(request)
        
        if "result" in response and "content" in response["result"]:
            return response["result"]["content"]
        else:
            logger.warning(f"No similar compounds found for SMILES: {smiles}")
            return []
    
    async def batch_compound_lookup(
        self,
        chembl_ids: List[str]
    ) -> List[Dict[str, Any]]:
        """Process multiple ChEMBL IDs efficiently."""
        request = {
            "jsonrpc": "2.0",
            "id": 1,
            "method": "tools/call",
            "params": {
                "name": "batch_compound_lookup",
                "arguments": {
                    "chembl_ids": chembl_ids
                }
            }
        }
        
        response = await self._send_request(request)
        
        if "result" in response and "content" in response["result"]:
            return response["result"]["content"]
        else:
            logger.warning(f"No compounds found for IDs: {chembl_ids}")
            return []
    
    async def get_compound_info(
        self,
        chembl_id: str
    ) -> Optional[Dict[str, Any]]:
        """Get detailed compound information by ChEMBL ID."""
        request = {
            "jsonrpc": "2.0",
            "id": 1,
            "method": "tools/call",
            "params": {
                "name": "get_compound_info",
                "arguments": {
                    "chembl_id": chembl_id
                }
            }
        }
        
        response = await self._send_request(request)
        
        if "result" in response and "content" in response["result"]:
            return response["result"]["content"]
        else:
            logger.warning(f"No compound info found for ID: {chembl_id}")
            return None


# Convenience functions for backward compatibility
async def fetch_target_actives_real(
    chembl_target_id: str, 
    min_pchembl: float = 6.0,
    activity_type: str = "IC50"
) -> List[Dict]:
    """Fetch active compounds for a ChEMBL target using real ChEMBL data."""
    client = ChEMBLMCPClient()
    
    try:
        # Get compounds tested against the target
        compounds = await client.get_target_compounds(
            target_chembl_id=chembl_target_id,
            activity_type=activity_type,
            limit=100
        )
        
        # Filter by pChEMBL value
        filtered_compounds = []
        for compound in compounds:
            pchembl_value = compound.get("pchembl_value")
            if pchembl_value and float(pchembl_value) >= min_pchembl:
                filtered_compounds.append({
                    "chembl_id": compound.get("molecule_chembl_id"),
                    "smiles": compound.get("canonical_smiles"),
                    "pchembl_value": float(pchembl_value),
                    "assay_type": compound.get("assay_type", activity_type),
                    "organism": compound.get("target_organism", "Homo sapiens")
                })
        
        logger.info(f"Found {len(filtered_compounds)} active compounds for {chembl_target_id}")
        return filtered_compounds
        
    except Exception as e:
        logger.error(f"Error fetching real ChEMBL data: {e}")
        return []
    finally:
        await client.stop_server()


async def search_similar_compounds_real(
    smiles: str,
    similarity_threshold: float = 0.7,
    limit: int = 25
) -> List[Dict]:
    """Find compounds similar to the given SMILES using real ChEMBL data."""
    client = ChEMBLMCPClient()
    
    try:
        similar_compounds = await client.search_similar_compounds(
            smiles=smiles,
            similarity_threshold=similarity_threshold,
            limit=limit
        )
        
        # Convert to our format
        formatted_compounds = []
        for compound in similar_compounds:
            formatted_compounds.append({
                "chembl_id": compound.get("molecule_chembl_id"),
                "smiles": compound.get("canonical_smiles"),
                "similarity": compound.get("similarity_score"),
                "name": compound.get("pref_name", ""),
                "molecular_weight": compound.get("molecular_weight")
            })
        
        logger.info(f"Found {len(formatted_compounds)} similar compounds for SMILES: {smiles}")
        return formatted_compounds
        
    except Exception as e:
        logger.error(f"Error searching similar compounds: {e}")
        return []
    finally:
        await client.stop_server()
