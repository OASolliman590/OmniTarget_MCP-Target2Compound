"""
Enhanced ChEMBL MCP Client with advanced pipeline integration.

This client leverages multiple ChEMBL MCP tools to provide comprehensive
drug discovery capabilities for the pipeline.
"""

from __future__ import annotations

import asyncio
import json
from typing import Any, Dict, List, Optional

from loguru import logger

from .chembl_mcp_client import ChEMBLMCPClient


class EnhancedChEMBLClient(ChEMBLMCPClient):
    """Enhanced ChEMBL client with advanced pipeline features."""
    
    async def discover_targets_by_disease(
        self,
        disease_terms: List[str],
        organism: str = "Homo sapiens",
        limit: int = 10
    ) -> List[Dict[str, Any]]:
        """Discover targets associated with disease terms."""
        discovered_targets = []
        
        for disease in disease_terms:
            try:
                # Search for targets related to the disease
                request = {
                    "jsonrpc": "2.0",
                    "id": 1,
                    "method": "tools/call",
                    "params": {
                        "name": "search_targets",
                        "arguments": {
                            "query": disease,
                            "organism": organism,
                            "limit": limit
                        }
                    }
                }
                
                response = await self._send_request(request)
                result = response.get("result", {}).get("content", [])
                
                if result:
                    data = json.loads(result[0]["text"])
                    if "targets" in data:
                        for target in data["targets"]:
                            discovered_targets.append({
                                "chembl_id": target.get("target_chembl_id"),
                                "uniprot_id": target.get("target_components", [{}])[0].get("accession"),
                                "gene_name": target.get("pref_name"),
                                "target_type": target.get("target_type"),
                                "organism": target.get("organism")
                            })
                            
            except Exception as e:
                logger.warning(f"Error discovering targets for {disease}: {e}")
        
        # Remove duplicates
        unique_targets = []
        seen_ids = set()
        for target in discovered_targets:
            if target["chembl_id"] not in seen_ids:
                unique_targets.append(target)
                seen_ids.add(target["chembl_id"])
        
        logger.info(f"Discovered {len(unique_targets)} unique targets for diseases: {disease_terms}")
        return unique_targets
    
    async def find_similar_compounds_enhanced(
        self,
        smiles_list: List[str],
        similarity_threshold: float = 0.7,
        limit_per_compound: int = 10
    ) -> Dict[str, List[Dict[str, Any]]]:
        """Find similar compounds for multiple SMILES using ChEMBL data."""
        similar_compounds = {}
        
        for smiles in smiles_list:
            try:
                request = {
                    "jsonrpc": "2.0",
                    "id": 1,
                    "method": "tools/call",
                    "params": {
                        "name": "search_similar_compounds",
                        "arguments": {
                            "smiles": smiles,
                            "similarity_threshold": similarity_threshold,
                            "limit": limit_per_compound
                        }
                    }
                }
                
                response = await self._send_request(request)
                result = response.get("result", {}).get("content", [])
                
                if result:
                    data = json.loads(result[0]["text"])
                    if "molecules" in data:
                        similar_compounds[smiles] = []
                        for molecule in data["molecules"]:
                            similar_compounds[smiles].append({
                                "chembl_id": molecule.get("molecule_chembl_id"),
                                "smiles": molecule.get("molecule_structures", {}).get("canonical_smiles"),
                                "name": molecule.get("pref_name"),
                                "molecular_weight": molecule.get("molecule_properties", {}).get("full_mwt"),
                                "logp": molecule.get("molecule_properties", {}).get("alogp"),
                                "similarity_score": molecule.get("similarity_score", 0.0)
                            })
                
            except Exception as e:
                logger.warning(f"Error finding similar compounds for {smiles}: {e}")
                similar_compounds[smiles] = []
        
        total_found = sum(len(compounds) for compounds in similar_compounds.values())
        logger.info(f"Found {total_found} similar compounds across {len(smiles_list)} input structures")
        return similar_compounds
    
    async def get_target_compounds_enhanced(
        self,
        target_chembl_id: str,
        activity_type: str = "IC50",
        min_pchembl: float = 6.0,
        limit: int = 50
    ) -> List[Dict[str, Any]]:
        """Get compounds tested against a target with enhanced filtering."""
        try:
            # First try the specific tool
            request = {
                "jsonrpc": "2.0",
                "id": 1,
                "method": "tools/call",
                "params": {
                    "name": "get_target_compounds",
                    "arguments": {
                        "target_chembl_id": target_chembl_id,
                        "activity_type": activity_type,
                        "limit": limit
                    }
                }
            }
            
            response = await self._send_request(request)
            result = response.get("result", {}).get("content", [])
            
            if result and "not yet implemented" not in result[0].get("text", ""):
                # Parse the response
                data = json.loads(result[0]["text"])
                compounds = []
                
                if "compounds" in data:
                    for compound in data["compounds"]:
                        pchembl_value = compound.get("pchembl_value")
                        if pchembl_value and float(pchembl_value) >= min_pchembl:
                            compounds.append({
                                "chembl_id": compound.get("molecule_chembl_id"),
                                "smiles": compound.get("canonical_smiles"),
                                "pchembl_value": float(pchembl_value),
                                "assay_type": compound.get("assay_type", activity_type),
                                "organism": compound.get("target_organism", "Homo sapiens")
                            })
                
                return compounds
            
            # Fallback to search_activities if get_target_compounds not implemented
            logger.info(f"Using search_activities as fallback for {target_chembl_id}")
            return await self.search_activities(
                target_chembl_id=target_chembl_id,
                activity_type=activity_type,
                limit=limit
            )
            
        except Exception as e:
            logger.error(f"Error getting target compounds for {target_chembl_id}: {e}")
            return []
    
    async def batch_compound_analysis(
        self,
        chembl_ids: List[str]
    ) -> Dict[str, Dict[str, Any]]:
        """Analyze multiple compounds in batch for comprehensive properties."""
        try:
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
            result = response.get("result", {}).get("content", [])
            
            if result:
                data = json.loads(result[0]["text"])
                compound_analysis = {}
                
                if "compounds" in data:
                    for compound in data["compounds"]:
                        chembl_id = compound.get("molecule_chembl_id")
                        compound_analysis[chembl_id] = {
                            "smiles": compound.get("molecule_structures", {}).get("canonical_smiles"),
                            "name": compound.get("pref_name"),
                            "molecular_weight": compound.get("molecule_properties", {}).get("full_mwt"),
                            "logp": compound.get("molecule_properties", {}).get("alogp"),
                            "hbd": compound.get("molecule_properties", {}).get("hbd"),
                            "hba": compound.get("molecule_properties", {}).get("hba"),
                            "psa": compound.get("molecule_properties", {}).get("psa"),
                            "lipinski_violations": compound.get("molecule_properties", {}).get("num_lipinski_ro5_violations"),
                            "drug_like": compound.get("molecule_properties", {}).get("num_lipinski_ro5_violations", 0) == 0
                        }
                
                logger.info(f"Analyzed {len(compound_analysis)} compounds in batch")
                return compound_analysis
                
        except Exception as e:
            logger.error(f"Error in batch compound analysis: {e}")
            return {}


# Convenience functions for pipeline integration
async def discover_disease_targets(
    disease_terms: List[str],
    organism: str = "Homo sapiens"
) -> List[Dict[str, Any]]:
    """Discover targets associated with disease terms using ChEMBL."""
    client = EnhancedChEMBLClient()
    
    try:
        await client.start_server()
        targets = await client.discover_targets_by_disease(disease_terms, organism)
        return targets
    except Exception as e:
        logger.error(f"Error discovering disease targets: {e}")
        return []
    finally:
        await client.stop_server()


async def find_chembl_similar_compounds(
    smiles_list: List[str],
    similarity_threshold: float = 0.7
) -> Dict[str, List[Dict[str, Any]]]:
    """Find ChEMBL compounds similar to input SMILES."""
    client = EnhancedChEMBLClient()
    
    try:
        await client.start_server()
        similar_compounds = await client.find_similar_compounds_enhanced(
            smiles_list, similarity_threshold
        )
        return similar_compounds
    except Exception as e:
        logger.error(f"Error finding similar compounds: {e}")
        return {}
    finally:
        await client.stop_server()
