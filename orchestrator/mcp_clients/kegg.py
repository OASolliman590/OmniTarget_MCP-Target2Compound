"""
KEGG MCP client for pathway and gene data.
"""

from typing import List, Dict, Any, Optional
from loguru import logger

from .base import BaseMCPClient, MCPError
from ..settings import settings


class KEGGClient(BaseMCPClient):
    """Client for KEGG MCP server."""
    
    def __init__(self):
        super().__init__(
            base_url=settings.mcp.kegg_url,
            timeout=settings.mcp.request_timeout,
            max_retries=settings.mcp.max_retries,
            rate_limit=settings.mcp.rate_limit
        )
    
    async def test_connection(self) -> bool:
        """Test KEGG MCP connection."""
        try:
            response = await self.get("info")
            return "kegg" in str(response).lower()
        except Exception as e:
            logger.warning(f"KEGG connection test failed: {e}")
            return False
    
    async def search_disease_pathways(
        self,
        disease_terms: List[str],
        organism: str = "hsa"
    ) -> List[Dict[str, Any]]:
        """Search for pathways associated with disease terms.
        
        Args:
            disease_terms: List of disease terms to search
            organism: KEGG organism code (default: hsa for human)
            
        Returns:
            List of pathway records
        """
        pathways = []
        
        for term in disease_terms:
            try:
                logger.info(f"Searching KEGG pathways for disease: {term}")
                
                # Search for disease-related pathways
                response = await self.get(
                    "search/pathway",
                    params={
                        "query": term,
                        "organism": organism
                    }
                )
                
                if "pathways" in response:
                    for pathway in response["pathways"]:
                        pathway_data = {
                            "pathway_id": pathway.get("id"),
                            "pathway_name": pathway.get("name"),
                            "description": pathway.get("description"),
                            "organism": organism,
                            "disease_term": term,
                            "source": "KEGG"
                        }
                        pathways.append(pathway_data)
                        
                logger.info(f"Found {len(response.get('pathways', []))} pathways for {term}")
                
            except Exception as e:
                logger.error(f"Error searching KEGG pathways for {term}: {e}")
                continue
        
        return pathways
    
    async def get_pathway_genes(
        self,
        pathway_id: str,
        organism: str = "hsa"
    ) -> List[Dict[str, Any]]:
        """Get genes associated with a pathway.
        
        Args:
            pathway_id: KEGG pathway identifier
            organism: KEGG organism code
            
        Returns:
            List of gene records
        """
        try:
            logger.debug(f"Getting genes for pathway {pathway_id}")
            
            response = await self.get(
                f"pathway/{pathway_id}/genes",
                params={"organism": organism}
            )
            
            genes = []
            if "genes" in response:
                for gene in response["genes"]:
                    gene_data = {
                        "gene_id": gene.get("id"),
                        "gene_symbol": gene.get("symbol"),
                        "gene_name": gene.get("name"),
                        "description": gene.get("description"),
                        "pathway_id": pathway_id,
                        "organism": organism,
                        "source": "KEGG"
                    }
                    genes.append(gene_data)
            
            logger.debug(f"Found {len(genes)} genes in pathway {pathway_id}")
            return genes
            
        except Exception as e:
            logger.error(f"Error getting genes for pathway {pathway_id}: {e}")
            return []
    
    async def get_gene_info(
        self,
        gene_id: str,
        organism: str = "hsa"
    ) -> Optional[Dict[str, Any]]:
        """Get detailed information for a gene.
        
        Args:
            gene_id: KEGG gene identifier
            organism: KEGG organism code
            
        Returns:
            Gene information dictionary or None
        """
        try:
            logger.debug(f"Getting info for gene {gene_id}")
            
            response = await self.get(
                f"gene/{gene_id}",
                params={"organism": organism}
            )
            
            if "gene" in response:
                gene_data = response["gene"]
                return {
                    "gene_id": gene_data.get("id"),
                    "gene_symbol": gene_data.get("symbol"),
                    "gene_name": gene_data.get("name"),
                    "description": gene_data.get("description"),
                    "position": gene_data.get("position"),
                    "organism": organism,
                    "pathways": gene_data.get("pathways", []),
                    "orthologs": gene_data.get("orthologs", []),
                    "source": "KEGG"
                }
            
            return None
            
        except Exception as e:
            logger.error(f"Error getting info for gene {gene_id}: {e}")
            return None
    
    async def convert_gene_ids(
        self,
        gene_ids: List[str],
        from_db: str = "kegg",
        to_db: str = "uniprot",
        organism: str = "hsa"
    ) -> Dict[str, str]:
        """Convert gene IDs between different databases.
        
        Args:
            gene_ids: List of gene identifiers
            from_db: Source database
            to_db: Target database
            organism: KEGG organism code
            
        Returns:
            Dictionary mapping input IDs to converted IDs
        """
        mapping = {}
        
        try:
            logger.info(f"Converting {len(gene_ids)} gene IDs from {from_db} to {to_db}")
            
            response = await self.post(
                "convert",
                data={
                    "ids": gene_ids,
                    "from": from_db,
                    "to": to_db,
                    "organism": organism
                }
            )
            
            if "mapping" in response:
                mapping = response["mapping"]
                logger.info(f"Successfully converted {len(mapping)} gene IDs")
            
        except Exception as e:
            logger.error(f"Error converting gene IDs: {e}")
        
        return mapping
    
    async def search_compounds(
        self,
        compound_name: str
    ) -> List[Dict[str, Any]]:
        """Search for compounds in KEGG.
        
        Args:
            compound_name: Compound name to search
            
        Returns:
            List of compound records
        """
        try:
            logger.debug(f"Searching KEGG compounds for: {compound_name}")
            
            response = await self.get(
                "search/compound",
                params={"query": compound_name}
            )
            
            compounds = []
            if "compounds" in response:
                for compound in response["compounds"]:
                    compound_data = {
                        "compound_id": compound.get("id"),
                        "name": compound.get("name"),
                        "formula": compound.get("formula"),
                        "molecular_weight": compound.get("mol_weight"),
                        "pathways": compound.get("pathways", []),
                        "source": "KEGG"
                    }
                    compounds.append(compound_data)
            
            logger.debug(f"Found {len(compounds)} compounds for {compound_name}")
            return compounds
            
        except Exception as e:
            logger.error(f"Error searching compounds for {compound_name}: {e}")
            return []
    
    async def get_organism_list(self) -> List[Dict[str, Any]]:
        """Get list of available organisms in KEGG.
        
        Returns:
            List of organism records
        """
        try:
            response = await self.get("organisms")
            
            organisms = []
            if "organisms" in response:
                for org in response["organisms"]:
                    organism_data = {
                        "code": org.get("code"),
                        "name": org.get("name"),
                        "scientific_name": org.get("scientific_name"),
                        "taxonomy_id": org.get("taxonomy_id")
                    }
                    organisms.append(organism_data)
            
            return organisms
            
        except Exception as e:
            logger.error(f"Error getting organism list: {e}")
            return []
    
    async def discover_disease_targets(
        self,
        disease_terms: List[str],
        organism: str = "hsa",
        max_pathways_per_disease: int = 50
    ) -> List[Dict[str, Any]]:
        """High-level method to discover targets for diseases.
        
        Args:
            disease_terms: List of disease terms
            organism: KEGG organism code
            max_pathways_per_disease: Maximum pathways to process per disease
            
        Returns:
            List of target gene records with pathway context
        """
        targets = []
        
        # Search for pathways related to each disease
        pathways = await self.search_disease_pathways(disease_terms, organism)
        
        # Limit pathways per disease to avoid excessive API calls
        if len(pathways) > max_pathways_per_disease:
            logger.warning(f"Too many pathways found ({len(pathways)}), limiting to {max_pathways_per_disease}")
            pathways = pathways[:max_pathways_per_disease]
        
        # Get genes for each pathway
        for pathway in pathways:
            pathway_id = pathway["pathway_id"]
            
            genes = await self.get_pathway_genes(pathway_id, organism)
            
            for gene in genes:
                # Add pathway context to gene
                target_data = {
                    **gene,
                    "discovery_pathway": pathway,
                    "discovery_source": "KEGG",
                    "target_id": gene["gene_id"]
                }
                targets.append(target_data)
        
        logger.info(f"Discovered {len(targets)} targets from {len(pathways)} pathways")
        return targets
