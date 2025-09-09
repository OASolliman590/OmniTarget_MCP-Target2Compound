"""
Target and biological data schemas.
"""

from typing import List, Optional, Dict, Any
from pydantic import BaseModel, Field
from datetime import datetime


class PathwayHit(BaseModel):
    """A pathway hit from KEGG or Reactome."""
    
    pathway_id: str = Field(description="Pathway identifier")
    pathway_name: str = Field(description="Pathway name/description")
    source: str = Field(description="Data source (KEGG/Reactome)")
    genes: List[str] = Field(description="Gene identifiers in pathway")
    p_value: Optional[float] = Field(default=None, description="Enrichment p-value")
    fdr: Optional[float] = Field(default=None, description="False discovery rate")
    
    class Config:
        extra = "allow"


class ExpressionRecord(BaseModel):
    """Protein expression data from Human Protein Atlas."""
    
    gene_name: str = Field(description="Gene symbol")
    tissue: str = Field(description="Tissue/cell type")
    cell_type: Optional[str] = Field(default=None, description="Specific cell type")
    expression_level: str = Field(description="Expression level (High/Medium/Low/Not detected)")
    reliability: str = Field(description="Reliability score")
    intensity: Optional[float] = Field(default=None, description="Staining intensity")
    quantity: Optional[float] = Field(default=None, description="Quantity score")
    location: Optional[str] = Field(default=None, description="Subcellular location")
    
    class Config:
        extra = "allow"


class PPIRecord(BaseModel):
    """Protein-protein interaction record from STRING."""
    
    protein_a: str = Field(description="First protein identifier")
    protein_b: str = Field(description="Second protein identifier")
    combined_score: float = Field(description="Combined confidence score")
    experimental_score: Optional[float] = Field(default=None, description="Experimental evidence")
    database_score: Optional[float] = Field(default=None, description="Database evidence")
    textmining_score: Optional[float] = Field(default=None, description="Text mining evidence")
    coexpression_score: Optional[float] = Field(default=None, description="Co-expression evidence")
    neighborhood_score: Optional[float] = Field(default=None, description="Neighborhood evidence")
    fusion_score: Optional[float] = Field(default=None, description="Gene fusion evidence")
    cooccurrence_score: Optional[float] = Field(default=None, description="Co-occurrence evidence")
    
    class Config:
        extra = "allow"


class SequenceRecord(BaseModel):
    """Protein sequence and annotation from UniProt."""
    
    accession: str = Field(description="UniProt accession")
    gene_name: Optional[str] = Field(default=None, description="Gene symbol")
    protein_name: str = Field(description="Protein name")
    organism: str = Field(description="Source organism")
    sequence: str = Field(description="Amino acid sequence")
    length: int = Field(description="Sequence length")
    mass: Optional[float] = Field(default=None, description="Molecular mass (Da)")
    
    # Functional annotations
    function: Optional[str] = Field(default=None, description="Protein function")
    pathway: Optional[List[str]] = Field(default=None, description="Associated pathways")
    go_terms: Optional[List[str]] = Field(default=None, description="Gene Ontology terms")
    
    # Active sites and features
    active_sites: Optional[List[Dict[str, Any]]] = Field(default=None, description="Active site annotations")
    binding_sites: Optional[List[Dict[str, Any]]] = Field(default=None, description="Binding site annotations")
    domains: Optional[List[Dict[str, Any]]] = Field(default=None, description="Domain annotations")
    
    class Config:
        extra = "allow"


class StructureRecord(BaseModel):
    """Protein structure information from PDB/AlphaFold."""
    
    structure_id: str = Field(description="Structure identifier (PDB ID or AF-XXX)")
    source: str = Field(description="Structure source (PDB/AlphaFold)")
    method: Optional[str] = Field(default=None, description="Experimental method")
    resolution: Optional[float] = Field(default=None, description="Resolution (Angstroms)")
    
    # File information
    file_path: Optional[str] = Field(default=None, description="Local file path")
    file_format: Optional[str] = Field(default=None, description="File format (pdb/cif)")
    
    # Quality metrics
    confidence: Optional[float] = Field(default=None, description="Confidence score (AlphaFold)")
    r_factor: Optional[float] = Field(default=None, description="R-factor (experimental)")
    r_free: Optional[float] = Field(default=None, description="R-free (experimental)")
    
    # Chain information
    chains: Optional[List[str]] = Field(default=None, description="Available chains")
    ligands: Optional[List[str]] = Field(default=None, description="Bound ligands")
    
    class Config:
        extra = "allow"


class Target(BaseModel):
    """A biological target with all associated data."""
    
    # Core identifiers
    target_id: str = Field(description="Primary target identifier")
    gene_name: Optional[str] = Field(default=None, description="Gene symbol")
    uniprot_id: Optional[str] = Field(default=None, description="UniProt accession")
    
    # Discovery metadata
    discovery_source: List[str] = Field(description="Sources that identified this target")
    discovery_pathways: List[PathwayHit] = Field(default_factory=list, description="Associated pathways")
    
    # Characterization data
    expression_data: List[ExpressionRecord] = Field(default_factory=list, description="Expression profiles")
    ppi_data: List[PPIRecord] = Field(default_factory=list, description="Protein interactions")
    sequence_data: Optional[SequenceRecord] = Field(default=None, description="Sequence information")
    structure_data: List[StructureRecord] = Field(default_factory=list, description="Structure information")
    
    # Quality flags
    has_structure: bool = Field(default=False, description="Has 3D structure available")
    has_sequence: bool = Field(default=False, description="Has sequence information")
    tissue_specific: bool = Field(default=False, description="Shows tissue-specific expression")
    
    # Metadata
    created_at: datetime = Field(default_factory=datetime.now)
    updated_at: datetime = Field(default_factory=datetime.now)
    
    class Config:
        extra = "allow"
        
    def update_flags(self) -> None:
        """Update quality flags based on available data."""
        self.has_structure = len(self.structure_data) > 0
        self.has_sequence = self.sequence_data is not None
        self.tissue_specific = len(self.expression_data) > 0
        self.updated_at = datetime.now()
