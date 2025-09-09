"""
Compound and molecular data schemas.
"""

from typing import List, Optional, Dict, Any
from pydantic import BaseModel, Field, validator
from datetime import datetime


class CompoundFeatures(BaseModel):
    """Molecular features and representations for compounds."""
    
    # Basic molecular properties
    molecular_weight: Optional[float] = Field(default=None, description="Molecular weight (g/mol)")
    logp: Optional[float] = Field(default=None, description="Log partition coefficient")
    tpsa: Optional[float] = Field(default=None, description="Topological polar surface area")
    hbd: Optional[int] = Field(default=None, description="Hydrogen bond donors")
    hba: Optional[int] = Field(default=None, description="Hydrogen bond acceptors")
    rotatable_bonds: Optional[int] = Field(default=None, description="Number of rotatable bonds")
    
    # Drug-likeness metrics
    lipinski_violations: Optional[int] = Field(default=None, description="Lipinski rule violations")
    qed_score: Optional[float] = Field(default=None, description="Quantitative drug-likeness score")
    
    # Molecular representations
    geminimol_embedding: Optional[List[float]] = Field(default=None, description="GeminiMol embedding vector")
    ouroboros_features: Optional[List[float]] = Field(default=None, description="Ouroboros feature vector")
    
    # Fingerprints
    morgan_fp: Optional[List[int]] = Field(default=None, description="Morgan fingerprint")
    maccs_fp: Optional[List[int]] = Field(default=None, description="MACCS fingerprint")
    
    class Config:
        extra = "allow"


class Compound(BaseModel):
    """A chemical compound with identifiers and features."""
    
    # Core identifiers
    compound_id: str = Field(description="Internal compound identifier")
    smiles: str = Field(description="Canonical SMILES representation")
    
    # External identifiers
    inchi: Optional[str] = Field(default=None, description="InChI identifier")
    inchi_key: Optional[str] = Field(default=None, description="InChI key")
    chembl_id: Optional[str] = Field(default=None, description="ChEMBL identifier")
    pubchem_cid: Optional[str] = Field(default=None, description="PubChem CID")
    
    # Names and synonyms
    name: Optional[str] = Field(default=None, description="Compound name")
    synonyms: List[str] = Field(default_factory=list, description="Alternative names")
    
    # Source information
    source_file: Optional[str] = Field(default=None, description="Source file path")
    source_line: Optional[int] = Field(default=None, description="Line number in source file")
    
    # Molecular features
    features: Optional[CompoundFeatures] = Field(default=None, description="Computed molecular features")
    
    # 3D structure information
    sdf_data: Optional[str] = Field(default=None, description="SDF format data with 3D coordinates")
    conformer_count: Optional[int] = Field(default=None, description="Number of generated conformers")
    
    # Quality flags
    is_valid: bool = Field(default=True, description="Valid SMILES and structure")
    is_drug_like: Optional[bool] = Field(default=None, description="Passes drug-likeness filters")
    has_3d_coords: bool = Field(default=False, description="Has 3D coordinates")
    
    # Metadata
    created_at: datetime = Field(default_factory=datetime.now)
    updated_at: datetime = Field(default_factory=datetime.now)
    
    @validator('smiles')
    def validate_smiles(cls, v):
        """Basic SMILES validation."""
        if not v or not v.strip():
            raise ValueError("SMILES cannot be empty")
        
        # Basic character validation
        invalid_chars = set(v) - set("CONSPFBrClIHc()[]#=-+1234567890@/\\%.")
        if invalid_chars:
            # Allow some flexibility for unusual atoms/formats
            pass
            
        return v.strip()
    
    @validator('compound_id')
    def validate_compound_id(cls, v):
        """Ensure compound ID is valid."""
        if not v or not v.strip():
            raise ValueError("Compound ID cannot be empty")
        return v.strip()
    
    class Config:
        extra = "allow"
        
    def update_flags(self) -> None:
        """Update quality flags based on available data."""
        if self.features:
            # Check drug-likeness based on Lipinski's rule
            violations = 0
            if self.features.molecular_weight and self.features.molecular_weight > 500:
                violations += 1
            if self.features.logp and self.features.logp > 5:
                violations += 1
            if self.features.hbd and self.features.hbd > 5:
                violations += 1
            if self.features.hba and self.features.hba > 10:
                violations += 1
                
            self.features.lipinski_violations = violations
            self.is_drug_like = violations <= 1
        
        self.has_3d_coords = self.sdf_data is not None
        self.updated_at = datetime.now()
    
    def to_dict(self, include_features: bool = True) -> Dict[str, Any]:
        """Convert to dictionary representation."""
        data = self.dict()
        if not include_features:
            data.pop('features', None)
        return data
