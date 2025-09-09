"""
MCP client modules for various biological databases.
"""

from .base import BaseMCPClient, MCPError, MCPConnectionError, MCPTimeoutError
from .kegg import KEGGClient
from .reactome import ReactomeClient
from .proteinatlas import ProteinAtlasClient
from .string import STRINGClient
from .uniprot import UniProtClient
from .pdb import PDBClient
from .chembl import ChEMBLClient

__all__ = [
    # Base classes
    "BaseMCPClient",
    "MCPError", 
    "MCPConnectionError",
    "MCPTimeoutError",
    
    # Specific clients
    "KEGGClient",
    "ReactomeClient", 
    "ProteinAtlasClient",
    "STRINGClient",
    "UniProtClient",
    "PDBClient",
    "ChEMBLClient",
]
