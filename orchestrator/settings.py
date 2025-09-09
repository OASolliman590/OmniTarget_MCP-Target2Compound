"""
Application settings and configuration management.
"""

from typing import List, Optional
from pydantic import Field
from pydantic_settings import BaseSettings


class MCPServerSettings(BaseSettings):
    """MCP Server connection settings."""
    
    kegg_url: str = Field(default="http://localhost:3000", env="MCP_KEGG_URL")
    reactome_url: str = Field(default="http://localhost:3001", env="MCP_REACTOME_URL")
    proteinatlas_url: str = Field(default="http://localhost:3002", env="MCP_PROTEINATLAS_URL")
    string_url: str = Field(default="http://localhost:3003", env="MCP_STRING_URL")
    uniprot_url: str = Field(default="http://localhost:3004", env="MCP_UNIPROT_URL")
    pdb_url: str = Field(default="http://localhost:3005", env="MCP_PDB_URL")
    chembl_url: str = Field(default="http://localhost:3006", env="MCP_CHEMBL_URL")
    
    # Request settings
    request_timeout: int = Field(default=30, env="MCP_REQUEST_TIMEOUT")
    max_retries: int = Field(default=3, env="MCP_MAX_RETRIES")
    rate_limit: int = Field(default=10, env="MCP_RATE_LIMIT")  # requests per second


class ModelSettings(BaseSettings):
    """ML Model settings."""
    
    # DeepDTA settings
    deepdta_dir: str = Field(default="third_party/DeepDTA", env="DEEPDTA_DIR")
    deepdta_model_path: Optional[str] = Field(default=None, env="DEEPDTA_MODEL_PATH")
    deepdta_weights_url: str = Field(
        default="https://github.com/hkmztrk/DeepDTA/releases/download/v1.0/model_weights.h5",
        env="DEEPDTA_WEIGHTS_URL"
    )
    
    # GeminiMol settings
    geminimol_model_path: Optional[str] = Field(default=None, env="GEMINIMOL_MODEL_PATH")
    
    # Vina settings
    vina_binary_path: str = Field(default="/Users/omara.soliman/opt/miniconda3/envs/mcp-drug-discovery/bin/vina_custom", env="VINA_BINARY_PATH")
    vina_exhaustiveness: int = Field(default=8, env="VINA_EXHAUSTIVENESS")
    vina_num_modes: int = Field(default=9, env="VINA_NUM_MODES")


class CacheSettings(BaseSettings):
    """Caching configuration."""
    
    redis_url: str = Field(default="redis://localhost:6379", env="REDIS_URL")
    cache_ttl: int = Field(default=3600, env="CACHE_TTL")  # 1 hour
    cache_dir: str = Field(default="./data/cache", env="CACHE_DIR")
    
    # Cache keys
    mcp_cache_prefix: str = "mcp:"
    pipeline_cache_prefix: str = "pipeline:"


class FileSettings(BaseSettings):
    """File and directory settings."""
    
    temp_dir: str = Field(default="./data/temp", env="TEMP_DIR")
    output_dir: str = Field(default="./data/outputs", env="OUTPUT_DIR")
    compounds_dir: str = Field(default="./data/compounds", env="COMPOUNDS_DIR")
    
    # Third party paths
    deepdta_dir: str = Field(default="./third_party/DeepDTA", env="DEEPDTA_DIR")
    geminimol_dir: str = Field(default="./third_party/GeminiMol", env="GEMINIMOL_DIR")
    vina_dir: str = Field(default="./third_party/autodock-vina", env="VINA_DIR")


class LoggingSettings(BaseSettings):
    """Logging configuration."""
    
    log_level: str = Field(default="INFO", env="LOG_LEVEL")
    log_format: str = Field(
        default="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level: <8}</level> | <cyan>{name}</cyan>:<cyan>{function}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>",
        env="LOG_FORMAT"
    )
    log_file: Optional[str] = Field(default=None, env="LOG_FILE")


class APISettings(BaseSettings):
    """FastAPI application settings."""
    
    title: str = "MCP Drug Discovery Orchestrator"
    description: str = "API for orchestrating drug discovery pipelines using MCP servers"
    version: str = "0.1.0"
    
    host: str = Field(default="0.0.0.0", env="API_HOST")
    port: int = Field(default=8000, env="API_PORT")
    debug: bool = Field(default=False, env="API_DEBUG")
    
    # Security
    api_key: Optional[str] = Field(default=None, env="API_KEY")
    allowed_hosts: List[str] = Field(default=["*"], env="API_ALLOWED_HOSTS")
    
    # Rate limiting
    rate_limit: int = Field(default=100, env="API_RATE_LIMIT")  # requests per minute


class Settings(BaseSettings):
    """Main application settings."""
    
    # Sub-settings
    mcp: MCPServerSettings = MCPServerSettings()
    models: ModelSettings = ModelSettings()
    cache: CacheSettings = CacheSettings()
    files: FileSettings = FileSettings()
    logging: LoggingSettings = LoggingSettings()
    api: APISettings = APISettings()
    
    # Global settings
    environment: str = Field(default="development", env="ENVIRONMENT")
    debug: bool = Field(default=False, env="DEBUG")
    
    class Config:
        env_file = ".env"
        env_file_encoding = "utf-8"
        case_sensitive = False


# Global settings instance
settings = Settings()
