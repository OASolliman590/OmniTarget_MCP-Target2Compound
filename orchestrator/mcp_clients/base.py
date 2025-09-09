"""
Base MCP client with common functionality.
"""

import asyncio
import hashlib
import json
from typing import Any, Dict, List, Optional, Union
from abc import ABC, abstractmethod

import httpx
from loguru import logger

from ..settings import settings


class MCPError(Exception):
    """Base exception for MCP client errors."""
    pass


class MCPConnectionError(MCPError):
    """Connection-related MCP errors."""
    pass


class MCPTimeoutError(MCPError):
    """Timeout-related MCP errors."""
    pass


class BaseMCPClient(ABC):
    """Base class for all MCP clients."""
    
    def __init__(
        self,
        base_url: str,
        timeout: float = None,
        max_retries: int = None,
        rate_limit: float = None
    ):
        """Initialize the MCP client.
        
        Args:
            base_url: Base URL for the MCP server
            timeout: Request timeout in seconds
            max_retries: Maximum number of retry attempts
            rate_limit: Maximum requests per second
        """
        self.base_url = base_url.rstrip('/')
        self.timeout = timeout or settings.mcp.request_timeout
        self.max_retries = max_retries or settings.mcp.max_retries
        self.rate_limit = rate_limit or settings.mcp.rate_limit
        
        # Rate limiting
        self._last_request_time = 0.0
        self._request_interval = 1.0 / self.rate_limit if self.rate_limit > 0 else 0.0
        
        # HTTP client configuration
        self._client = httpx.AsyncClient(
            timeout=httpx.Timeout(self.timeout),
            limits=httpx.Limits(max_connections=10, max_keepalive_connections=5),
            headers={
                "User-Agent": "MCP-Drug-Discovery-Orchestrator/0.1.0",
                "Accept": "application/json",
                "Content-Type": "application/json"
            }
        )
        
        logger.info(f"Initialized {self.__class__.__name__} for {self.base_url}")
    
    async def __aenter__(self):
        """Async context manager entry."""
        return self
    
    async def __aexit__(self, exc_type, exc_val, exc_tb):
        """Async context manager exit."""
        await self.close()
    
    async def close(self) -> None:
        """Close the HTTP client."""
        await self._client.aclose()
    
    def _cache_key(self, endpoint: str, params: Dict[str, Any]) -> str:
        """Generate cache key for request."""
        key_data = {
            "endpoint": endpoint,
            "params": params,
            "client": self.__class__.__name__
        }
        key_str = json.dumps(key_data, sort_keys=True)
        return hashlib.md5(key_str.encode()).hexdigest()
    
    async def _rate_limit(self) -> None:
        """Apply rate limiting."""
        if self._request_interval > 0:
            current_time = asyncio.get_event_loop().time()
            time_since_last = current_time - self._last_request_time
            
            if time_since_last < self._request_interval:
                await asyncio.sleep(self._request_interval - time_since_last)
            
            self._last_request_time = asyncio.get_event_loop().time()
    
    async def _make_request(
        self,
        method: str,
        endpoint: str,
        params: Optional[Dict[str, Any]] = None,
        data: Optional[Dict[str, Any]] = None,
        headers: Optional[Dict[str, str]] = None
    ) -> Dict[str, Any]:
        """Make HTTP request with retry logic.
        
        Args:
            method: HTTP method (GET, POST, etc.)
            endpoint: API endpoint
            params: Query parameters
            data: Request body data
            headers: Additional headers
            
        Returns:
            Response data as dictionary
            
        Raises:
            MCPConnectionError: On connection issues
            MCPTimeoutError: On timeout
            MCPError: On other errors
        """
        url = f"{self.base_url}/{endpoint.lstrip('/')}"
        params = params or {}
        data = data or {}
        headers = headers or {}
        
        # Apply rate limiting
        await self._rate_limit()
        
        last_exception = None
        
        for attempt in range(self.max_retries + 1):
            try:
                logger.debug(f"Making {method} request to {url} (attempt {attempt + 1})")
                
                response = await self._client.request(
                    method=method,
                    url=url,
                    params=params,
                    json=data if method.upper() in ['POST', 'PUT', 'PATCH'] else None,
                    headers=headers
                )
                
                # Check for HTTP errors
                if response.status_code >= 400:
                    error_msg = f"HTTP {response.status_code}: {response.text}"
                    logger.warning(f"Request failed: {error_msg}")
                    
                    if response.status_code >= 500:
                        # Server error - retry
                        last_exception = MCPError(error_msg)
                        continue
                    else:
                        # Client error - don't retry
                        raise MCPError(error_msg)
                
                # Parse response
                try:
                    result = response.json()
                except json.JSONDecodeError:
                    result = {"text": response.text}
                
                logger.debug(f"Request successful: {len(str(result))} chars")
                return result
                
            except httpx.TimeoutException as e:
                last_exception = MCPTimeoutError(f"Request timeout: {e}")
                logger.warning(f"Request timeout (attempt {attempt + 1}): {e}")
                
            except httpx.ConnectError as e:
                last_exception = MCPConnectionError(f"Connection error: {e}")
                logger.warning(f"Connection error (attempt {attempt + 1}): {e}")
                
            except httpx.RequestError as e:
                last_exception = MCPError(f"Request error: {e}")
                logger.warning(f"Request error (attempt {attempt + 1}): {e}")
            
            # Wait before retry (exponential backoff)
            if attempt < self.max_retries:
                wait_time = 2 ** attempt
                logger.debug(f"Waiting {wait_time}s before retry")
                await asyncio.sleep(wait_time)
        
        # All retries failed
        raise last_exception or MCPError("Request failed after all retries")
    
    async def get(
        self,
        endpoint: str,
        params: Optional[Dict[str, Any]] = None,
        **kwargs
    ) -> Dict[str, Any]:
        """Make GET request."""
        return await self._make_request("GET", endpoint, params=params, **kwargs)
    
    async def post(
        self,
        endpoint: str,
        data: Optional[Dict[str, Any]] = None,
        **kwargs
    ) -> Dict[str, Any]:
        """Make POST request."""
        return await self._make_request("POST", endpoint, data=data, **kwargs)
    
    async def health_check(self) -> bool:
        """Check if the MCP server is healthy.
        
        Returns:
            True if server is responsive, False otherwise
        """
        try:
            # Try to make a simple request to check connectivity
            response = await self._client.get(
                f"{self.base_url}/health",
                timeout=5.0
            )
            return response.status_code == 200
        except Exception as e:
            logger.warning(f"Health check failed for {self.base_url}: {e}")
            return False
    
    @abstractmethod
    async def test_connection(self) -> bool:
        """Test connection with a simple API call.
        
        Each client should implement this with a lightweight test request.
        
        Returns:
            True if connection works, False otherwise
        """
        pass
    
    # Utility methods for common operations
    
    def _normalize_gene_id(self, gene_id: str) -> str:
        """Normalize gene identifier."""
        return gene_id.strip().upper()
    
    def _normalize_protein_id(self, protein_id: str) -> str:
        """Normalize protein identifier."""
        return protein_id.strip().upper()
    
    def _normalize_species(self, species: str) -> str:
        """Normalize species name."""
        return species.strip().lower()
    
    def _validate_response(self, response: Dict[str, Any], required_fields: List[str]) -> bool:
        """Validate that response contains required fields.
        
        Args:
            response: Response dictionary
            required_fields: List of required field names
            
        Returns:
            True if all required fields present
        """
        for field in required_fields:
            if field not in response:
                logger.warning(f"Missing required field '{field}' in response")
                return False
        return True
