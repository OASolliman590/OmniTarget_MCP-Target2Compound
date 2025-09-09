"""Score normalization utilities."""

import numpy as np
from typing import List, Dict, Any
from loguru import logger


class ScoreNormalizer:
    """Normalize scores using different methods."""
    
    @staticmethod
    def z_score_normalize(scores: List[float]) -> List[float]:
        """Z-score normalization."""
        if not scores:
            return []
        
        scores_array = np.array(scores)
        mean = np.mean(scores_array)
        std = np.std(scores_array)
        
        if std == 0:
            return [0.0] * len(scores)
        
        normalized = (scores_array - mean) / std
        return normalized.tolist()
    
    @staticmethod
    def min_max_normalize(scores: List[float]) -> List[float]:
        """Min-max normalization to [0, 1]."""
        if not scores:
            return []
        
        scores_array = np.array(scores)
        min_val = np.min(scores_array)
        max_val = np.max(scores_array)
        
        if max_val == min_val:
            return [0.5] * len(scores)
        
        normalized = (scores_array - min_val) / (max_val - min_val)
        return normalized.tolist()
    
    @staticmethod
    def rank_normalize(scores: List[float]) -> List[float]:
        """Rank-based normalization."""
        if not scores:
            return []
        
        # Convert to ranks
        ranks = np.argsort(np.argsort(scores))
        normalized = ranks / (len(scores) - 1) if len(scores) > 1 else [0.5]
        return normalized.tolist()
