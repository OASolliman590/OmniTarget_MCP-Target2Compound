"""Score fusion utilities."""

from typing import Dict, List, Any
import numpy as np
from loguru import logger

from .normalize import ScoreNormalizer


class ScoreFusion:
    """Combine multiple scores into integrated rankings."""
    
    def __init__(self, weights: Dict[str, float]):
        """Initialize with scoring weights."""
        self.weights = weights
        self.normalizer = ScoreNormalizer()
    
    def combine_scores(
        self,
        score_data: List[Dict[str, Any]],
        normalization_method: str = "zscore"
    ) -> List[Dict[str, Any]]:
        """Combine multiple score types into integrated scores.
        
        Args:
            score_data: List of records with different score types
            normalization_method: Method for score normalization
            
        Returns:
            List of records with combined scores
        """
        try:
            # Extract individual score types
            docking_scores = []
            evidence_scores = []
            
            for record in score_data:
                docking_scores.append(record.get("docking_score", 0.0))
                evidence_scores.append(record.get("evidence_score", 0.0))
            
            # Normalize scores
            if normalization_method == "zscore":
                docking_norm = self.normalizer.z_score_normalize(docking_scores)
                evidence_norm = self.normalizer.z_score_normalize(evidence_scores)
            elif normalization_method == "minmax":
                docking_norm = self.normalizer.min_max_normalize(docking_scores)
                evidence_norm = self.normalizer.min_max_normalize(evidence_scores)
            else:
                docking_norm = docking_scores
                evidence_norm = evidence_scores
            
            # Combine scores
            combined_scores = []
            for i, record in enumerate(score_data):
                combined_score = (
                    self.weights.get("docking", 0) * docking_norm[i] +
                    self.weights.get("evidence", 0) * evidence_norm[i]
                )
                
                record_with_combined = record.copy()
                record_with_combined.update({
                    "combined_score": combined_score,
                    "docking_normalized": docking_norm[i],
                    "evidence_normalized": evidence_norm[i]
                })
                combined_scores.append(record_with_combined)
            
            # Sort by combined score (descending)
            combined_scores.sort(key=lambda x: x["combined_score"], reverse=True)
            
            # Add ranks
            for i, record in enumerate(combined_scores):
                record["rank"] = i + 1
            
            return combined_scores
            
        except Exception as e:
            logger.error(f"Error combining scores: {e}")
            return score_data
