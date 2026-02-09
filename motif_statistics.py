"""
Statistics Module
Calculate statistical significance of motif enrichment
"""

import numpy as np
from scipy import stats
from typing import List, Dict
import math


class MotifStatistics:
    """Statistical analysis for motif enrichment"""
    
    def __init__(self, base_probabilities: Dict[str, float] = None):
        """
        Initialize statistics calculator
        
        Args:
            base_probabilities: Dictionary of base probabilities (default: equal)
        """
        if base_probabilities is None:
            self.base_probs = {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25, 'U': 0.25}
        else:
            self.base_probs = base_probabilities
    
    def calculate_motif_probability(self, motif: str) -> float:
        """
        Calculate probability of a motif occurring by chance
        
        Args:
            motif: Sequence motif
            
        Returns:
            Probability of this exact motif at any position
        """
        prob = 1.0
        for base in motif:
            prob *= self.base_probs.get(base, 0.25)
        return prob
    
    def calculate_binomial_pvalue(self, motif: str, num_sequences_with_motif: int,
                                   total_sequences: int, 
                                   random_region_length: int) -> float:
        """
        Calculate p-value using binomial test
        
        Tests if motif appears in more sequences than expected by chance
        
        Args:
            motif: The sequence motif
            num_sequences_with_motif: Number of sequences containing the motif
            total_sequences: Total number of sequences
            random_region_length: Average length of random region
            
        Returns:
            P-value (probability of seeing this or more extreme by chance)
        """
        motif_len = len(motif)
        
        # Probability of this specific motif at one position
        p_motif_per_position = self.calculate_motif_probability(motif)
        
        # Number of possible positions in random region
        n_positions = max(1, random_region_length - motif_len + 1)
        
        # Probability that motif appears at least once in one sequence
        # P(at least one occurrence) = 1 - P(zero occurrences)
        p_in_sequence = 1 - (1 - p_motif_per_position) ** n_positions
        
        # Binomial test: probability of seeing motif in k or more sequences
        # when the per-sequence probability is p_in_sequence
        return stats.binomtest(
            num_sequences_with_motif,
            total_sequences,
            p_in_sequence,
            alternative='greater'
        ).pvalue
    
    def add_pvalues_to_motifs(self, motifs: List[Dict], 
                               mean_region_length: float) -> List[Dict]:
        """
        Add p-values to motif dictionaries
        
        Args:
            motifs: List of motif dictionaries
            mean_region_length: Mean length of random regions
            
        Returns:
            Updated motif list with p-values
        """
        for motif_dict in motifs:
            pvalue = self.calculate_binomial_pvalue(
                motif_dict['motif'],
                motif_dict['num_sequences'],
                motif_dict['total_sequences'],
                int(mean_region_length)
            )
            motif_dict['p_value'] = pvalue
            
        return motifs
    
    def apply_fdr_correction(self, motifs: List[Dict], 
                             method: str = 'fdr_bh') -> List[Dict]:
        """
        Apply False Discovery Rate correction to p-values
        
        Args:
            motifs: List of motif dictionaries with p_values
            method: Correction method ('bonferroni' or 'fdr_bh')
            
        Returns:
            Updated motif list with corrected p-values
        """
        if not motifs:
            return motifs
        
        # Extract p-values
        pvalues = [m.get('p_value', 1.0) for m in motifs]
        
        if method == 'bonferroni':
            corrected = self._bonferroni_correction(pvalues)
        elif method == 'fdr_bh':
            corrected = self._benjamini_hochberg_correction(pvalues)
        else:
            raise ValueError(f"Unknown correction method: {method}")
        
        # Add corrected p-values
        for motif_dict, adj_pval in zip(motifs, corrected):
            motif_dict['adjusted_p_value'] = adj_pval
            motif_dict['significant'] = adj_pval < 0.05
            
        return motifs
    
    def _bonferroni_correction(self, pvalues: List[float]) -> List[float]:
        """Bonferroni correction"""
        n = len(pvalues)
        return [min(p * n, 1.0) for p in pvalues]
    
    def _benjamini_hochberg_correction(self, pvalues: List[float]) -> List[float]:
        """Benjamini-Hochberg FDR correction"""
        n = len(pvalues)
        
        # Create list of (original_index, p-value) and sort by p-value
        indexed_pvalues = [(i, p) for i, p in enumerate(pvalues)]
        indexed_pvalues.sort(key=lambda x: x[1])
        
        # Calculate adjusted p-values
        adjusted = [0] * n
        prev_bh = 0
        
        for rank, (original_idx, pval) in enumerate(indexed_pvalues, 1):
            bh_value = pval * n / rank
            # Ensure monotonicity
            bh_value = max(bh_value, prev_bh)
            adjusted[original_idx] = min(bh_value, 1.0)
            prev_bh = bh_value
            
        return adjusted
    
    def calculate_enrichment_score(self, motif: str, 
                                    num_sequences_with_motif: int,
                                    total_sequences: int,
                                    random_region_length: int) -> float:
        """
        Calculate fold-enrichment over expected
        
        Args:
            motif: The sequence motif
            num_sequences_with_motif: Observed number of sequences
            total_sequences: Total number of sequences
            random_region_length: Average length of random region
            
        Returns:
            Fold enrichment (observed/expected)
        """
        motif_len = len(motif)
        p_motif_per_position = self.calculate_motif_probability(motif)
        n_positions = max(1, random_region_length - motif_len + 1)
        p_in_sequence = 1 - (1 - p_motif_per_position) ** n_positions
        
        expected_sequences = p_in_sequence * total_sequences
        observed_sequences = num_sequences_with_motif
        
        if expected_sequences > 0:
            return observed_sequences / expected_sequences
        else:
            return float('inf')
    
    def add_enrichment_scores(self, motifs: List[Dict], 
                              mean_region_length: float) -> List[Dict]:
        """
        Add enrichment scores to motif dictionaries
        
        Args:
            motifs: List of motif dictionaries
            mean_region_length: Mean length of random regions
            
        Returns:
            Updated motif list with enrichment scores
        """
        for motif_dict in motifs:
            enrichment = self.calculate_enrichment_score(
                motif_dict['motif'],
                motif_dict['num_sequences'],
                motif_dict['total_sequences'],
                int(mean_region_length)
            )
            motif_dict['fold_enrichment'] = enrichment
            
        return motifs
    
    def calculate_base_composition(self, sequences: List[str]) -> Dict[str, float]:
        """
        Calculate actual base composition from sequences
        
        Args:
            sequences: List of sequence strings
            
        Returns:
            Dictionary of base frequencies
        """
        all_bases = ''.join(sequences).upper()
        total = len(all_bases)
        
        if total == 0:
            return {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}
        
        composition = {
            'A': all_bases.count('A') / total,
            'C': all_bases.count('C') / total,
            'G': all_bases.count('G') / total,
            'T': all_bases.count('T') / total,
            'U': all_bases.count('U') / total
        }
        
        return composition
    
    def get_significance_summary(self, motifs: List[Dict], 
                                  alpha: float = 0.05) -> Dict:
        """
        Get summary statistics of significant motifs
        
        Args:
            motifs: List of motif dictionaries with p-values
            alpha: Significance threshold
            
        Returns:
            Summary dictionary
        """
        if not motifs:
            return {
                'total_motifs': 0,
                'significant_raw': 0,
                'significant_adjusted': 0,
                'min_p_value': None,
                'min_adjusted_p_value': None
            }
        
        significant_raw = sum(1 for m in motifs if m.get('p_value', 1) < alpha)
        significant_adj = sum(1 for m in motifs if m.get('adjusted_p_value', 1) < alpha)
        
        p_values = [m['p_value'] for m in motifs if 'p_value' in m]
        adj_p_values = [m['adjusted_p_value'] for m in motifs if 'adjusted_p_value' in m]
        
        return {
            'total_motifs': len(motifs),
            'significant_raw': significant_raw,
            'significant_adjusted': significant_adj,
            'min_p_value': min(p_values) if p_values else None,
            'min_adjusted_p_value': min(adj_p_values) if adj_p_values else None,
            'mean_p_value': sum(p_values) / len(p_values) if p_values else None
        }
