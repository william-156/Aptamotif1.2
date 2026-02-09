"""
Motif Finder Module
Discovers enriched k-mer motifs in aptamer sequences
"""

from typing import List, Dict, Set, Tuple
from collections import defaultdict, Counter
import itertools


class MotifFinder:
    """Find enriched motifs in aptamer sequences"""
    
    def __init__(self, min_length: int = 5, max_length: int = 15, 
                 min_sequences: int = 2):
        """
        Initialize motif finder
        
        Args:
            min_length: Minimum motif length (bp)
            max_length: Maximum motif length (bp)
            min_sequences: Minimum number of sequences motif must appear in
        """
        self.min_length = min_length
        self.max_length = max_length
        self.min_sequences = min_sequences
        
    def find_all_kmers(self, sequences: List[str], k: int) -> Dict[str, Set[int]]:
        """
        Find all k-mers and which sequences they appear in
        
        Args:
            sequences: List of sequence strings
            k: k-mer length
            
        Returns:
            Dictionary mapping k-mer -> set of sequence indices
        """
        kmer_to_seqs = defaultdict(set)
        
        for seq_idx, seq in enumerate(sequences):
            seq = seq.upper()
            # Find all k-mers in this sequence
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i+k]
                # Only consider k-mers with standard bases
                if all(base in 'ACGTU' for base in kmer):
                    kmer_to_seqs[kmer].add(seq_idx)
                    
        return kmer_to_seqs
    
    def find_enriched_motifs(self, sequences: List[Dict]) -> List[Dict]:
        """
        Find all enriched motifs across different k-mer lengths
        
        Args:
            sequences: List of sequence dictionaries with 'random_region'
            
        Returns:
            List of motif dictionaries with statistics
        """
        # Extract just the random regions
        seq_strings = [s['random_region'] for s in sequences if 'random_region' in s]
        
        if not seq_strings:
            return []
        
        all_motifs = []
        
        # Search for motifs of each length
        for k in range(self.min_length, self.max_length + 1):
            kmer_dict = self.find_all_kmers(seq_strings, k)
            
            # Filter by minimum occurrence
            for kmer, seq_indices in kmer_dict.items():
                if len(seq_indices) >= self.min_sequences:
                    # Calculate statistics
                    motif_info = self._calculate_motif_stats(
                        kmer, seq_indices, seq_strings
                    )
                    all_motifs.append(motif_info)
        
        # Sort by number of sequences (descending) then by p-value
        all_motifs.sort(key=lambda x: (-x['num_sequences'], x.get('p_value', 1)))
        
        return all_motifs
    
    def _calculate_motif_stats(self, kmer: str, seq_indices: Set[int], 
                                sequences: List[str]) -> Dict:
        """
        Calculate statistics for a motif
        
        Args:
            kmer: The k-mer sequence
            seq_indices: Set of sequence indices containing this k-mer
            sequences: All sequences
            
        Returns:
            Dictionary with motif statistics
        """
        num_seqs_with_motif = len(seq_indices)
        total_seqs = len(sequences)
        
        # Count occurrences per sequence
        occurrences_per_seq = []
        total_occurrences = 0
        
        for seq_idx in seq_indices:
            seq = sequences[seq_idx]
            count = seq.count(kmer)
            occurrences_per_seq.append(count)
            total_occurrences += count
        
        # Calculate positions where motif occurs
        positions = []
        for seq_idx in seq_indices:
            seq = sequences[seq_idx]
            start = 0
            while True:
                pos = seq.find(kmer, start)
                if pos == -1:
                    break
                positions.append(pos)
                start = pos + 1
        
        return {
            'motif': kmer,
            'length': len(kmer),
            'num_sequences': num_seqs_with_motif,
            'total_sequences': total_seqs,
            'frequency': num_seqs_with_motif / total_seqs,
            'total_occurrences': total_occurrences,
            'mean_occurrences_per_seq': total_occurrences / num_seqs_with_motif,
            'sequence_indices': list(seq_indices),
            'positions': positions,
            'gc_content': self._calculate_gc_content(kmer)
        }
    
    def _calculate_gc_content(self, sequence: str) -> float:
        """Calculate GC content of sequence"""
        gc_count = sequence.count('G') + sequence.count('C')
        return gc_count / len(sequence) if len(sequence) > 0 else 0
    
    def find_reverse_complements(self, motifs: List[Dict]) -> List[Dict]:
        """
        Identify motifs that are reverse complements of each other
        
        Args:
            motifs: List of motif dictionaries
            
        Returns:
            Updated motif list with reverse complement info
        """
        complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'U': 'A'}
        
        # Create mapping of motifs to their indices
        motif_to_idx = {m['motif']: i for i, m in enumerate(motifs)}
        
        for motif_dict in motifs:
            seq = motif_dict['motif']
            # Calculate reverse complement
            rev_comp = ''.join(complement.get(base, base) for base in reversed(seq))
            
            if rev_comp in motif_to_idx and rev_comp != seq:
                motif_dict['reverse_complement'] = rev_comp
                motif_dict['has_reverse_complement'] = True
            else:
                motif_dict['has_reverse_complement'] = False
                
        return motifs
    
    def filter_redundant_motifs(self, motifs: List[Dict]) -> List[Dict]:
        """
        Remove redundant motifs that are substrings of longer motifs
        with the same sequence distribution
        
        Strategy: Keep longer motifs that appear in the same sequences as shorter ones
        
        Args:
            motifs: List of motif dictionaries
            
        Returns:
            Filtered list of non-redundant motifs
        """
        if not motifs:
            return []
        
        # Sort by length (longest first) then by num_sequences
        sorted_motifs = sorted(motifs, 
                             key=lambda x: (-len(x['motif']), -x['num_sequences']))
        
        kept_motifs = []
        
        for candidate in sorted_motifs:
            candidate_seq = candidate['motif']
            candidate_indices = set(candidate['sequence_indices'])
            
            # Check if this motif is redundant with any kept motif
            is_redundant = False
            
            for kept in kept_motifs:
                kept_seq = kept['motif']
                kept_indices = set(kept['sequence_indices'])
                
                # Check if candidate is substring of kept motif
                if candidate_seq in kept_seq:
                    # Check if they appear in >80% the same sequences
                    overlap = len(candidate_indices & kept_indices)
                    candidate_coverage = overlap / len(candidate_indices) if candidate_indices else 0
                    
                    if candidate_coverage > 0.8:
                        is_redundant = True
                        break
                
                # Check if kept is substring of candidate (shouldn't happen due to sorting)
                elif kept_seq in candidate_seq:
                    # Check if they appear in >80% the same sequences
                    overlap = len(candidate_indices & kept_indices)
                    kept_coverage = overlap / len(kept_indices) if kept_indices else 0
                    
                    if kept_coverage > 0.8:
                        # Replace kept with candidate (longer version)
                        kept_motifs.remove(kept)
                        break
            
            if not is_redundant:
                kept_motifs.append(candidate)
        
        # Re-sort by significance
        kept_motifs.sort(key=lambda x: (-x['num_sequences'], 
                                        x.get('adjusted_p_value', 1)))
        
        return kept_motifs
    
    def find_merged_motifs(self, motifs: List[Dict], sequences: List[Dict]) -> List[Dict]:
        """
        Attempt to merge overlapping short motifs into longer consensus motifs
        
        For example, if TATGG and CCTAT both appear in the same sequences,
        check if CCTATGG exists in those sequences
        
        Args:
            motifs: List of motif dictionaries
            sequences: List of sequence dictionaries
            
        Returns:
            List of potential merged motifs
        """
        merged_candidates = []
        seq_strings = [s['random_region'] for s in sequences if 'random_region' in s]
        
        # Look for motif pairs that share most of their sequences
        for i, motif1 in enumerate(motifs):
            for motif2 in motifs[i+1:]:
                # Check overlap in sequences
                indices1 = set(motif1['sequence_indices'])
                indices2 = set(motif2['sequence_indices'])
                overlap = indices1 & indices2
                
                # If >70% overlap, try to find merged version
                if len(overlap) / max(len(indices1), len(indices2)) > 0.7:
                    seq1 = motif1['motif']
                    seq2 = motif2['motif']
                    
                    # Try different merge patterns
                    merge_patterns = []
                    
                    # Check if seq2 overlaps with end of seq1
                    for shift in range(1, min(len(seq1), len(seq2))):
                        if seq1[-shift:] == seq2[:shift]:
                            merged = seq1 + seq2[shift:]
                            merge_patterns.append(merged)
                    
                    # Check if seq1 overlaps with end of seq2
                    for shift in range(1, min(len(seq1), len(seq2))):
                        if seq2[-shift:] == seq1[:shift]:
                            merged = seq2 + seq1[shift:]
                            merge_patterns.append(merged)
                    
                    # Check if one is substring of potential larger motif
                    # (e.g., CCTAT and TATGG might both be part of CCTATGG)
                    for seq_idx in overlap:
                        if seq_idx < len(seq_strings):
                            full_seq = seq_strings[seq_idx]
                            # Find if both motifs appear near each other
                            pos1 = full_seq.find(seq1)
                            pos2 = full_seq.find(seq2)
                            
                            if pos1 != -1 and pos2 != -1 and abs(pos1 - pos2) < 10:
                                start = min(pos1, pos2)
                                end = max(pos1 + len(seq1), pos2 + len(seq2))
                                potential_merged = full_seq[start:end]
                                if len(potential_merged) <= 20:  # Reasonable size
                                    merge_patterns.append(potential_merged)
                    
                    # Test each merge pattern
                    for merged_seq in set(merge_patterns):
                        # Count occurrences
                        merged_indices = []
                        for idx, seq in enumerate(seq_strings):
                            if merged_seq in seq:
                                merged_indices.append(idx)
                        
                        # If merged version appears in most of the overlapping sequences
                        if len(merged_indices) >= len(overlap) * 0.7:
                            merged_candidates.append({
                                'motif': merged_seq,
                                'length': len(merged_seq),
                                'num_sequences': len(merged_indices),
                                'sequence_indices': merged_indices,
                                'merged_from': [seq1, seq2],
                                'is_merged': True
                            })
        
        return merged_candidates
    
    def create_motif_matrix(self, motifs: List[Dict], sequences: List[Dict]) -> Tuple:
        """
        Create a binary matrix of motif presence/absence
        
        Args:
            motifs: List of motif dictionaries
            sequences: List of sequence dictionaries
            
        Returns:
            Tuple of (matrix, motif_names, sequence_names)
        """
        import numpy as np
        
        n_seqs = len(sequences)
        n_motifs = len(motifs)
        
        matrix = np.zeros((n_seqs, n_motifs), dtype=int)
        
        for j, motif_dict in enumerate(motifs):
            for seq_idx in motif_dict['sequence_indices']:
                if seq_idx < n_seqs:
                    matrix[seq_idx, j] = 1
        
        motif_names = [m['motif'] for m in motifs]
        sequence_names = [s['id'] for s in sequences]
        
        return matrix, motif_names, sequence_names
    
    def create_pairwise_similarity_matrix(self, motifs: List[Dict], 
                                         sequences: List[Dict]) -> Tuple:
        """
        Create a matrix showing shared significant motifs between each pair of sequences
        
        Args:
            motifs: List of significant motif dictionaries
            sequences: List of sequence dictionaries
            
        Returns:
            Tuple of (similarity_matrix, sequence_names)
        """
        import numpy as np
        
        n_seqs = len(sequences)
        similarity_matrix = np.zeros((n_seqs, n_seqs), dtype=int)
        
        # For each motif, increment count for each pair of sequences that share it
        for motif_dict in motifs:
            seq_indices = motif_dict['sequence_indices']
            
            # For each pair of sequences sharing this motif
            for i in seq_indices:
                for j in seq_indices:
                    if i < n_seqs and j < n_seqs:
                        similarity_matrix[i, j] += 1
        
        sequence_names = [s['id'] for s in sequences]
        
        return similarity_matrix, sequence_names
