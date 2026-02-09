"""
Sequence Parser Module
Handles input sequences, primer trimming, and format conversion
"""

from Bio import SeqIO
from Bio.Seq import Seq
from typing import List, Dict, Tuple
import re


class SequenceParser:
    """Parse and process aptamer sequences"""
    
    def __init__(self, forward_primer: str, reverse_primer_region: str, 
                 molecule_type: str = "RNA"):
        """
        Initialize parser with primer sequences
        
        Args:
            forward_primer: 5' to 3' forward primer sequence
            reverse_primer_region: The sequence region that appears after N-region
                                  (reverse complement of reverse primer)
            molecule_type: "RNA" or "ssDNA"
        """
        self.forward_primer = forward_primer.upper().replace(" ", "").replace("-", "")
        self.reverse_primer_region = reverse_primer_region.upper().replace(" ", "").replace("-", "")
        self.molecule_type = molecule_type
        
    def parse_sequences(self, input_text: str, file_format: str = "auto") -> List[Dict]:
        """
        Parse input sequences from text or file
        
        Args:
            input_text: Raw sequence data
            file_format: "fasta", "text", or "auto"
            
        Returns:
            List of dictionaries with sequence info
        """
        sequences = []
        
        # Detect format
        if file_format == "auto":
            if input_text.strip().startswith(">"):
                file_format = "fasta"
            else:
                file_format = "text"
        
        if file_format == "fasta":
            sequences = self._parse_fasta(input_text)
        else:
            sequences = self._parse_text(input_text)
            
        return sequences
    
    def _parse_fasta(self, fasta_text: str) -> List[Dict]:
        """Parse FASTA format"""
        sequences = []
        current_header = None
        current_seq = []
        
        for line in fasta_text.strip().split('\n'):
            line = line.strip()
            if line.startswith('>'):
                if current_header is not None:
                    seq_str = ''.join(current_seq).upper()
                    sequences.append({
                        'id': current_header,
                        'full_sequence': seq_str
                    })
                current_header = line[1:].strip()
                current_seq = []
            else:
                current_seq.append(line)
        
        # Don't forget last sequence
        if current_header is not None:
            seq_str = ''.join(current_seq).upper()
            sequences.append({
                'id': current_header,
                'full_sequence': seq_str
            })
            
        return sequences
    
    def _parse_text(self, text: str) -> List[Dict]:
        """Parse plain text (one sequence per line)"""
        sequences = []
        for i, line in enumerate(text.strip().split('\n'), 1):
            line = line.strip()
            if line and not line.startswith('#'):  # Skip comments
                sequences.append({
                    'id': f'Seq_{i}',
                    'full_sequence': line.upper()
                })
        return sequences
    
    def extract_random_region(self, sequences: List[Dict]) -> List[Dict]:
        """
        Extract the N-region (random region) from sequences
        
        Args:
            sequences: List of sequence dictionaries
            
        Returns:
            Updated list with extracted random regions
        """
        processed = []
        
        for seq_dict in sequences:
            full_seq = seq_dict['full_sequence']
            
            # Find forward primer
            fwd_idx = full_seq.find(self.forward_primer)
            if fwd_idx == -1:
                # Try fuzzy matching (allow 1-2 mismatches)
                fwd_idx = self._fuzzy_find(full_seq, self.forward_primer, max_mismatches=2)
            
            # Find reverse primer region
            rev_idx = full_seq.find(self.reverse_primer_region)
            if rev_idx == -1:
                rev_idx = self._fuzzy_find(full_seq, self.reverse_primer_region, max_mismatches=2)
            
            if fwd_idx != -1 and rev_idx != -1:
                # Extract random region
                start = fwd_idx + len(self.forward_primer)
                end = rev_idx
                random_region = full_seq[start:end]
                
                seq_dict['random_region'] = random_region
                seq_dict['random_region_length'] = len(random_region)
                seq_dict['primers_found'] = True
            else:
                seq_dict['random_region'] = full_seq  # Use full sequence if primers not found
                seq_dict['random_region_length'] = len(full_seq)
                seq_dict['primers_found'] = False
                
            processed.append(seq_dict)
        
        return processed
    
    def _fuzzy_find(self, sequence: str, pattern: str, max_mismatches: int = 2) -> int:
        """
        Find pattern in sequence allowing mismatches
        
        Returns:
            Index of match or -1 if not found
        """
        pattern_len = len(pattern)
        seq_len = len(sequence)
        
        for i in range(seq_len - pattern_len + 1):
            mismatches = sum(1 for a, b in zip(sequence[i:i+pattern_len], pattern) if a != b)
            if mismatches <= max_mismatches:
                return i
        return -1
    
    def convert_to_rna(self, sequence: str) -> str:
        """Convert DNA to RNA (T -> U)"""
        return sequence.replace('T', 'U')
    
    def prepare_for_folding(self, sequences: List[Dict]) -> List[Dict]:
        """
        Prepare sequences for structure prediction
        
        Args:
            sequences: List of sequence dictionaries with full_sequence
            
        Returns:
            Updated list with folding-ready sequences
        """
        for seq_dict in sequences:
            # Use full sequence for folding (not just random region)
            if 'full_sequence' in seq_dict:
                dna_seq = seq_dict['full_sequence']
                
                if self.molecule_type == "RNA":
                    seq_dict['folding_sequence'] = self.convert_to_rna(dna_seq)
                else:
                    seq_dict['folding_sequence'] = dna_seq
                    
        return sequences
    
    def get_stats(self, sequences: List[Dict]) -> Dict:
        """Get statistics about parsed sequences"""
        if not sequences:
            return {}
        
        total = len(sequences)
        with_primers = sum(1 for s in sequences if s.get('primers_found', False))
        
        lengths = [s.get('random_region_length', 0) for s in sequences 
                   if 'random_region_length' in s]
        
        return {
            'total_sequences': total,
            'sequences_with_primers': with_primers,
            'sequences_without_primers': total - with_primers,
            'mean_random_region_length': sum(lengths) / len(lengths) if lengths else 0,
            'min_random_region_length': min(lengths) if lengths else 0,
            'max_random_region_length': max(lengths) if lengths else 0
        }
    
    def find_identical_sequences(self, sequences: List[Dict]) -> List[Dict]:
        """
        Find groups of identical sequences
        
        Args:
            sequences: List of sequence dictionaries
            
        Returns:
            List of dictionaries describing identical sequence groups
        """
        from collections import defaultdict
        
        seq_groups = defaultdict(list)
        
        for seq_dict in sequences:
            seq = seq_dict.get('full_sequence', '')
            seq_groups[seq].append(seq_dict['id'])
        
        # Only keep groups with duplicates
        duplicates = []
        for seq, ids in seq_groups.items():
            if len(ids) > 1:
                duplicates.append({
                    'sequence': seq,
                    'count': len(ids),
                    'sequence_ids': ids
                })
        
        # Sort by count (most duplicates first)
        duplicates.sort(key=lambda x: -x['count'])
        
        return duplicates
    
    def find_similar_sequences(self, sequences: List[Dict], 
                               max_mismatches: int = 2) -> List[Dict]:
        """
        Find groups of highly similar sequences (allowing mismatches)
        
        Args:
            sequences: List of sequence dictionaries
            max_mismatches: Maximum number of mismatches allowed
            
        Returns:
            List of dictionaries describing similar sequence groups
        """
        similar_groups = []
        processed = set()
        
        for i, seq_dict1 in enumerate(sequences):
            if seq_dict1['id'] in processed:
                continue
                
            seq1 = seq_dict1.get('full_sequence', '')
            group = [seq_dict1['id']]
            
            for j, seq_dict2 in enumerate(sequences[i+1:], i+1):
                if seq_dict2['id'] in processed:
                    continue
                    
                seq2 = seq_dict2.get('full_sequence', '')
                
                # Only compare sequences of same length
                if len(seq1) == len(seq2):
                    mismatches = sum(1 for a, b in zip(seq1, seq2) if a != b)
                    
                    if mismatches <= max_mismatches:
                        group.append(seq_dict2['id'])
                        processed.add(seq_dict2['id'])
            
            if len(group) > 1:
                processed.add(seq_dict1['id'])
                similar_groups.append({
                    'representative': seq_dict1['id'],
                    'count': len(group),
                    'sequence_ids': group,
                    'sequence': seq1
                })
        
        # Sort by count
        similar_groups.sort(key=lambda x: -x['count'])
        
        return similar_groups
