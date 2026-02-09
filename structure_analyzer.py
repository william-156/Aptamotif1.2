"""
Structure Analysis Module
RNA/DNA secondary structure prediction using ViennaRNA
"""

import RNA
from typing import List, Dict, Tuple
from collections import Counter
import re


class StructureAnalyzer:
    """Predict and analyze RNA/DNA secondary structures"""
    
    def __init__(self, molecule_type: str = "RNA", temperature: float = 37.0,
                 salt_concentration: float = 1.021, magnesium_concentration: float = 0.0):
        """
        Initialize structure analyzer
        
        Args:
            molecule_type: "RNA" or "ssDNA"
            temperature: Temperature in Celsius for folding (default: 37.0)
            salt_concentration: Monovalent salt concentration in mol/L (default: 1.021 M)
            magnesium_concentration: Divalent Mg2+ concentration in mol/L (default: 0.0 M)
        """
        self.molecule_type = molecule_type
        self.temperature = temperature
        self.salt_concentration = salt_concentration
        self.magnesium_concentration = magnesium_concentration
        
        # Set RNA parameters
        self.model = RNA.md()
        self.model.temperature = temperature
        
        # Set salt concentrations (in mol/L)
        # saltMLions: monovalent salt concentration
        # saltDPions: divalent salt concentration (Mg2+)
        try:
            self.model.saltMLions = salt_concentration
            self.model.saltDPions = magnesium_concentration
        except AttributeError:
            # Older ViennaRNA versions may not support salt parameters
            print(f"Warning: Salt parameters not supported in this ViennaRNA version")
        
        if molecule_type == "ssDNA":
            # For DNA, adjust parameters
            self.model.uniq_ML = 1
    
    def fold_sequence(self, sequence: str) -> Tuple[str, float]:
        """
        Fold a single sequence and return structure and MFE
        
        Args:
            sequence: RNA or DNA sequence
            
        Returns:
            Tuple of (dot-bracket structure, minimum free energy)
        """
        # Create fold compound
        fc = RNA.fold_compound(sequence, self.model)
        
        # Calculate MFE structure
        structure, mfe = fc.mfe()
        
        return structure, mfe
    
    def _build_vienna_fasta(self, seq_id: str, sequence: str, structure: str | None = None) -> str:
        # Sanitize header
        seq_id = re.sub(r"\s+", "_", seq_id.strip())

        # Normalize sequence
        sequence = sequence.upper().replace("T", "U")
        sequence = re.sub(r"[^ACGUN]", "", sequence)

        if structure:
            structure = structure.strip()
            if len(structure) != len(sequence):
                raise ValueError("Structure length must match sequence length")

            return f">{seq_id}\n{sequence}\n{structure}\n"

        return f">{seq_id}\n{sequence}\n"
    
    def generate_structure_svg(self, sequence: str, structure: str,
                               filename: str = None, layout: int = 1) -> str:
        """
        Multi-tiered structure visualization
        
        Tries multiple methods in priority order:
        1. RNAplot (best quality, requires system install)
        2. forgi (good quality, pure Python)
        3. Text fallback (always works)
        
        Args:
            sequence: Nucleotide sequence
            structure: Dot-bracket structure
            filename: Optional filename to save SVG
            layout: Layout type (for RNAplot)
            
        Returns:
            SVG content as string
        """
        # Method 1: Try RNAplot (works on HuggingFace, local machines)
        svg = self._try_rnaplot(sequence, structure, layout)
        if svg and '<svg' in svg.lower():
            print("✓ Using RNAplot (best quality)")
            if filename:
                with open(filename, 'w') as f:
                    f.write(svg)
            return svg
        
        # Method 2: Try forgi (works everywhere, pure Python)
        svg = self._try_forgi(sequence, structure)
        if svg and '<svg' in svg.lower():
            print("✓ Using forgi visualization (good quality)")
            if filename:
                with open(filename, 'w') as f:
                    f.write(svg)
            return svg
        
        # Method 3: Fallback to text-based SVG
        print("⚠ Using text fallback")
        return self._generate_fallback_svg(sequence, structure)
    
    def _try_rnaplot(self, sequence: str, structure: str, layout: int) -> str:
        """
        Try to use local RNAplot command
        
        Returns:
            SVG content or None if failed
        """
        import subprocess
        import tempfile
        import os
        import shutil
        import time
        
        if not shutil.which("RNAplot"):
            return None
        
        with tempfile.TemporaryDirectory() as tmpdir:
            try:
                # Method A: Try stdin approach
                input_data = f"{sequence}\n{structure}\n"
                
                result = subprocess.run(
                    ["RNAplot", "--output-format=svg", f"--layout-type={layout}"],
                    input=input_data,
                    cwd=tmpdir,
                    capture_output=True,
                    text=True,
                    timeout=15,
                )
                
                time.sleep(0.1)
                
                svg_files = [f for f in os.listdir(tmpdir) if f.endswith('.svg')]
                
                if svg_files:
                    with open(os.path.join(tmpdir, svg_files[0])) as f:
                        return f.read()
                
            except Exception as e:
                print(f"RNAplot method failed: {e}")
        
        return None
    
    def _try_forgi(self, sequence: str, structure: str) -> str:
        """
        Try to use forgi library for visualization
        
        Returns:
            SVG content or None if failed
        """
        try:
            import matplotlib
            matplotlib.use('Agg')  # Non-interactive backend
            import matplotlib.pyplot as plt
            import forgi.visual.mplotlib as fvm
            import forgi.graph.bulge_graph as fgb
            from io import BytesIO
            
            # Create RNA structure graph from dot-bracket
            bg = fgb.BulgeGraph.from_dotbracket(structure, seq=sequence)
            
            # Create figure with good size
            fig = plt.figure(figsize=(10, 10))
            ax = plt.gca()
            
            # Plot RNA structure
            fvm.plot_rna(bg, ax=ax, 
                        text_kwargs={"fontweight": "bold", "fontsize": 8},
                        lighten=0.7,
                        backbone_kwargs={"linewidth": 3})
            
            # Clean up axes
            ax.set_axis_off()
            ax.set_aspect('equal')
            plt.tight_layout()
            
            # Save to SVG string
            buf = BytesIO()
            plt.savefig(buf, format='svg', bbox_inches='tight', 
                       pad_inches=0.1, dpi=150)
            plt.close(fig)
            
            # Get SVG content
            buf.seek(0)
            svg_content = buf.read().decode('utf-8')
            
            return svg_content
            
        except ImportError:
            print("forgi not installed (pip install forgi)")
            return None
        except Exception as e:
            print(f"forgi visualization failed: {e}")
            import traceback
            traceback.print_exc()
            return None
    
    def _generate_svg_via_ps2svg(self, sequence: str, structure: str) -> str:
        """
        Legacy fallback method - no longer used
        """
        return self._generate_fallback_svg(sequence, structure)
    
    def test_svg_capability(self) -> dict:
        """
        Test which SVG generation method is available
        
        Returns:
            Dictionary with capability info
        """
        import subprocess
        import shutil
        
        result = {
            'rnaplot_available': False,
            'rnaplot_path': None,
            'method': 'fallback',
            'message': ''
        }
        
        # Test if RNAplot command is available
        rnaplot_path = shutil.which('RNAplot')
        
        if rnaplot_path:
            result['rnaplot_available'] = True
            result['rnaplot_path'] = rnaplot_path
            
            # Try to run it to verify it works
            try:
                test_result = subprocess.run(
                    ['RNAplot', '--help'],
                    capture_output=True,
                    timeout=5
                )
                
                if test_result.returncode == 0:
                    result['method'] = 'RNAplot_command'
                    result['message'] = f'✅ RNAplot found at: {rnaplot_path}'
                    return result
                else:
                    result['message'] = f'⚠️ RNAplot found but may not work properly'
                    return result
                    
            except Exception as e:
                result['message'] = f'⚠️ RNAplot found but failed to execute: {str(e)}'
                return result
        
        # If RNAplot not found
        result['message'] = '❌ RNAplot not found. Install ViennaRNA package.'
        result['install_help'] = {
            'macos': 'brew install viennarna',
            'ubuntu': 'sudo apt-get install vienna-rna',
            'conda': 'conda install -c bioconda viennarna'
        }
        
        return result
    
    def _generate_fallback_svg(self, sequence: str, structure: str) -> str:
        """
        Generate a simple text-based SVG if other methods fail
        """
        # Create simple SVG with structure text
        svg_width = max(600, len(sequence) * 8)
        svg_height = 300
        
        svg = f'''<svg width="{svg_width}" height="{svg_height}" xmlns="http://www.w3.org/2000/svg">
            <rect width="100%" height="100%" fill="white"/>
            <text x="10" y="30" font-family="monospace" font-size="12" fill="black">Sequence:</text>
            <text x="10" y="50" font-family="monospace" font-size="11" fill="blue">{sequence}</text>
            <text x="10" y="80" font-family="monospace" font-size="12" fill="black">Structure:</text>
            <text x="10" y="100" font-family="monospace" font-size="11" fill="green">{structure}</text>
            <text x="10" y="130" font-family="monospace" font-size="10" fill="gray" style="font-style: italic;">
                Note: Install ViennaRNA with RNAplot for graphical structure diagrams
            </text>
            <text x="10" y="150" font-family="monospace" font-size="10" fill="gray" style="font-style: italic;">
                Run: sudo apt-get install vienna-rna (Ubuntu) or brew install viennarna (macOS)
            </text>
        </svg>'''
        
        return svg
    
    def generate_all_structure_svgs(self, sequences: List[Dict],
                                    output_dir: str) -> List[str]:
        """
        Generate SVG files for all sequences
        
        Args:
            sequences: List of sequence dicts with structures
            output_dir: Directory to save SVG files
            
        Returns:
            List of generated SVG file paths
        """
        import os
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        svg_files = []
        
        for seq_dict in sequences:
            if 'structure' not in seq_dict or 'folding_sequence' not in seq_dict:
                continue
            
            seq_id = seq_dict['id']
            safe_id = seq_id.replace('/', '_').replace(' ', '_')
            
            svg_filename = os.path.join(output_dir, f"{safe_id}_structure.svg")
            
            svg_content = self.generate_structure_svg(
                seq_dict['folding_sequence'],
                seq_dict['structure'],
                svg_filename
            )
            
            if svg_content:
                svg_files.append(svg_filename)
                seq_dict['structure_svg'] = svg_filename
        
        return svg_files
    
    def fold_all_sequences(self, sequences: List[Dict]) -> List[Dict]:
        """
        Fold all sequences and add structure information
        
        Args:
            sequences: List of sequence dictionaries with 'folding_sequence'
            
        Returns:
            Updated sequence list with structure predictions
        """
        for seq_dict in sequences:
            if 'folding_sequence' in seq_dict:
                seq = seq_dict['folding_sequence']
                structure, mfe = self.fold_sequence(seq)
                
                seq_dict['structure'] = structure
                seq_dict['mfe'] = mfe
                seq_dict['structure_elements'] = self._parse_structure(structure)
                
        return sequences
    
    def _parse_structure(self, structure: str) -> Dict:
        """
        Parse dot-bracket structure to identify elements
        
        Args:
            structure: Dot-bracket notation string
            
        Returns:
            Dictionary with structure element counts
        """
        # Count structural elements
        unpaired = structure.count('.')
        paired = structure.count('(') + structure.count(')')
        
        # Identify hairpin loops (regions like "(((...)))")
        hairpins = len(re.findall(r'\(+\.+\)+', structure))
        
        # Identify bulges and internal loops
        bulges = len(re.findall(r'\(+\.+\(+', structure)) + \
                 len(re.findall(r'\)+\.+\)+', structure))
        
        return {
            'unpaired_bases': unpaired,
            'paired_bases': paired,
            'total_bases': len(structure),
            'percent_paired': (paired / len(structure) * 100) if len(structure) > 0 else 0,
            'hairpin_loops': hairpins,
            'bulges_internal_loops': bulges
        }
    
    def calculate_ensemble_diversity(self, sequence: str, num_samples: int = 100) -> float:
        """
        Calculate structural ensemble diversity
        
        Args:
            sequence: RNA/DNA sequence
            num_samples: Number of structures to sample
            
        Returns:
            Ensemble diversity score
        """
        fc = RNA.fold_compound(sequence, self.model)
        
        # Sample structures from ensemble
        structures = []
        for _ in range(num_samples):
            structure = fc.pbacktrack()
            structures.append(structure)
        
        # Calculate diversity (number of unique structures / total)
        unique_structures = len(set(structures))
        diversity = unique_structures / num_samples
        
        return diversity
    
    def find_consensus_structure(self, sequences: List[Dict]) -> Dict:
        """
        Find consensus structural elements across sequences
        
        Args:
            sequences: List of sequence dicts with structure predictions
            
        Returns:
            Dictionary with consensus structure information
        """
        if not sequences:
            return {}
        
        all_structures = [s['structure'] for s in sequences if 'structure' in s]
        
        if not all_structures:
            return {}
        
        # Analyze structure elements
        element_counts = {
            'hairpins': [],
            'percent_paired': [],
            'mfe_values': []
        }
        
        for seq_dict in sequences:
            if 'structure_elements' in seq_dict:
                elements = seq_dict['structure_elements']
                element_counts['hairpins'].append(elements.get('hairpin_loops', 0))
                element_counts['percent_paired'].append(elements.get('percent_paired', 0))
            if 'mfe' in seq_dict:
                element_counts['mfe_values'].append(seq_dict['mfe'])
        
        # Calculate consensus statistics
        consensus = {
            'mean_hairpins': sum(element_counts['hairpins']) / len(element_counts['hairpins'])
                           if element_counts['hairpins'] else 0,
            'mean_percent_paired': sum(element_counts['percent_paired']) / len(element_counts['percent_paired'])
                                 if element_counts['percent_paired'] else 0,
            'mean_mfe': sum(element_counts['mfe_values']) / len(element_counts['mfe_values'])
                       if element_counts['mfe_values'] else 0,
            'min_mfe': min(element_counts['mfe_values']) if element_counts['mfe_values'] else 0,
            'max_mfe': max(element_counts['mfe_values']) if element_counts['mfe_values'] else 0
        }
        
        return consensus
    
    def identify_structural_motifs(self, sequences: List[Dict]) -> List[Dict]:
        """
        Identify common structural motifs (e.g., hairpin loops with specific sequences)
        
        Args:
            sequences: List of sequence dicts with structures
            
        Returns:
            List of structural motif dictionaries
        """
        structural_motifs = []
        
        # Look for common hairpin loop sequences
        hairpin_sequences = []
        
        for seq_dict in sequences:
            if 'structure' not in seq_dict or 'folding_sequence' not in seq_dict:
                continue
                
            structure = seq_dict['structure']
            sequence = seq_dict['folding_sequence']
            
            # Find hairpin loops (pattern: opening stem, unpaired loop, closing stem)
            # Look for patterns like "(((...)))"
            i = 0
            while i < len(structure):
                if structure[i] == '(':
                    # Count opening brackets
                    open_count = 0
                    j = i
                    while j < len(structure) and structure[j] == '(':
                        open_count += 1
                        j += 1
                    
                    # Look for unpaired region
                    loop_start = j
                    while j < len(structure) and structure[j] == '.':
                        j += 1
                    loop_end = j
                    
                    # Check for closing brackets
                    close_count = 0
                    while j < len(structure) and structure[j] == ')':
                        close_count += 1
                        j += 1
                    
                    # If we found a hairpin (balanced brackets with loop)
                    if open_count >= 2 and close_count >= 2 and loop_end > loop_start:
                        loop_seq = sequence[loop_start:loop_end]
                        if len(loop_seq) >= 3:  # Minimum loop size
                            hairpin_sequences.append({
                                'loop_sequence': loop_seq,
                                'stem_length': min(open_count, close_count),
                                'loop_length': loop_end - loop_start,
                                'seq_id': seq_dict['id']
                            })
                    
                    i = j
                else:
                    i += 1
        
        # Count recurring hairpin loop sequences
        loop_counter = Counter([h['loop_sequence'] for h in hairpin_sequences])
        
        for loop_seq, count in loop_counter.most_common(20):  # Top 20
            if count >= 2:  # Appears in at least 2 sequences
                # Get all instances
                instances = [h for h in hairpin_sequences if h['loop_sequence'] == loop_seq]
                
                structural_motifs.append({
                    'type': 'hairpin_loop',
                    'loop_sequence': loop_seq,
                    'count': count,
                    'mean_stem_length': sum(h['stem_length'] for h in instances) / len(instances),
                    'sequences': [h['seq_id'] for h in instances]
                })
        
        return structural_motifs
    
    def generate_structure_annotation(self, sequence: str, structure: str,
                                       mfe: float) -> str:
        """
        Generate a formatted structure annotation
        
        Args:
            sequence: Nucleotide sequence
            structure: Dot-bracket structure
            mfe: Minimum free energy
            
        Returns:
            Formatted string with sequence, structure, and MFE
        """
        annotation = f"Sequence:  {sequence}\n"
        annotation += f"Structure: {structure}\n"
        annotation += f"MFE:       {mfe:.2f} kcal/mol\n"
        
        return annotation
    
    def export_structures(self, sequences: List[Dict], output_format: str = "text") -> str:
        """
        Export all structures in specified format
        
        Args:
            sequences: List of sequence dicts with structures
            output_format: "text", "fasta", or "ct"
            
        Returns:
            Formatted structure data
        """
        output = []
        
        for seq_dict in sequences:
            if 'structure' not in seq_dict:
                continue
            
            seq_id = seq_dict['id']
            seq = seq_dict.get('folding_sequence', '')
            structure = seq_dict['structure']
            mfe = seq_dict.get('mfe', 0)
            
            if output_format == "text":
                output.append(f">{seq_id}")
                output.append(self.generate_structure_annotation(seq, structure, mfe))
                output.append("")
            
            elif output_format == "fasta":
                output.append(f">{seq_id} MFE={mfe:.2f}")
                output.append(seq)
                output.append(structure)
                output.append("")
        
        return '\n'.join(output)
