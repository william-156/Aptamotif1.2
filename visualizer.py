"""
Visualization Module
Create heatmaps, sequence logos, and other visualizations
"""

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import logomaker
import plotly.graph_objects as go
import plotly.express as px
from typing import List, Dict, Tuple
import io
import base64


class MotifVisualizer:
    """Create visualizations for motif analysis"""
    
    def __init__(self, style: str = "seaborn-v0_8-darkgrid"):
        """
        Initialize visualizer
        
        Args:
            style: Matplotlib style
        """
        self.style = style
        plt.style.use('default')  # Use default style for compatibility
        
    def create_motif_heatmap(self, matrix: np.ndarray, motif_names: List[str],
                             sequence_names: List[str],
                             interactive: bool = True) -> go.Figure:
        """
        Create heatmap showing motif presence across sequences
        
        Args:
            matrix: Binary matrix (sequences x motifs)
            motif_names: List of motif names
            sequence_names: List of sequence names
            interactive: If True, return Plotly figure; else matplotlib
            
        Returns:
            Plotly or matplotlib figure
        """
        if interactive:
            fig = go.Figure(data=go.Heatmap(
                z=matrix,
                x=motif_names,
                y=sequence_names,
                colorscale=[[0, 'white'], [1, '#2E86AB']],
                showscale=True,
                hovertemplate='Sequence: %{y}<br>Motif: %{x}<br>Present: %{z}<extra></extra>'
            ))
            
            fig.update_layout(
                title='Motif Presence Across Sequences',
                xaxis_title='Motifs',
                yaxis_title='Sequences',
                height=max(400, len(sequence_names) * 20),
                width=max(600, len(motif_names) * 40),
                xaxis={'tickangle': -45}
            )
            
            return fig
        else:
            fig, ax = plt.subplots(figsize=(max(10, len(motif_names) * 0.5),
                                           max(8, len(sequence_names) * 0.3)))
            sns.heatmap(matrix, xticklabels=motif_names, yticklabels=sequence_names,
                       cmap='Blues', cbar_kws={'label': 'Presence'},
                       ax=ax, linewidths=0.5)
            ax.set_xlabel('Motifs')
            ax.set_ylabel('Sequences')
            ax.set_title('Motif Presence Across Sequences')
            plt.xticks(rotation=45, ha='right')
            plt.tight_layout()
            return fig
    
    def create_sequence_logo(self, sequences: List[str], title: str = "Sequence Logo") -> plt.Figure:
        """
        Create sequence logo from aligned sequences
        
        Args:
            sequences: List of sequence strings (should be aligned/same length)
            title: Plot title
            
        Returns:
            Matplotlib figure
        """
        # Filter to same length sequences
        if not sequences:
            return None
        
        seq_length = len(sequences[0])
        filtered_seqs = [s for s in sequences if len(s) == seq_length]
        
        if not filtered_seqs:
            return None
        
        # Create counts matrix
        counts_matrix = np.zeros((seq_length, 4))  # A, C, G, T/U
        base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'U': 3}
        
        for seq in filtered_seqs:
            for pos, base in enumerate(seq):
                if base in base_to_idx:
                    counts_matrix[pos, base_to_idx[base]] += 1
        
        # Convert to dataframe
        counts_df = pd.DataFrame(counts_matrix, columns=['A', 'C', 'G', 'T'])
        
        # Create information matrix
        info_df = logomaker.transform_matrix(counts_df,
                                             from_type='counts',
                                             to_type='information')
        
        # Create logo
        fig, ax = plt.subplots(figsize=(seq_length * 0.5, 3))
        logo = logomaker.Logo(info_df, ax=ax)
        ax.set_ylabel('Information (bits)')
        ax.set_xlabel('Position')
        ax.set_title(title)
        
        plt.tight_layout()
        return fig
    
    def create_enrichment_barplot(self, motifs: List[Dict], top_n: int = 20,
                                  interactive: bool = True):
        """
        Create bar plot of top enriched motifs
        
        Args:
            motifs: List of motif dictionaries with enrichment scores
            top_n: Number of top motifs to show
            interactive: If True, return Plotly figure
            
        Returns:
            Plotly or matplotlib figure
        """
        # Sort by fold enrichment
        sorted_motifs = sorted(motifs, key=lambda x: x.get('fold_enrichment', 0),
                              reverse=True)[:top_n]
        
        motif_names = [m['motif'] for m in sorted_motifs]
        enrichments = [m.get('fold_enrichment', 0) for m in sorted_motifs]
        pvalues = [m.get('adjusted_p_value', 1) for m in sorted_motifs]
        
        # Color by significance
        colors = ['red' if p < 0.05 else 'blue' for p in pvalues]
        
        if interactive:
            fig = go.Figure(data=[
                go.Bar(
                    x=enrichments,
                    y=motif_names,
                    orientation='h',
                    marker_color=colors,
                    hovertemplate='Motif: %{y}<br>Enrichment: %{x:.2f}x<br>Adj. p-value: %{text}<extra></extra>',
                    text=[f"{p:.2e}" for p in pvalues]
                )
            ])
            
            fig.update_layout(
                title=f'Top {top_n} Enriched Motifs',
                xaxis_title='Fold Enrichment',
                yaxis_title='Motif',
                height=max(400, len(motif_names) * 25),
                showlegend=False
            )
            
            return fig
        else:
            fig, ax = plt.subplots(figsize=(10, max(6, len(motif_names) * 0.3)))
            y_pos = np.arange(len(motif_names))
            ax.barh(y_pos, enrichments, color=colors)
            ax.set_yticks(y_pos)
            ax.set_yticklabels(motif_names)
            ax.set_xlabel('Fold Enrichment')
            ax.set_ylabel('Motif')
            ax.set_title(f'Top {top_n} Enriched Motifs')
            plt.tight_layout()
            return fig
    
    def create_pvalue_volcano_plot(self, motifs: List[Dict],
                                   interactive: bool = True):
        """
        Create volcano plot of motif enrichment vs significance
        
        Args:
            motifs: List of motif dictionaries
            interactive: If True, return Plotly figure
            
        Returns:
            Plotly or matplotlib figure
        """
        enrichments = [np.log2(m.get('fold_enrichment', 1)) for m in motifs]
        neg_log_pvals = [-np.log10(max(m.get('adjusted_p_value', 1), 1e-300))
                        for m in motifs]
        motif_names = [m['motif'] for m in motifs]
        
        # Determine significance
        significant = [m.get('adjusted_p_value', 1) < 0.05 for m in motifs]
        colors_list = ['red' if sig else 'gray' for sig in significant]
        
        if interactive:
            df = pd.DataFrame({
                'log2_enrichment': enrichments,
                'neg_log10_pvalue': neg_log_pvals,
                'motif': motif_names,
                'significant': significant
            })
            
            fig = px.scatter(df, x='log2_enrichment', y='neg_log10_pvalue',
                           color='significant',
                           hover_data=['motif'],
                           color_discrete_map={True: 'red', False: 'gray'},
                           labels={
                               'log2_enrichment': 'Log2(Fold Enrichment)',
                               'neg_log10_pvalue': '-Log10(Adjusted p-value)'
                           },
                           title='Motif Enrichment Volcano Plot')
            
            # Add significance threshold line
            fig.add_hline(y=-np.log10(0.05), line_dash="dash",
                         line_color="blue", annotation_text="p=0.05")
            
            return fig
        else:
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.scatter(enrichments, neg_log_pvals, c=colors_list, alpha=0.6)
            ax.axhline(-np.log10(0.05), color='blue', linestyle='--',
                      label='p=0.05')
            ax.set_xlabel('Log2(Fold Enrichment)')
            ax.set_ylabel('-Log10(Adjusted p-value)')
            ax.set_title('Motif Enrichment Volcano Plot')
            ax.legend()
            plt.tight_layout()
            return fig
    
    def create_gc_content_distribution(self, motifs: List[Dict],
                                       interactive: bool = True):
        """
        Create distribution plot of GC content in motifs
        
        Args:
            motifs: List of motif dictionaries
            interactive: If True, return Plotly figure
            
        Returns:
            Plotly or matplotlib figure
        """
        gc_contents = [m.get('gc_content', 0) * 100 for m in motifs]
        
        if interactive:
            fig = go.Figure(data=[go.Histogram(x=gc_contents, nbinsx=20)])
            fig.update_layout(
                title='GC Content Distribution of Motifs',
                xaxis_title='GC Content (%)',
                yaxis_title='Count',
                showlegend=False
            )
            return fig
        else:
            fig, ax = plt.subplots(figsize=(8, 5))
            ax.hist(gc_contents, bins=20, edgecolor='black')
            ax.set_xlabel('GC Content (%)')
            ax.set_ylabel('Count')
            ax.set_title('GC Content Distribution of Motifs')
            plt.tight_layout()
            return fig
    
    def create_structure_comparison(self, sequences: List[Dict],
                                   interactive: bool = True):
        """
        Create visualization comparing structure properties
        
        Args:
            sequences: List of sequence dicts with structure info
            interactive: If True, return Plotly figure
            
        Returns:
            Plotly or matplotlib figure
        """
        mfe_values = [s.get('mfe', 0) for s in sequences if 'mfe' in s]
        percent_paired = [s.get('structure_elements', {}).get('percent_paired', 0)
                         for s in sequences if 'structure_elements' in s]
        seq_names = [s['id'] for s in sequences if 'mfe' in s]
        
        if interactive:
            fig = go.Figure()
            
            fig.add_trace(go.Scatter(
                x=percent_paired,
                y=mfe_values,
                mode='markers',
                marker=dict(size=10, color='blue', opacity=0.6),
                text=seq_names,
                hovertemplate='%{text}<br>Paired: %{x:.1f}%<br>MFE: %{y:.2f} kcal/mol<extra></extra>'
            ))
            
            fig.update_layout(
                title='Secondary Structure Properties',
                xaxis_title='Percent Paired (%)',
                yaxis_title='Minimum Free Energy (kcal/mol)',
                height=500
            )
            
            return fig
        else:
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.scatter(percent_paired, mfe_values, alpha=0.6)
            ax.set_xlabel('Percent Paired (%)')
            ax.set_ylabel('Minimum Free Energy (kcal/mol)')
            ax.set_title('Secondary Structure Properties')
            plt.tight_layout()
            return fig
    
    def fig_to_base64(self, fig) -> str:
        """Convert matplotlib figure to base64 string"""
        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
        buf.seek(0)
        img_str = base64.b64encode(buf.read()).decode()
        plt.close(fig)
        return img_str
    
    def create_pairwise_similarity_heatmap(self, similarity_matrix: np.ndarray,
                                          sequence_names: List[str],
                                          interactive: bool = True) -> go.Figure:
        """
        Create heatmap showing number of shared significant motifs between sequence pairs
        
        Args:
            similarity_matrix: Matrix of shared motif counts (sequences x sequences)
            sequence_names: List of sequence names
            interactive: If True, return Plotly figure
            
        Returns:
            Plotly or matplotlib figure
        """
        if interactive:
            # Use rainbow/turbo colorscale for better visibility
            # Options: 'Turbo', 'Jet', 'Rainbow', 'Portland', 'Spectral'
            fig = go.Figure(data=go.Heatmap(
                z=similarity_matrix,
                x=sequence_names,
                y=sequence_names,
                colorscale='Turbo',  # Vibrant rainbow scale
                showscale=True,
                colorbar=dict(
                    title="Shared<br>Motifs",
                    titleside='right',
                    tickmode='linear',
                    tick0=0,
                    dtick=1
                ),
                hovertemplate='<b>%{y}</b> vs <b>%{x}</b><br>Shared motifs: <b>%{z}</b><extra></extra>',
                text=similarity_matrix,
                texttemplate='%{text}',
                textfont={"size": 10}
            ))
            
            fig.update_layout(
                title='Pairwise Sequence Similarity (Shared Significant Motifs)',
                xaxis_title='Sequence',
                yaxis_title='Sequence',
                height=max(600, len(sequence_names) * 25),
                width=max(600, len(sequence_names) * 25),
                xaxis={'tickangle': -45, 'side': 'bottom'},
                yaxis={'tickangle': 0}
            )
            
            return fig
        else:
            fig, ax = plt.subplots(figsize=(max(10, len(sequence_names) * 0.5),
                                           max(10, len(sequence_names) * 0.5)))
            
            # Use rainbow colormap for matplotlib
            sns.heatmap(similarity_matrix,
                       xticklabels=sequence_names,
                       yticklabels=sequence_names,
                       cmap='turbo',  # Rainbow colormap
                       cbar_kws={'label': 'Shared Motifs'},
                       ax=ax,
                       annot=True,
                       fmt='d',
                       linewidths=0.5,
                       linecolor='white')
            
            ax.set_xlabel('Sequence')
            ax.set_ylabel('Sequence')
            ax.set_title('Pairwise Sequence Similarity (Shared Significant Motifs)')
            plt.xticks(rotation=45, ha='right')
            plt.yticks(rotation=0)
            plt.tight_layout()
            return fig
