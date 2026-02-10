"""
AptaMotif - Aptamer Motif Analysis Tool
Main Streamlit Web Application
"""

import streamlit as st
import pandas as pd
import numpy as np
from sequence_parser import SequenceParser
from motif_finder import MotifFinder
from motif_statistics import MotifStatistics
from structure_analyzer import StructureAnalyzer
from visualizer import MotifVisualizer
import io

# Page configuration
st.set_page_config(
    page_title="AptaMotif",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Initialize session state
if 'analysis_complete' not in st.session_state:
    st.session_state.analysis_complete = False
if 'results' not in st.session_state:
    st.session_state.results = {}

# Custom CSS
st.markdown("""
    <style>
    .main-header {
        font-size: 2.5rem;
        font-weight: bold;
        color: #2E86AB;
        text-align: center;
        margin-bottom: 0.5rem;
    }
    .sub-header {
        font-size: 1.2rem;
        color: #666;
        text-align: center;
        margin-bottom: 2rem;
    }
    .stAlert {
        margin-top: 1rem;
    }
    </style>
""", unsafe_allow_html=True)

# Header
st.markdown('<div class="main-header">🧬 AptaMotif</div>', unsafe_allow_html=True)
st.markdown('<div class="sub-header">Aptamer Motif Discovery & Enrichment Analysis</div>',
           unsafe_allow_html=True)

# Sidebar - Configuration
st.sidebar.header("⚙️ Configuration")

# Molecule type
molecule_type = st.sidebar.selectbox(
    "Molecule Type",
    options=["RNA", "ssDNA"],
    help="Select RNA for RNA aptamers or ssDNA for DNA aptamers"
)

# Primer configuration
st.sidebar.subheader("Primer Configuration")

with st.sidebar.expander("📝 Edit Primers", expanded=False):
    forward_primer = st.text_input(
        "Forward Primer (5' → 3')",
        value="GGGAGATACCAGCTTATTCAATT",
        help="Enter the forward primer sequence"
    )
    
    reverse_primer_region = st.text_input(
        "Reverse Primer Region (appears after N-region)",
        value="AGATAGTAAGTGCAATCT",
        help="The sequence that appears after the random region"
    )
    
    n_region_length = st.number_input(
        "Expected N-Region Length",
        min_value=10,
        max_value=200,
        value=71,
        help="Expected length of the random region"
    )

# Analysis parameters
st.sidebar.subheader("Analysis Parameters")

min_motif_length = st.sidebar.slider(
    "Minimum Motif Length (bp)",
    min_value=4,
    max_value=20,
    value=5,
    help="Minimum length of motifs to search for"
)

max_motif_length = st.sidebar.slider(
    "Maximum Motif Length (bp)",
    min_value=min_motif_length,
    max_value=25,
    value=15,
    help="Maximum length of motifs to search for"
)

min_sequences = st.sidebar.slider(
    "Minimum Sequences",
    min_value=2,
    max_value=20,
    value=2,
    help="Minimum number of sequences a motif must appear in"
)

fdr_threshold = st.sidebar.slider(
    "FDR Threshold",
    min_value=0.01,
    max_value=0.10,
    value=0.05,
    step=0.01,
    help="False Discovery Rate threshold for significance"
)

# Structure analysis options
st.sidebar.subheader("Structure Analysis")
perform_structure = st.sidebar.checkbox(
    "Perform Secondary Structure Analysis",
    value=True,
    help="Predict secondary structures using ViennaRNA (may take longer)"
)

if perform_structure:
    with st.sidebar.expander("⚙️ Folding Parameters", expanded=False):
        folding_temperature = st.number_input(
            "Temperature (°C)",
            min_value=0.0,
            max_value=100.0,
            value=37.0,
            step=1.0,
            help="Temperature for RNA/DNA folding (default: 37°C)"
        )
        
        salt_concentration = st.number_input(
            "Monovalent Salt [Na+] (M)",
            min_value=0.0,
            max_value=5.0,
            value=1.021,
            step=0.1,
            format="%.3f",
            help="Monovalent salt concentration in mol/L (default: 1.021 M ≈ physiological)"
        )
        
        magnesium_concentration = st.number_input(
            "Divalent Salt [Mg2+] (mM)",
            min_value=0.0,
            max_value=100.0,
            value=0.0,
            step=0.5,
            format="%.1f",
            help="Magnesium concentration in mmol/L (default: 0 mM)"
        )
        
        st.caption("💡 Typical conditions:")
        st.caption("• Physiological: 37°C, 150 mM Na+, 2 mM Mg2+")
        st.caption("• SELEX buffer: 25°C, variable salts")
        st.caption("• In vitro transcription: 37°C, 1 M Na+")
else:
    folding_temperature = 37.0
    salt_concentration = 1.021
    magnesium_concentration = 0.0

# Main content area
tab1, tab2, tab3 = st.tabs(["📊 Input & Analysis", "📈 Results", "ℹ️ Help"])

with tab1:
    st.header("Sequence Input")
    
    col1, col2 = st.columns([2, 1])
    
    with col1:
        input_method = st.radio(
            "Input Method",
            options=["Paste Sequences", "Upload File"],
            horizontal=True
        )
        
        if input_method == "Paste Sequences":
            sequence_input = st.text_area(
                "Enter Sequences",
                height=300,
                key="sequence_text_area",
                placeholder="Paste sequences here (one per line or FASTA format)...\n\n"
                           "Example:\n"
                           ">Seq1\n"
                           "TTCTAATACGACTCACTATAGGGAGATACCAGCTTATTCAATT"
                           "GGATCCGGATCCGGATCCAGATAGTAAGTGCAATCT\n"
                           ">Seq2\n"
                           "TTCTAATACGACTCACTATAGGGAGATACCAGCTTATTCAATT"
                           "CCAATTCCAATTCCAATTCCAGATAGTAAGTGCAATCT",
                help="Enter sequences in FASTA format or plain text (one per line)"
            )
        else:
            uploaded_file = st.file_uploader(
                "Upload Sequence File",
                type=['txt', 'fasta', 'fa'],
                help="Upload a text or FASTA file with sequences"
            )
            
            if uploaded_file is not None:
                # Store uploaded content in session state
                sequence_input = uploaded_file.read().decode('utf-8')
                st.session_state.uploaded_sequences = sequence_input
            elif 'uploaded_sequences' in st.session_state:
                # Use previously uploaded content
                sequence_input = st.session_state.uploaded_sequences
            else:
                sequence_input = ""
    
    with col2:
        st.info("""
        **Input Requirements:**
        - FASTA or plain text format
        - One sequence per line
        - Sequences should include primer regions
        - Minimum 2 sequences recommended
        
        **What happens:**
        1. Primers are identified and trimmed
        2. Random regions are extracted
        3. Motifs are discovered
        4. Statistical significance is calculated
        5. Structures are predicted (optional)
        """)
    
    # Analysis button
    st.markdown("---")
    col_btn1, col_btn2, col_btn3 = st.columns([1, 1, 2])
    
    with col_btn1:
        analyze_button = st.button("🔬 Run Analysis", type="primary", use_container_width=True)
    
    with col_btn2:
        if st.session_state.analysis_complete:
            clear_button = st.button("🗑️ Clear Results", use_container_width=True)
            if clear_button:
                st.session_state.analysis_complete = False
                st.session_state.results = {}
                st.rerun()
    
    # Run analysis
    if analyze_button and sequence_input:
        with st.spinner("🔬 Analyzing sequences..."):
            try:
                # Initialize components
                parser = SequenceParser(
                    forward_primer=forward_primer,
                    reverse_primer_region=reverse_primer_region,
                    molecule_type=molecule_type
                )
                
                motif_finder = MotifFinder(
                    min_length=min_motif_length,
                    max_length=max_motif_length,
                    min_sequences=min_sequences
                )
                
                stats_analyzer = MotifStatistics()
                visualizer = MotifVisualizer()
                
                # Step 1: Parse sequences
                st.info("📖 Parsing sequences...")
                sequences = parser.parse_sequences(sequence_input)
                sequences = parser.extract_random_region(sequences)
                sequences = parser.prepare_for_folding(sequences)
                
                seq_stats = parser.get_stats(sequences)
                
                # Check for identical/similar sequences
                st.info("🔍 Checking for duplicate sequences...")
                identical_seqs = parser.find_identical_sequences(sequences)
                similar_seqs = parser.find_similar_sequences(sequences, max_mismatches=2)
                
                # Step 2: Find motifs
                st.info("🔍 Discovering motifs...")
                motifs = motif_finder.find_enriched_motifs(sequences)
                
                if not motifs:
                    st.error("No motifs found! Try adjusting parameters.")
                    st.stop()
                
                # Step 3: Calculate statistics
                st.info("📊 Calculating statistical significance...")
                motifs = stats_analyzer.add_pvalues_to_motifs(
                    motifs,
                    seq_stats['mean_random_region_length']
                )
                motifs = stats_analyzer.add_enrichment_scores(
                    motifs,
                    seq_stats['mean_random_region_length']
                )
                motifs = stats_analyzer.apply_fdr_correction(motifs, method='fdr_bh')
                
                # Find merged motifs
                st.info("🔗 Identifying merged motif patterns...")
                merged_motifs = motif_finder.find_merged_motifs(motifs, sequences)
                
                # Filter redundant motifs
                st.info("🧹 Filtering redundant motifs...")
                filtered_motifs = motif_finder.filter_redundant_motifs(motifs)
                
                # Filter by FDR threshold
                significant_motifs = [m for m in filtered_motifs
                                    if m.get('adjusted_p_value', 1) < fdr_threshold]
                
                # Step 4: Structure analysis (optional)
                structures = None
                structural_motifs = None
                consensus_structure = None
                structure_svgs = []
                
                if perform_structure:
                    st.info("🧬 Predicting secondary structures...")
                    structure_analyzer = StructureAnalyzer(
                        molecule_type=molecule_type,
                        temperature=folding_temperature,
                        salt_concentration=salt_concentration / 1000.0 if salt_concentration > 10 else salt_concentration,  # Convert mM to M if needed
                        magnesium_concentration=magnesium_concentration / 1000.0  # Convert mM to M
                    )
                    sequences = structure_analyzer.fold_all_sequences(sequences)
                    consensus_structure = structure_analyzer.find_consensus_structure(sequences)
                    structural_motifs = structure_analyzer.identify_structural_motifs(sequences)
                    
                    structures = sequences
                
                # Step 5: Create visualizations
                st.info("📈 Generating visualizations...")
                
                # Motif matrix for heatmap (use all filtered motifs, not just top 50)
                matrix, motif_names, seq_names = motif_finder.create_motif_matrix(
                    filtered_motifs,
                    sequences
                )
                
                # Pairwise similarity matrix
                pairwise_matrix, pairwise_names = motif_finder.create_pairwise_similarity_matrix(
                    significant_motifs,
                    sequences
                )
                
                # Store results in session state
                st.session_state.results = {
                    'sequences': sequences,
                    'seq_stats': seq_stats,
                    'motifs': filtered_motifs,  # All filtered motifs
                    'significant_motifs': significant_motifs,
                    'merged_motifs': merged_motifs,
                    'identical_seqs': identical_seqs,
                    'similar_seqs': similar_seqs,
                    'matrix': matrix,
                    'motif_names': motif_names,
                    'seq_names': seq_names,
                    'pairwise_matrix': pairwise_matrix,
                    'pairwise_names': pairwise_names,
                    'structures': structures,
                    'structural_motifs': structural_motifs,
                    'consensus_structure': consensus_structure,
                    'molecule_type': molecule_type,  # Store for on-demand SVG generation
                    'folding_temperature': folding_temperature,
                    'salt_concentration': salt_concentration,
                    'magnesium_concentration': magnesium_concentration
                }
                
                st.session_state.analysis_complete = True
                st.success("✅ Analysis complete! Switch to the Results tab to view findings.")
                
            except Exception as e:
                st.error(f"❌ Error during analysis: {str(e)}")
                st.exception(e)
    
    elif analyze_button and not sequence_input:
        st.warning("⚠️ Please enter or upload sequences first!")

with tab2:
    if not st.session_state.analysis_complete:
        st.info("👈 Please run an analysis first using the Input & Analysis tab.")
    else:
        results = st.session_state.results
        
        # Summary statistics
        st.header("📊 Analysis Summary")
        
        col1, col2, col3, col4 = st.columns(4)
        
        with col1:
            st.metric("Total Sequences", results['seq_stats']['total_sequences'])
        with col2:
            st.metric("Unique Motifs", len(results['motifs']))
        with col3:
            st.metric("Significant Motifs", len(results['significant_motifs']))
        with col4:
            st.metric("Mean Region Length",
                     f"{results['seq_stats']['mean_random_region_length']:.1f} bp")
        
        # Sequence identity warnings
        if results.get('identical_seqs') or results.get('similar_seqs'):
            st.warning("⚠️ Duplicate or highly similar sequences detected!")
            
            if results.get('identical_seqs'):
                with st.expander("🔴 Identical Sequences", expanded=True):
                    for group in results['identical_seqs']:
                        st.write(f"**{group['count']} identical copies:** {', '.join(group['sequence_ids'])}")
            
            if results.get('similar_seqs'):
                with st.expander("🟡 Highly Similar Sequences (≤2 mismatches)"):
                    for group in results['similar_seqs']:
                        st.write(f"**{group['count']} similar sequences:** {', '.join(group['sequence_ids'])}")
        
        # Merged motifs info
        if results.get('merged_motifs'):
            st.info(f"ℹ️ Found {len(results['merged_motifs'])} potential merged motif patterns")
            with st.expander("🔗 View Merged Motifs"):
                merged_df = pd.DataFrame([{
                    'Merged Motif': m['motif'],
                    'Length': m['length'],
                    'Sequences': m['num_sequences'],
                    'Merged From': ' + '.join(m['merged_from'])
                } for m in results['merged_motifs']])
                st.dataframe(merged_df, use_container_width=True)
        
        st.markdown("---")
        
        # Motif table
        st.header("🔬 Enriched Motifs")
        
        # Add filter options
        filter_col1, filter_col2 = st.columns(2)
        with filter_col1:
            show_only_significant = st.checkbox("Show only significant motifs (adj. p < 0.05)", value=False)
        with filter_col2:
            min_enrichment = st.slider("Minimum fold enrichment", 1.0, 10.0, 1.0, 0.5)
        
        # Filter motifs based on selection
        display_motifs = results['motifs']
        if show_only_significant:
            display_motifs = results['significant_motifs']
        display_motifs = [m for m in display_motifs if m.get('fold_enrichment', 0) >= min_enrichment]
        
        # Create DataFrame with sequence IDs
        motif_df = pd.DataFrame([{
            'Motif': m['motif'],
            'Length': m['length'],
            'Sequences': m['num_sequences'],
            'Frequency': f"{m['frequency']:.2%}",
            'Enrichment': f"{m.get('fold_enrichment', 0):.2f}x",
            'p-value': f"{m.get('p_value', 1):.2e}",
            'Adj. p-value': f"{m.get('adjusted_p_value', 1):.2e}",
            'GC%': f"{m.get('gc_content', 0) * 100:.1f}%",
            'Sequence IDs': ', '.join([results['sequences'][i]['id'] for i in m['sequence_indices']
                                       if i < len(results['sequences'])])
        } for m in display_motifs])
        
        st.dataframe(motif_df, use_container_width=True, height=600)
        
        st.caption(f"Showing {len(display_motifs)} motifs")
        
        # Download button for motifs
        csv = motif_df.to_csv(index=False)
        st.download_button(
            label="📥 Download Motif Table (CSV)",
            data=csv,
            file_name="motif_enrichment.csv",
            mime="text/csv"
        )
        
        st.markdown("---")
        
        # Visualizations
        st.header("📈 Visualizations")
        
        viz_tabs = st.tabs(["Motif Heatmap", "Pairwise Similarity", "Enrichment Plot",
                           "Volcano Plot", "GC Distribution", "Sequence Logos"])
        
        visualizer = MotifVisualizer()
        
        with viz_tabs[0]:
            st.subheader("Motif Presence Heatmap")
            st.caption("Shows which motifs appear in which sequences")
            fig_heatmap = visualizer.create_motif_heatmap(
                results['matrix'],
                results['motif_names'],
                results['seq_names'],
                interactive=True
            )
            st.plotly_chart(fig_heatmap, use_container_width=True)
        
        with viz_tabs[1]:
            st.subheader("Pairwise Sequence Similarity")
            st.caption("Number of significant motifs shared between each pair of sequences")
            
            # Colorscale selector
            col_color1, col_color2 = st.columns([1, 3])
            with col_color1:
                colorscale_option = st.selectbox(
                    "Color Scale",
                    options=['Turbo', 'Jet', 'Rainbow', 'Portland', 'Spectral', 'Viridis', 'Plasma', 'Hot'],
                    index=0,
                    help="Choose colorscale for better visibility"
                )
            
            with col_color2:
                st.caption("**Color Guide:** Blue/Purple = Few shared motifs, Red/Yellow = Many shared motifs")
            
            # Create heatmap with selected colorscale
            import plotly.graph_objects as go
            
            fig_pairwise = go.Figure(data=go.Heatmap(
                z=results['pairwise_matrix'],
                x=results['pairwise_names'],
                y=results['pairwise_names'],
                colorscale=colorscale_option,
                showscale=True,
                colorbar=dict(
                title=dict(
                text="Shared<br>Motifs",
                side="right"   # ← THIS is valid
                ),
                tickmode='linear',
                tick0=0,
                dtick=1
                ),
                hovertemplate='<b>%{y}</b> vs <b>%{x}</b><br>Shared motifs: <b>%{z}</b><extra></extra>',
                text=results['pairwise_matrix'],
                texttemplate='%{text}',
                textfont={"size": 10}
            ))
            
            fig_pairwise.update_layout(
                title='Pairwise Sequence Similarity (Shared Significant Motifs)',
                xaxis_title='Sequence',
                yaxis_title='Sequence',
                height=max(600, len(results['pairwise_names']) * 25),
                width=max(600, len(results['pairwise_names']) * 25),
                xaxis={'tickangle': -45, 'side': 'bottom'},
                yaxis={'tickangle': 0}
            )
            
            st.plotly_chart(fig_pairwise, use_container_width=True)
        
        with viz_tabs[2]:
            st.subheader("Top Enriched Motifs")
            num_to_show = st.slider("Number of motifs to display", 10, 50, 20)
            fig_enrichment = visualizer.create_enrichment_barplot(
                display_motifs,
                top_n=num_to_show,
                interactive=True
            )
            st.plotly_chart(fig_enrichment, use_container_width=True)
        
        with viz_tabs[3]:
            st.subheader("Volcano Plot")
            st.caption("Fold enrichment vs. statistical significance")
            fig_volcano = visualizer.create_pvalue_volcano_plot(
                results['motifs'],
                interactive=True
            )
            st.plotly_chart(fig_volcano, use_container_width=True)
        
        with viz_tabs[4]:
            st.subheader("GC Content Distribution")
            fig_gc = visualizer.create_gc_content_distribution(
                results['motifs'],
                interactive=True
            )
            st.plotly_chart(fig_gc, use_container_width=True)
        
        with viz_tabs[5]:
            st.subheader("Sequence Logos for Top Motifs")
            
            # Show logos for significant motifs or top motifs
            logo_motifs = results['significant_motifs'][:20] if results['significant_motifs'] else results['motifs'][:20]
            
            # Create expandable sections for each motif
            for i, motif_dict in enumerate(logo_motifs, 1):
                motif_seq = motif_dict['motif']
                seq_ids = ', '.join([results['sequences'][idx]['id'] for idx in motif_dict['sequence_indices']
                                    if idx < len(results['sequences'])])
                
                with st.expander(f"{i}. {motif_seq} (in {motif_dict['num_sequences']} sequences: {seq_ids})"):
                    try:
                        fig_logo = visualizer.create_sequence_logo(
                            [motif_seq],
                            title=f"Motif: {motif_seq}"
                        )
                        if fig_logo:
                            st.pyplot(fig_logo)
                    except Exception as e:
                        st.warning(f"Could not generate logo: {str(e)}")
        
        # Structure analysis results
        if results['structures'] is not None:
            st.markdown("---")
            st.header("🧬 Secondary Structure Analysis")

            struct_col1, struct_col2 = st.columns(2)
            
            with struct_col1:
                st.subheader("Consensus Structure Properties")
                consensus = results['consensus_structure']
                
                # Display folding conditions
                temp = results.get('folding_temperature', 37.0)
                salt = results.get('salt_concentration', 1.021)
                mg = results.get('magnesium_concentration', 0.0)
                
                st.info(f"**Folding Conditions:**\n\n"
                       f"Temperature: {temp}°C\n\n"
                       f"[Na+]: {salt*1000 if salt < 10 else salt:.0f} mM\n\n"
                       f"[Mg2+]: {mg*1000 if mg < 1 else mg:.1f} mM")
                
                if consensus:
                    st.metric("Mean % Paired", f"{consensus['mean_percent_paired']:.1f}%")
                    st.metric("Mean MFE", f"{consensus['mean_mfe']:.2f} kcal/mol")
                    st.metric("Mean Hairpins", f"{consensus['mean_hairpins']:.1f}")
            
            with struct_col2:
                st.subheader("Structure Comparison")
                fig_struct = visualizer.create_structure_comparison(
                    results['structures'],
                    interactive=True
                )
                st.plotly_chart(fig_struct, use_container_width=True)
            
            # Structural motifs
            if results['structural_motifs']:
                st.subheader("Structural Motifs")
                struct_df = pd.DataFrame([{
                    'Type': sm['type'],
                    'Loop Sequence': sm['loop_sequence'],
                    'Count': sm['count'],
                    'Mean Stem Length': f"{sm['mean_stem_length']:.1f} bp",
                    'Sequences': ', '.join(sm['sequences'])
                } for sm in results['structural_motifs']])
                
                st.dataframe(struct_df, use_container_width=True)
            
            # Individual structures - show all in expandable sections
            st.subheader("Individual Structures")
            
            for seq_dict in results['structures']:
                if 'structure' in seq_dict:
                    with st.expander(f"📋 {seq_dict['id']} (MFE: {seq_dict.get('mfe', 0):.2f} kcal/mol)"):
                        col_struct1, col_struct2 = st.columns([1, 1])
                        
                        with col_struct1:
                            st.text(f"Sequence:  {seq_dict.get('folding_sequence', '')}")
                            st.text(f"Structure: {seq_dict['structure']}")
                            st.text(f"MFE: {seq_dict.get('mfe', 0):.2f} kcal/mol")
                        
                        with col_struct2:
                            # Button to generate SVG on demand
                            button_key = f"svg_btn_{seq_dict['id']}"
                            
                            layout = 1 #default
                            
                            if st.button(f"🎨 Generate Structure Diagram", key=button_key):
                                with st.spinner("Generating SVG..."):
                                    try:
                                        from structure_analyzer import StructureAnalyzer
                                        
                                        # Get folding parameters from results
                                        temp = results.get('folding_temperature', 37.0)
                                        salt = results.get('salt_concentration', 1.021)
                                        mg = results.get('magnesium_concentration', 0.0)
                                        
                                        # Convert mM to M if needed
                                        if salt > 10:
                                            salt = salt / 1000.0
                                        mg = mg / 1000.0
                                        
                                        analyzer = StructureAnalyzer(
                                            molecule_type=results.get('molecule_type', 'RNA'),
                                            temperature=temp,
                                            salt_concentration=salt,
                                            magnesium_concentration=mg
                                        )
                                        
                                        svg_content = analyzer.generate_structure_svg(
                                            seq_dict.get('folding_sequence', ''),
                                            seq_dict['structure'],
                                            filename=None,
                                            layout = layout
                                        )
                                        
                                        st.download_button(
                                            label="⬇️ Download structure SVG",
                                            data=svg_content,
                                            file_name=f"{seq_dict['id']}_structure.svg",
                                            mime="image/svg+xml",
                                        )
                                        
                                        if svg_content and '<svg' in svg_content:
                                            # Display folding conditions
                                            st.caption(f"Folding: {temp}°C, {salt*1000:.0f} mM Na+, {mg*1000:.1f} mM Mg2+")
                                            # Display SVG
                                            st.components.v1.html(svg_content, height=400, scrolling=True)
                                        else:
                                            st.error("Could not generate structure diagram. Make sure RNAplot is installed.")
                                            st.code("Install with: sudo apt-get install vienna-rna")
                                    except Exception as e:
                                        st.error(f"Error generating diagram: {str(e)}")
                                        st.info("Make sure ViennaRNA is installed: sudo apt-get install vienna-rna")

            

with tab3:
    st.header("ℹ️ Help & Documentation")
    
    st.markdown("""
    ## About AptaMotif
    
    AptaMotif is a comprehensive tool for discovering and analyzing enriched sequence motifs 
    in aptamer libraries from SELEX experiments.
    
    ### Features
    
    - **Motif Discovery**: Identifies k-mers (5-15 bp by default) that appear in multiple sequences
    - **Statistical Analysis**: Calculates p-values using binomial test and applies FDR correction
    - **Redundancy Filtering**: Automatically removes overlapping/redundant motifs
    - **Merged Motif Detection**: Identifies when short motifs may be fragments of longer patterns
    - **Sequence Identity Check**: Detects duplicate or near-duplicate sequences
    - **Secondary Structure**: Predicts RNA/DNA structures using ViennaRNA with SVG diagrams
    - **Pairwise Similarity**: Shows shared motifs between sequence pairs
    - **Visualization**: Interactive heatmaps, sequence logos, and enrichment plots
    - **Flexible Configuration**: Customizable primers, parameters, and analysis options
    
    ### Workflow
    
    1. **Configure**: Set your primer sequences and analysis parameters in the sidebar
    2. **Input**: Paste or upload your sequenced aptamer clones
    3. **Analyze**: Click "Run Analysis" to process your sequences
    4. **Explore**: View results, visualizations, and download data
    
    ### Statistical Methods
    
    **Binomial Test**: Tests whether each motif appears more frequently than expected by random chance,
    assuming equal base probabilities (25% each for A, C, G, T).
    
    **Fold Enrichment Calculation**:
    ```
    P(motif at one position) = (0.25)^k  where k = motif length
    Number of positions in N-region = L - k + 1  where L = region length
    P(motif in sequence) = 1 - (1 - P_position)^(number of positions)
    Expected sequences with motif = P(motif in sequence) × Total sequences
    
    Fold Enrichment = Observed sequences / Expected sequences
    ```
    
    Example: For motif "GGATCC" (6 bp) in 71 bp N-region with 20 total sequences:
    - P(at position) = 0.25^6 = 0.000244
    - Positions = 71 - 6 + 1 = 66
    - P(in sequence) = 1 - (1 - 0.000244)^66 = 0.0160
    - Expected = 0.0160 × 20 = 0.32 sequences
    - If observed in 10 sequences: Enrichment = 10 / 0.32 = 31.25x
    
    **FDR Correction (Benjamini-Hochberg)**: Controls for multiple testing when examining thousands of motifs.
    Adjusted p-value < 0.05 indicates significant enrichment after correction.
    
    ### Redundancy Filtering
    
    **Problem**: Short overlapping motifs (e.g., "CCTAT" and "TATGG") may both be detected when they're
    actually part of a longer motif ("CCTATGG").
    
    **Solution**: The tool:
    1. Identifies longer motifs that contain shorter ones
    2. Checks if they appear in >80% the same sequences
    3. Keeps only the longer, more informative motif
    4. Attempts to merge overlapping motifs when patterns suggest a longer consensus
    
    ### Sequence Identity Detection
    
    The tool automatically checks for:
    - **Identical sequences**: Exact duplicates (common in cloning artifacts)
    - **Near-identical sequences**: Up to 2 mismatches (PCR errors or sequencing errors)
    
    These are flagged so you can decide whether to remove duplicates before analysis.
    
    ### Structure Prediction
    
    - Uses ViennaRNA RNAfold for minimum free energy (MFE) structure prediction
    - **IMPORTANT**: Folds the FULL sequence (~112 nt), not just the random region
    - Generates SVG diagrams for visual inspection
    - Identifies common structural motifs (hairpin loops, stems)
    - Provides consensus structural statistics
    
    **Folding Parameters** (configurable in sidebar):
    - **Temperature**: Affects thermodynamic calculations (default: 37°C)
      - Lower temp → more stable structures (more negative MFE)
      - Higher temp → less stable structures
    - **Monovalent Salt [Na+]**: Affects electrostatic screening (default: 1.021 M)
      - Higher salt → stabilizes structures (screens phosphate repulsion)
      - Physiological: ~150 mM (0.15 M)
      - SELEX buffers: varies, typically 50-500 mM
    - **Divalent Salt [Mg2+]**: Strong stabilizer (default: 0 mM)
      - Mg2+ is much more effective than Na+ at stabilizing RNA
      - Typical: 0-10 mM in biological conditions
      - SELEX: often 2-5 mM
    
    **Common Conditions**:
    - Physiological: 37°C, 150 mM Na+, 2 mM Mg2+
    - In vitro transcription: 37°C, 1 M Na+, 0 mM Mg2+
    - SELEX selection: Variable (match your buffer!)
    - Room temperature assays: 25°C, buffer-dependent salts
    
    ### Input Format
    
    **FASTA Format**:
    ```
    >Seq1
    GGGAGATACCAGCTTATTCAATT...
    >Seq2
    GGGAGATACCAGCTTATTCAATT...
    ```
    
    **Plain Text** (one sequence per line):
    ```
    GGGAGATACCAGCTTATTCAATT...
    GGGAGATACCAGCTTATTCAATT...
    ```
    
    ### Interpreting Results
    
    **Motif Table Columns**:
    - **Sequence IDs**: NEW - Shows which sequences contain each motif
    - **Adj. p-value < 0.05**: Significantly enriched (after FDR correction)
    - **Fold Enrichment > 2x**: Appears at least twice as often as expected
    - Higher enrichment = stronger selection
    
    **Pairwise Similarity Heatmap**: 
    - Shows number of significant motifs shared between each pair of sequences
    - High values (bright colors) = sequences share many motifs
    - Useful for identifying sequence families or clusters
    
    **Structure Diagrams**:
    - All structures shown in expandable dropdowns (not truncated)
    - Click to view SVG diagram and dot-bracket notation
    - MFE values indicate thermodynamic stability (more negative = more stable)
    
    ### Tips
    
    - Use quality-checked sequences with clean chromatograms
    - Ensure sequences contain both forward and reverse primer regions
    - Check for duplicates in the warnings section
    - Review merged motifs to understand longer patterns
    - Use filters to focus on significant and highly enriched motifs
    - For large datasets (>100 sequences), structure analysis may take several minutes
    - Significant motifs (adjusted p-value < 0.05) are most reliable
    
    ### Citation
    
    If you use AptaMotif in your research, please cite relevant tools:
    - ViennaRNA Package for structure prediction
    - Biopython for sequence analysis
    
    ### Support
    
    For questions or issues, please contact your bioinformatics core facility or Will Euston.
    """)


# Footer
st.markdown("---")
st.markdown("""
<div style='text-align: center; color: #666; padding: 1rem;'>
    <p>AptaMotif v1.2 | Built with Streamlit, Biopython, and ViennaRNA</p>
    <p>🧬 For research use only</p>
</div>
""", unsafe_allow_html=True)
