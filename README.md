# 🧬 AptaMotif - Aptamer Motif Discovery & Analysis

A comprehensive web-based tool for discovering enriched sequence motifs and analyzing secondary structures in aptamer libraries from SELEX experiments.

**Live Demo**: https://aptamotif.streamlit.app/

---

## ✨ Features

- **Motif Discovery**: Identifies enriched k-mers (5-25 bp) shared across sequences
- **Statistical Analysis**: Binomial test with FDR (Benjamini-Hochberg) correction
- **Smart Filtering**: Removes redundant motifs, identifies merged patterns
- **Duplicate Detection**: Finds identical and near-identical sequences
- **Secondary Structure**: RNA/DNA structure prediction using ViennaRNA
- **Folding Parameters**: Customizable temperature, Na+, and Mg2+ concentrations
- **Multi-Tier Visualization**: 
  - **RNAplot** (best quality, HuggingFace/local)
  - **forgi** (good quality, works on Streamlit Cloud!)
  - Text fallback (always works)
- **Interactive Plots**: Heatmaps, enrichment plots, volcano plots, sequence logos
- **Export Everything**: CSV tables, SVG structures, publication figures

---

## 🚀 Quick Start

### Using the Live App

1. Go to https://aptamotif.streamlit.app/
2. Paste your sequences (FASTA format)
3. Configure settings in sidebar
4. Click "Run Analysis"
5. Download results!

### Example Input

```
>Clone1
GGGAGATACCAGCTTATTCAATTGGATCCGGATCCAGATAGTAAGTGCAATCT
>Clone2
GGGAGATACCAGCTTATTCAATTCCAATTCCAATTCCAGATAGTAAGTGCAATCT
```

---

## 📦 Installation (Local/HuggingFace)

### Option 1: Streamlit Cloud (Current Deployment)

**What works:**
- ✅ Motif discovery
- ✅ Statistical analysis
- ✅ Structure prediction (ViennaRNA Python)
- ✅ **forgi visualization** (pure Python)
- ❌ RNAplot (requires system packages)

**Deployed at**: https://aptamotif.streamlit.app/

### Option 2: HuggingFace Spaces (Full Features)

**What works:**
- ✅ Everything from Streamlit Cloud
- ✅ **RNAplot visualization** (best quality!)
- ✅ Full ViennaRNA support

**To deploy on HuggingFace:**

1. Create new Space at https://huggingface.co/spaces
2. Choose Streamlit SDK
3. Upload all files from this folder
4. Add `packages.txt` with:
   ```
   vienna-rna
   ```
5. Wait for build (~5-10 min)
6. Done!

### Option 3: Run Locally

```bash
# Clone repo
git clone https://github.com/YOUR_USERNAME/aptamotif.git
cd aptamotif

# Install dependencies
pip install -r requirements.txt

# Optional: Install ViennaRNA for RNAplot
# macOS:
brew install viennarna
# Ubuntu:
sudo apt-get install vienna-rna
# Conda:
conda install -c bioconda viennarna

# Run app
streamlit run app.py
```

---

## 🎨 Visualization Methods

### Tier 1: RNAplot (Best Quality ⭐⭐⭐⭐⭐)
- Spring-embedded layout
- Publication quality
- **Requires**: System ViennaRNA install
- **Works on**: HuggingFace Spaces, local machines
- **Not on**: Streamlit Cloud

### Tier 2: forgi (Good Quality ⭐⭐⭐⭐)
- Graph-based layout  
- Professional appearance
- **Requires**: Only Python packages (forgi, matplotlib)
- **Works on**: Everywhere including Streamlit Cloud!
- **NEW**: Just added in this version

### Tier 3: Text Fallback (Always Works ⭐⭐)
- Shows dot-bracket notation
- All structure data visible
- No graphical diagram

**The app automatically uses the best available method!**

---

## 📋 Project Structure

```
aptamotif_production/
├── app.py                      # Main Streamlit application
├── sequence_parser.py          # Sequence parsing & primer extraction
├── motif_finder.py            # Motif discovery engine
├── motif_statistics.py        # Statistical analysis
├── structure_analyzer.py      # Structure prediction + visualization (NEW!)
├── visualizer.py              # Heatmaps, plots, logos
├── requirements.txt           # Python dependencies (includes forgi!)
├── example_sequences.fasta    # Test data
└── README.md                  # This file
```

---

## 🔧 What's New in This Version

### ✅ forgi Integration (Major Update!)

**Problem Solved**: SVG structure diagrams now work on Streamlit Cloud!

**Changes**:
1. Added `forgi` library to requirements.txt
2. Updated `structure_analyzer.py` with three-tier approach:
   - Try RNAplot (if available)
   - Fall back to forgi (pure Python)
   - Fall back to text (always works)
3. Automatic method selection

**Result**: Structure visualization now works everywhere!

---

## 🚢 Deployment Instructions

### Deploy to Your Streamlit Cloud

1. **Push to GitHub:**
   ```bash
   git init
   git add .
   git commit -m "Initial commit with forgi support"
   git remote add origin https://github.com/YOUR_USERNAME/aptamotif.git
   git push -u origin main
   ```

2. **Connect to Streamlit:**
   - Go to https://share.streamlit.io/
   - Click "New app"
   - Connect your GitHub repo
   - Select `app.py` as main file
   - Click "Deploy"

3. **Wait ~5 minutes** for build

4. **Test structure visualization** - should work with forgi!

### Deploy to HuggingFace Spaces

1. **Create `packages.txt`:**
   ```
   vienna-rna
   ```

2. **Push to HuggingFace:**
   ```bash
   git clone https://huggingface.co/spaces/YOUR_USERNAME/aptamotif
   cd aptamotif
   # Copy all files here
   git add .
   git commit -m "Deploy with forgi + RNAplot support"
   git push
   ```

3. **Wait ~10 minutes** for build (ViennaRNA takes time)

4. **Test** - should have both RNAplot AND forgi available!

---

## 📊 Requirements

### Python Packages (requirements.txt)
```
streamlit>=1.28.0
biopython>=1.81
pandas>=2.0.0
numpy>=1.24.0
scipy>=1.10.0
matplotlib>=3.7.0
seaborn>=0.12.0
logomaker>=0.8
plotly>=5.14.0
statsmodels>=0.14.0
ViennaRNA>=2.6.0
forgi>=2.0.3          # NEW!
networkx>=2.8.0       # NEW! (required by forgi)
```

### System Packages (packages.txt - HuggingFace only)
```
vienna-rna
```

---

## 🧪 Testing

### Test Locally Before Deploying

```bash
# Run app
streamlit run app.py

# Test with example sequences
# Click "Load Example Sequences"
# Run analysis
# Try "Generate Structure Diagram" button

# Check terminal output:
# Should see: "✓ Using forgi visualization"
```

### Verify Each Method Works

**Test RNAplot (if available):**
```bash
which RNAplot
echo -e "GGGAAACCC\n(((...)))" | RNAplot --output-format=svg
```

**Test forgi:**
```python
import forgi
import forgi.visual.mplotlib as fvm
bg = forgi.graph.bulge_graph.BulgeGraph.from_dotbracket("(((...)))")
fvm.plot_rna(bg)
```

---

## 🎯 Usage Tips

### For Your Lab

1. **Streamlit Cloud** (current): Good for daily use, forgi visualization
2. **Local install** with ViennaRNA: Best quality for publications
3. **HuggingFace**: Best of both worlds (if you want to migrate)

### For Publications

- Use local version with RNAplot for best quality figures
- Or use HuggingFace Space (has RNAplot)
- forgi quality is still publication-acceptable!

---

## 🐛 Troubleshooting

### Structure Visualization Not Working

**Check which method is being used:**
- Look at terminal/logs
- Should see: `✓ Using forgi visualization`

**If you see errors:**
```bash
# Test forgi installation
python -c "import forgi; print('forgi works!')"

# Reinstall if needed
pip install --upgrade forgi networkx
```

### Streamlit Cloud Build Fails

**Common issue**: Missing dependencies

**Solution**: Verify requirements.txt includes:
- `forgi>=2.0.3`
- `networkx>=2.8.0`

---

## 📖 Documentation

### Structure Prediction

- **Method**: ViennaRNA `RNA.fold_compound()`
- **Parameters**: Temperature, Na+, Mg2+ (customizable)
- **Output**: Dot-bracket notation, MFE, structural elements

### Structure Visualization

- **Method 1 (Best)**: RNAplot - requires system install
- **Method 2 (Good)**: forgi - pure Python, works everywhere
- **Method 3 (Basic)**: Text fallback

**The app chooses automatically!**

---

## 🤝 Contributing

Want to improve AptaMotif?

1. Fork the repo
2. Make changes
3. Test locally
4. Submit pull request

---

## 📝 License

MIT License - Free for academic and commercial use

---

## 🙏 Acknowledgments

Built with:
- **Streamlit** - Web framework
- **ViennaRNA** - Structure prediction
- **forgi** - Pure Python structure visualization (NEW!)
- **Biopython** - Sequence handling
- **Plotly** - Interactive plots
- **Logomaker** - Sequence logos

---

## 📞 Support

- **Issues**: Report on GitHub Issues
- **Questions**: Check documentation in Help tab
- **Lab use**: Clone and customize for your needs!

---

**Version**: 2.0 (forgi integration)  
**Last Updated**: 2026-02-08  
**Status**: Production Ready ✅

---

## 🚀 Ready to Deploy!

This version includes:
- ✅ All core functionality
- ✅ forgi support for Streamlit Cloud
- ✅ RNAplot support for HuggingFace/local
- ✅ Clean, documented code
- ✅ Example sequences
- ✅ Complete documentation

**Just push to Git and deploy!** 🎉
