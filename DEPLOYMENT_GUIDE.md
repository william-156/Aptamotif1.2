# 🚀 Deployment Guide - Step by Step

## 📦 What You Have

A production-ready AptaMotif package with:
- ✅ forgi support (works on Streamlit Cloud!)
- ✅ RNAplot support (works on HuggingFace/local)
- ✅ Clean, tested code
- ✅ All features working

---

## 🎯 Quick Deploy to Streamlit Cloud

### Step 1: Push to Your GitHub

```bash
# Navigate to this folder
cd /path/to/aptamotif_production

# Initialize git (if not already)
git init

# Add all files
git add .

# Commit
git commit -m "Add forgi support for structure visualization"

# Connect to your GitHub repo
# (Replace with your actual GitHub username and repo)
git remote add origin https://github.com/YOUR_USERNAME/aptamotif.git

# Push
git branch -M main
git push -u origin main
```

### Step 2: Deploy on Streamlit Cloud

1. Go to https://share.streamlit.io/
2. Sign in with GitHub
3. Click "New app"
4. Select:
   - **Repository**: YOUR_USERNAME/aptamotif
   - **Branch**: main
   - **Main file path**: app.py
5. Click "Deploy!"

### Step 3: Wait for Build (~5 minutes)

Streamlit will:
- Install Python packages
- Install forgi
- Build your app

### Step 4: Test!

1. Once deployed, visit your app URL
2. Load example sequences
3. Run analysis
4. Try "Generate Structure Diagram"
5. Should see forgi visualization! ✓

---

## 📝 Files to Commit

Make sure these files are in your repo:

### Core Files (Required)
- ✅ `app.py` - Main application
- ✅ `sequence_parser.py`
- ✅ `motif_finder.py`
- ✅ `motif_statistics.py`
- ✅ `structure_analyzer.py` - **Updated with forgi!**
- ✅ `visualizer.py`
- ✅ `requirements.txt` - **Includes forgi!**

### Documentation & Examples
- ✅ `README.md` - User documentation
- ✅ `example_sequences.fasta` - Test data
- ✅ `DEPLOYMENT_GUIDE.md` - This file

### Optional (Not needed for deployment)
- `.gitignore` - Recommended
- `LICENSE` - If making public

---

## 📋 Pre-Deployment Checklist

Before pushing, verify:

### ✓ requirements.txt includes forgi
```bash
grep "forgi" requirements.txt
# Should show: forgi>=2.0.3
```

### ✓ structure_analyzer.py has new functions
```bash
grep "_try_forgi" structure_analyzer.py
# Should find the function
```

### ✓ All imports work locally
```bash
python -c "from structure_analyzer import StructureAnalyzer; print('OK')"
```

### ✓ Git status is clean
```bash
git status
# Check all files are added
```

---

## 🔄 Updating Your Deployed App

After initial deployment, to update:

```bash
# Make your changes
# Then:

git add .
git commit -m "Description of changes"
git push

# Streamlit Cloud auto-detects and rebuilds!
```

---

## 🌐 Deployment Options Comparison

### Option 1: Streamlit Cloud (Current)

**Pros:**
- ✅ Free tier
- ✅ Easy deployment
- ✅ Auto-updates from Git
- ✅ **forgi works!** (NEW)

**Cons:**
- ❌ No system packages (no RNAplot)
- ⚠️ CPU performance (free tier)

**Best for**: Lab daily use, quick analysis

### Option 2: HuggingFace Spaces

**Pros:**
- ✅ Free tier
- ✅ System package support
- ✅ Both forgi AND RNAplot work!
- ✅ Better documentation platform

**Cons:**
- ⚠️ Slower cold starts
- ⚠️ Slightly more setup

**Best for**: Public sharing, best quality output

### Option 3: Local Installation

**Pros:**
- ✅ Fastest performance
- ✅ Full control
- ✅ All features
- ✅ Private data

**Cons:**
- ❌ Requires installation
- ❌ Not shareable

**Best for**: Publications, sensitive data

---

## 🎨 Customization Before Deploy

### Update App Title/Branding

Edit `app.py`:
```python
st.set_page_config(
    page_title="YourLab AptaMotif",  # Change this
    page_icon="🧬",
    ...
)
```

### Add Your Lab Info

Edit `app.py` header:
```python
st.markdown("**Your Lab Name**")
st.markdown("Contact: email@university.edu")
```

### Change Default Settings

Edit sidebar defaults in `app.py`:
```python
forward_primer = st.sidebar.text_input(
    "Forward Primer",
    value="YOUR_DEFAULT_PRIMER",  # Change
    ...
)
```

---

## 🐛 Troubleshooting Deployment

### Build Fails on Streamlit

**Check logs:**
- Click "Manage app" in Streamlit Cloud
- Check "Logs" tab
- Look for error messages

**Common issues:**
1. **Typo in requirements.txt** → Fix spelling
2. **Import error** → Check all files pushed
3. **Syntax error** → Test locally first

**Solution:**
```bash
# Fix the issue
git add .
git commit -m "Fix build error"
git push
# Streamlit auto-rebuilds
```

### forgi Not Working

**Symptoms**: Falls back to text visualization

**Check:**
```bash
# In Streamlit logs, should see:
"✓ Using forgi visualization"

# If you see:
"forgi not installed"
# → requirements.txt missing forgi
```

**Fix:**
```bash
# Add to requirements.txt:
forgi>=2.0.3
networkx>=2.8.0

# Push update
git add requirements.txt
git commit -m "Add forgi dependency"
git push
```

### App Works But Slow

**Causes:**
- Large sequence files
- Many sequences (>100)
- Complex visualizations

**Solutions:**
1. Optimize settings (reduce motif length range)
2. Disable structure prediction if not needed
3. Use local installation for large jobs
4. Consider upgrading Streamlit tier

---

## ✅ Verification Steps After Deploy

### 1. Basic Functionality
- [ ] App loads without errors
- [ ] Can input sequences
- [ ] Run analysis works
- [ ] Results display correctly

### 2. Motif Discovery
- [ ] Finds motifs in example data
- [ ] Statistical analysis works
- [ ] Visualizations appear
- [ ] Can download CSV

### 3. Structure Prediction
- [ ] Structures fold correctly
- [ ] MFE values reasonable
- [ ] Dot-bracket notation shown

### 4. Structure Visualization (KEY!)
- [ ] "Generate Structure Diagram" button appears
- [ ] Click generates SVG
- [ ] Check logs: "✓ Using forgi visualization"
- [ ] SVG displays in app
- [ ] No fallback text message

---

## 📊 Expected Build Time

### Streamlit Cloud
- First build: ~5-7 minutes
- Subsequent builds: ~2-3 minutes
- Installing forgi: ~30 seconds

### HuggingFace Spaces
- First build: ~10-15 minutes (ViennaRNA takes time)
- Subsequent builds: ~5-8 minutes

---

## 🎉 Success Criteria

You know it's working when:

1. **App loads** at your Streamlit URL
2. **Example sequences work** without errors
3. **Structure visualization shows** forgi SVG (not fallback text)
4. **Logs confirm**: `✓ Using forgi visualization`
5. **Users can download** results

---

## 📞 Getting Help

### If Stuck:

1. **Check Streamlit logs** for errors
2. **Test locally first** (`streamlit run app.py`)
3. **Verify all files pushed** (`git status`)
4. **Check requirements.txt** has all packages
5. **Search Streamlit forums** for similar issues

### Common Resources:

- Streamlit Docs: https://docs.streamlit.io/
- forgi GitHub: https://github.com/ViennaRNA/forgi
- This README: All usage details

---

## 🚀 Ready to Deploy!

**Quick checklist:**
- [ ] All files in this folder
- [ ] Git repo initialized
- [ ] Pushed to GitHub
- [ ] Connected to Streamlit Cloud
- [ ] Deployed and building

**Once deployed:**
- [ ] Test with example sequences
- [ ] Verify forgi visualization works
- [ ] Share URL with your lab!

---

## 🎯 Next Steps After Successful Deployment

1. **Share with lab** → Send them the URL
2. **Get feedback** → What features do they want?
3. **Monitor usage** → Check Streamlit analytics
4. **Iterate** → Add features, fix bugs
5. **Consider HuggingFace** → If want RNAplot too

---

**Good luck with deployment! 🚀**

The app is ready to go. Just push to Git and deploy to Streamlit Cloud!

Questions? Check the logs. Most issues show clear error messages that are easy to fix.
