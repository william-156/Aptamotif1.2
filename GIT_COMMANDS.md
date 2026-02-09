# 🚀 Quick Git Deployment Commands

## First Time Setup

```bash
# Navigate to this folder
cd /path/to/aptamotif_production

# Initialize git
git init

# Add all files
git add .

# Make first commit
git commit -m "Initial commit: AptaMotif with forgi support"

# Add your GitHub repo (create it first on GitHub!)
git remote add origin https://github.com/YOUR_USERNAME/aptamotif.git

# Push to GitHub
git branch -M main
git push -u origin main
```

---

## Updating After Changes

```bash
# Check what changed
git status

# Add changed files
git add .

# Commit with message
git commit -m "Update: describe your changes here"

# Push to GitHub (triggers Streamlit rebuild)
git push
```

---

## Common Git Commands

```bash
# See all files to be committed
git status

# See recent commits
git log --oneline

# Undo changes to a file
git checkout -- filename.py

# Create new branch
git checkout -b feature-name

# Switch branches
git checkout main
```

---

## Troubleshooting

### "fatal: remote origin already exists"
```bash
git remote remove origin
git remote add origin https://github.com/YOUR_USERNAME/aptamotif.git
```

### "fatal: not a git repository"
```bash
cd /path/to/correct/folder
git init
```

### Push rejected (conflicts)
```bash
git pull origin main --rebase
# Fix conflicts if any
git push
```

---

## Before First Push - Checklist

- [ ] Created GitHub repo at github.com
- [ ] Have repo URL ready
- [ ] All files added (`git status` shows green)
- [ ] Committed with good message
- [ ] Remote added correctly

Then: `git push -u origin main`

---

## Your GitHub Repo URL

Replace YOUR_USERNAME with your GitHub username:

```
https://github.com/YOUR_USERNAME/aptamotif.git
```

Don't have a repo yet? Create one:
1. Go to https://github.com/new
2. Name it "aptamotif"
3. Don't initialize with README (we have one)
4. Create repository
5. Copy the URL it gives you
6. Use in `git remote add origin` command

---

## After Pushing to GitHub

1. **Go to Streamlit Cloud**: https://share.streamlit.io/
2. **Click "New app"**
3. **Select** your aptamotif repo
4. **Deploy!**

That's it! 🎉
