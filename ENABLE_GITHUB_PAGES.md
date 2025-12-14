# Enable GitHub Pages - Quick Guide

## Current Status

✅ GitHub Actions workflow configured (`.github/workflows/pages.yml`)
✅ Jekyll site configured (`_config.yml`)
✅ TypeScript and TailwindCSS build scripts ready
✅ All files pushed to repository

## Enable GitHub Pages

### Option 1: Via GitHub Web Interface (Recommended)

1. **Navigate to repository settings:**
   - Go to: https://github.com/Sapana-Micro-Software/solve-genomics/settings/pages

2. **Configure Pages:**
   - Under "Source", select **"GitHub Actions"**
   - Click "Save"

3. **Verify:**
   - Go to the "Actions" tab
   - You should see the "Build and Deploy GitHub Pages" workflow running
   - Once complete, your site will be available at:
     `https://sapana-micro-software.github.io/solve-genomics/`

### Option 2: Via GitHub CLI

```bash
gh api repos/Sapana-Micro-Software/solve-genomics/pages \
  -X POST \
  -f source[type]=actions
```

### Option 3: Manual Trigger

If Pages is already enabled but the workflow hasn't run:

1. Go to: https://github.com/Sapana-Micro-Software/solve-genomics/actions
2. Select "Build and Deploy GitHub Pages"
3. Click "Run workflow" → "Run workflow"

## What Happens Next

1. **GitHub Actions will:**
   - Checkout the code
   - Setup Ruby and Node.js
   - Install dependencies (npm, bundler)
   - Build TypeScript → JavaScript
   - Build TailwindCSS → CSS
   - Build Jekyll site → `_site/`
   - Deploy to GitHub Pages

2. **Site will be available at:**
   - `https://sapana-micro-software.github.io/solve-genomics/`

3. **Future updates:**
   - Every push to `main` branch will automatically rebuild and deploy

## Troubleshooting

- **Workflow not running?** Check repository Settings → Actions → General → Workflow permissions
- **Build failing?** Check the Actions tab for error logs
- **Site not updating?** Wait a few minutes for deployment to complete

## Repository Links

- **Repository:** https://github.com/Sapana-Micro-Software/solve-genomics
- **Settings:** https://github.com/Sapana-Micro-Software/solve-genomics/settings
- **Pages Settings:** https://github.com/Sapana-Micro-Software/solve-genomics/settings/pages
- **Actions:** https://github.com/Sapana-Micro-Software/solve-genomics/actions

