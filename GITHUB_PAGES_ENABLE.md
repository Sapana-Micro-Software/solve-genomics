# Enable GitHub Pages

## Automatic Setup via GitHub Actions

Your repository is already configured with GitHub Actions for automatic deployment. Follow these steps:

### Step 1: Enable GitHub Pages

1. Go to your repository: https://github.com/Sapana-Micro-Software/solve-genomics
2. Click on **Settings** (top menu)
3. Scroll down to **Pages** (left sidebar)
4. Under **Source**, select **"GitHub Actions"**
5. Click **Save**

### Step 2: Verify Deployment

1. Go to the **Actions** tab in your repository
2. You should see a workflow run: "Build and Deploy GitHub Pages"
3. Wait for it to complete (usually 2-5 minutes)
4. Once complete, your site will be available at:
   **https://sapana-micro-software.github.io/solve-genomics/**

### Step 3: Check Build Status

If the workflow fails:
- Click on the failed workflow run
- Check the logs for errors
- Common issues:
  - Missing dependencies (should be handled automatically)
  - TypeScript compilation errors
  - Jekyll build errors

## Manual Verification

To test locally before pushing:

```bash
# Install dependencies
npm install
bundle install

# Build
npm run build

# Serve locally
npm run serve
# Or
bundle exec jekyll serve
```

## Repository Structure

The GitHub Pages setup includes:

- ✅ Jekyll configuration (`_config.yml`)
- ✅ TypeScript source (`_typescript/`)
- ✅ TailwindCSS configuration
- ✅ GitHub Actions workflow (`.github/workflows/pages.yml`)
- ✅ Homepage (`index.md`)
- ✅ Layouts (`_layouts/default.html`)
- ✅ Compiled assets (`assets/css/main.css`, `assets/js/`)

## Custom Domain (Optional)

To use a custom domain:

1. Add a `CNAME` file in the root with your domain
2. Update DNS settings as per GitHub's instructions
3. Enable custom domain in Pages settings

## Troubleshooting

### Build Fails
- Check Actions logs for specific errors
- Verify all dependencies are in `package.json` and `Gemfile`
- Ensure TypeScript compiles without errors

### Site Not Updating
- Clear browser cache
- Check Actions tab for latest deployment
- Verify the workflow completed successfully

### 404 Errors
- Check `_config.yml` baseurl setting
- Verify file paths are correct
- Ensure `index.md` or `index.html` exists

## Support

For GitHub Pages documentation:
- https://docs.github.com/en/pages

For Jekyll documentation:
- https://jekyllrb.com/docs/

