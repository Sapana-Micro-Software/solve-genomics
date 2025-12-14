# GitHub Pages Build Fix

## Issues Fixed

### 1. Bundle Install Failure (Exit Code 4)
**Problem:** The `bundle install` command was failing with exit code 4, typically caused by:
- Jekyll version incompatibility with GitHub Pages
- Missing or incorrect Gemfile configuration
- Bundler cache conflicts

**Solution:**
- Switched from individual Jekyll gems to the `github-pages` gem, which ensures compatibility
- Updated workflow to:
  - Disable bundler-cache initially to avoid conflicts
  - Explicitly install bundler v2.0
  - Set bundle config to use local path (`vendor/bundle`)
  - Add retry logic (`--retry 3`) for network issues

### 2. Gemfile Updates
**Before:**
```ruby
gem "jekyll", "~> 4.3"
gem "jekyll-feed", "~> 0.15"
# ... other individual gems
```

**After:**
```ruby
gem "github-pages", group: :jekyll_plugins
gem "jekyll-feed", "~> 0.15"
# ... additional plugins
```

The `github-pages` gem includes:
- Compatible Jekyll version
- All necessary plugins
- GitHub Pages-specific configurations

### 3. Workflow Improvements
- Added explicit bundler installation
- Configured bundle to use local path (avoids permission issues)
- Added retry logic for network failures
- Proper environment variables for Jekyll build

### 4. Configuration Updates
- Updated `_config.yml` baseurl to `/solve-genomics`
- Removed theme override (handled by github-pages gem)

## Files Changed

1. **Gemfile**: Switched to `github-pages` gem
2. **.github/workflows/pages.yml**: Improved bundle install process
3. **_config.yml**: Updated baseurl and removed theme override

## Testing

The workflow should now:
1. ✅ Successfully install Ruby dependencies
2. ✅ Build TypeScript and TailwindCSS
3. ✅ Build Jekyll site without errors
4. ✅ Deploy to GitHub Pages

## Next Steps

1. Enable GitHub Pages in repository settings:
   - Go to: https://github.com/Sapana-Micro-Software/solve-genomics/settings/pages
   - Select "GitHub Actions" as source
   - Save

2. Monitor the workflow:
   - Check Actions tab for build status
   - Site will be available at: https://sapana-micro-software.github.io/solve-genomics/

## Troubleshooting

If build still fails:
1. Check Actions logs for specific error messages
2. Verify Ruby version compatibility (currently using 3.1)
3. Ensure all dependencies are listed in Gemfile
4. Check that `_config.yml` is valid YAML

