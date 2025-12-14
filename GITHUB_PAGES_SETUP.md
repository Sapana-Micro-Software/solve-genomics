# GitHub Pages Setup Guide

This repository is configured for GitHub Pages using Jekyll, TypeScript, TailwindCSS, and GitHub Actions.

## Quick Start

1. **Enable GitHub Pages:**
   - Go to repository Settings → Pages
   - Source: GitHub Actions
   - The site will auto-deploy on push to main/master

2. **Local Development:**
   ```bash
   # Install dependencies
   npm install
   bundle install
   
   # Build and serve
   npm run dev
   ```

## File Structure

```
├── _config.yml              # Jekyll configuration
├── _layouts/                # Jekyll layouts
│   └── default.html         # Default page layout
├── _typescript/             # TypeScript source files
│   ├── main.ts             # Main application logic
│   └── charts.ts           # Chart utilities
├── assets/
│   ├── css/
│   │   ├── input.css       # TailwindCSS input
│   │   └── main.css        # Compiled CSS (generated)
│   └── js/
│       ├── main.js         # Compiled JS (generated)
│       └── charts.js       # Compiled charts (generated)
├── .github/workflows/
│   ├── pages.yml           # GitHub Pages deployment
│   └── ci.yml              # CI workflow
├── index.md                # Homepage (Markdown)
├── package.json            # Node.js dependencies
├── Gemfile                 # Ruby/Jekyll dependencies
├── tsconfig.json           # TypeScript configuration
├── tailwind.config.js      # TailwindCSS configuration
└── postcss.config.js       # PostCSS configuration
```

## Technologies

- **Jekyll 4.3** - Static site generator
- **TypeScript 5.3** - Type-safe JavaScript
- **TailwindCSS 3.3** - Utility-first CSS
- **GitHub Actions** - CI/CD automation
- **Node.js 20** - JavaScript runtime
- **Ruby 3.1** - Jekyll runtime

## Build Process

1. TypeScript → JavaScript (`_typescript/` → `assets/js/`)
2. TailwindCSS → CSS (`assets/css/input.css` → `assets/css/main.css`)
3. Jekyll → Static HTML (`*.md`, `_layouts/` → `_site/`)

## GitHub Actions Workflows

### Pages Workflow (`pages.yml`)
- Triggers on push to main/master
- Builds TypeScript
- Builds TailwindCSS
- Builds Jekyll site
- Deploys to GitHub Pages

### CI Workflow (`ci.yml`)
- Runs on push/PR
- Tests TypeScript compilation
- Tests CSS build
- Validates builds

## Customization

### Update Site Info
Edit `_config.yml`:
```yaml
title: Your Title
author: Your Name
url: https://yourusername.github.io
baseurl: /repository-name
```

### Add Pages
Create `.md` files in root or `_pages/` directory with front matter:
```yaml
---
layout: default
title: Page Title
---
```

### Modify Styles
Edit `assets/css/input.css` and run `npm run build:css`

### Add TypeScript
Add files to `_typescript/` and import in `main.ts`

## Deployment

The site automatically deploys when you push to the main branch. Check Actions tab for build status.

## Troubleshooting

- **Build fails:** Check GitHub Actions logs
- **Styles not loading:** Run `npm run build:css`
- **JS not working:** Run `npm run build:ts`
- **Jekyll errors:** Check `_config.yml` syntax

