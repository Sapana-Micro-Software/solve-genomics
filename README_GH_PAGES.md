# GitHub Pages Setup

This repository uses Jekyll with TypeScript, TailwindCSS, and GitHub Actions for automated deployment.

## Setup

1. **Install dependencies:**
   ```bash
   npm install
   bundle install
   ```

2. **Build the site:**
   ```bash
   npm run build
   ```

3. **Development server:**
   ```bash
   npm run dev
   ```

## Structure

- `_config.yml` - Jekyll configuration
- `_layouts/` - Jekyll layouts
- `_typescript/` - TypeScript source files
- `assets/css/input.css` - TailwindCSS input file
- `assets/css/main.css` - Compiled CSS (generated)
- `assets/js/` - Compiled JavaScript (generated)
- `.github/workflows/` - GitHub Actions workflows

## GitHub Actions

The repository includes two workflows:

1. **CI** (`ci.yml`) - Runs on push/PR to test builds
2. **Pages** (`pages.yml`) - Builds and deploys to GitHub Pages

## Technologies

- **Jekyll** - Static site generator
- **TypeScript** - Type-safe JavaScript
- **TailwindCSS** - Utility-first CSS framework
- **Chart.js/D3** - Data visualization (optional)
- **GitHub Actions** - CI/CD automation

## Deployment

The site is automatically deployed to GitHub Pages when changes are pushed to the main branch.

