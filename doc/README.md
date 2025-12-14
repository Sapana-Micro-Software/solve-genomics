# Documentation

This directory contains LaTeX documentation for the DNA Sequence Alignment project.

## Files

- **`paper.tex`**: Comprehensive research paper with complete analysis
- **`presentation.tex`**: Beamer presentation slides
- **`Makefile`**: Build system for LaTeX documents

## Building Documents

### Prerequisites

Install LaTeX distribution:
- **macOS**: `brew install --cask mactex` or `brew install basictex`
- **Linux**: `sudo apt-get install texlive-full` (Ubuntu/Debian) or `sudo yum install texlive-scheme-full` (RHEL/CentOS)
- **Windows**: Install MiKTeX or TeX Live

### Build Commands

```bash
# Build both paper and presentation
make all

# Build paper only
make paper

# Build presentation only
make presentation

# Build and view paper
make view-paper

# Build and view presentation
make view-presentation

# Clean build artifacts
make clean
```

### Manual Build

```bash
# Paper
pdflatex paper.tex
pdflatex paper.tex  # Run twice for references

# Presentation
pdflatex presentation.tex
pdflatex presentation.tex  # Run twice for references
```

## Paper Contents

The paper (`paper.tex`) includes:

1. **Introduction**: Background, motivation, contributions
2. **Related Work**: Key algorithms and references
3. **Methodology**: Algorithm categories and implementations
4. **Experimental Setup**: Test sequences and benchmark configuration
5. **Results and Analysis**: 
   - Performance comparisons
   - Complexity analysis
   - Scalability analysis
   - Accuracy analysis
   - Compression effectiveness
6. **Discussion**: Algorithm selection guidelines, trade-offs, limitations
7. **Conclusion**: Summary and future work
8. **References**: Bibliographic citations

## Presentation Contents

The presentation (`presentation.tex`) includes:

1. **Introduction**: Human DNA overview, problem statement
2. **Algorithms Overview**: All algorithm categories
3. **Exact Matching**: Performance and algorithm details
4. **Dynamic Programming**: Smith-Waterman and Needleman-Wunsch
5. **Approximate Matching**: Fuzzy search and WARP-CTC
6. **Compression Methods**: Grammar-based and lossy compression
7. **Modern Approaches**: Embeddings, deep learning, MCMC
8. **Parallel Methods**: Scaling and distributed search
9. **Benchmark Results**: Comprehensive performance data
10. **Analysis and Discussion**: Guidelines and trade-offs
11. **Conclusion**: Summary and future directions

## Customization

### Adding Figures

1. Place figure files in `figures/` directory
2. Use `\includegraphics{figures/filename}` in LaTeX
3. Supported formats: PDF, PNG, JPG

### Modifying Content

- Edit `paper.tex` for paper content
- Edit `presentation.tex` for presentation content
- Both use standard LaTeX/Beamer syntax

### Changing Style

- Paper: Modify document class options in `\documentclass[...]{article}`
- Presentation: Change Beamer theme in `\usetheme{...}`

## Output

Built PDFs are placed in `build/` directory:
- `build/paper.pdf`
- `build/presentation.pdf`

## Notes

- Run LaTeX twice to resolve cross-references
- Bibliography requires BibTeX if using `.bib` files
- Some LaTeX distributions may require additional packages
- For best results, use a complete LaTeX distribution (TeX Live Full or MacTeX)

