# Makefile for building LaTeX documents

.PHONY: all paper presentation clean view-paper view-presentation

# LaTeX compiler
LATEX = pdflatex
BIBTEX = bibtex

# Output directory
OUTPUT_DIR = build

# Source files
PAPER_SRC = paper.tex
PRESENTATION_SRC = presentation.tex

# Output files
PAPER_PDF = $(OUTPUT_DIR)/paper.pdf
PRESENTATION_PDF = $(OUTPUT_DIR)/presentation.pdf

all: paper presentation

paper: $(PAPER_PDF)

presentation: $(PRESENTATION_PDF)

# Create output directory
$(OUTPUT_DIR):
	mkdir -p $(OUTPUT_DIR)

# Build paper
$(PAPER_PDF): $(PAPER_SRC) | $(OUTPUT_DIR)
	$(LATEX) -output-directory=$(OUTPUT_DIR) $(PAPER_SRC)
	$(LATEX) -output-directory=$(OUTPUT_DIR) $(PAPER_SRC)
	@echo "Paper built: $(PAPER_PDF)"

# Build presentation
$(PRESENTATION_PDF): $(PRESENTATION_SRC) | $(OUTPUT_DIR)
	$(LATEX) -output-directory=$(OUTPUT_DIR) $(PRESENTATION_SRC)
	$(LATEX) -output-directory=$(OUTPUT_DIR) $(PRESENTATION_SRC)
	@echo "Presentation built: $(PRESENTATION_PDF)"

# View documents (platform-specific)
view-paper: $(PAPER_PDF)
	@if command -v open > /dev/null; then \
		open $(PAPER_PDF); \
	elif command -v xdg-open > /dev/null; then \
		xdg-open $(PAPER_PDF); \
	elif command -v evince > /dev/null; then \
		evince $(PAPER_PDF) &; \
	else \
		echo "Please open $(PAPER_PDF) manually"; \
	fi

view-presentation: $(PRESENTATION_PDF)
	@if command -v open > /dev/null; then \
		open $(PRESENTATION_PDF); \
	elif command -v xdg-open > /dev/null; then \
		xdg-open $(PRESENTATION_PDF); \
	elif command -v evince > /dev/null; then \
		evince $(PRESENTATION_PDF) &; \
	else \
		echo "Please open $(PRESENTATION_PDF) manually"; \
	fi

# Clean build artifacts
clean:
	rm -rf $(OUTPUT_DIR)
	rm -f *.aux *.log *.out *.toc *.nav *.snm *.vrb *.bbl *.blg

# Help
help:
	@echo "Available targets:"
	@echo "  all              - Build both paper and presentation"
	@echo "  paper            - Build paper only"
	@echo "  presentation     - Build presentation only"
	@echo "  view-paper       - Build and view paper"
	@echo "  view-presentation - Build and view presentation"
	@echo "  clean            - Remove build artifacts"
	@echo "  help             - Show this help message"

