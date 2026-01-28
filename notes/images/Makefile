.PHONY: help all clean distclean

PDFLATEXOPTS = -file-line-error -interaction=nonstopmode -halt-on-error
ALL = levels.pdf liouville_levels.pdf

help:   ## Show this help
	@grep -E '^([a-zA-Z_-]+):.*## ' $(MAKEFILE_LIST) | awk -F ':.*## ' '{printf "%-20s %s\n", $$1, $$2}'

all: $(ALL)  ## Compile all files

%.pdf: %.tex %.tikz
	@echo "\n** generate figure: $< -> $@"
	@lualatex $(PDFLATEXOPTS) $<

clean:
	@rm -f *.runlog
	@rm -f *.aux
	@rm -f *.log
	@rm -f *.toc
	@rm -f *.blg
	@rm -rf *.out
	@rm -f *.bak
	@rm -f *.ilg
	@rm -f *.snm
	@rm -f *.nav
	@rm -f *.table
	@rm -f *.dvi
	@rm -f *.fls
	@rm -f *.bbl
	@rm -f *~
	@rm -f *Notes.bib
	@rm -f *-eps-converted-to.pdf
	@rm -f *.fdb_latexmk
	@rm -f *.synctex.gz*
	@rm -f ./latexrun
	@rm -rf __pycache__
	@rm -f results/*.pdf

distclean: clean
	rm -f $(ALL)
