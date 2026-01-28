.PHONY: help all docx pdf preview preview-pdf preview-all clean distclean
.DEFAULT_GOAL := pdf

NAME = hamiltonian

ALL = $(NAME).docx $(NAME).pdf $(NAME).html

help:   ## Show this help
	@grep -E '^([a-zA-Z_-]+):.*## ' $(MAKEFILE_LIST) | awk -F ':.*## ' '{printf "%-20s %s\n", $$1, $$2}'

%.docx: %.qmd  $(wildcard images/*.tikz images/*.tex images/*.pdf images/*.png)
	make -C images all
	quarto render $< --to docx

%.pdf: %.qmd $(wildcard images/*.tikz images/*.tex images/*.pdf images/*.png)
	make -C images all
	quarto render $< --to pdf

%.html: %.qmd $(wildcard images/*.tikz images/*.tex images/*.pdf images/*.png)
	make -C images all
	quarto render $< --to html

all:  $(ALL)  ## Generate all output files

docx: $(NAME).docx  ## Generate the docx output file

pdf: $(NAME).pdf  ## Generate the pdf output file

preview:  ## Open an html preview
	quarto preview $(NAME).qmd --to html

preview-pdf:  ## Open a preview of the pdf version
	quarto preview $(NAME).qmd --to pdf

preview-all:  ## Open a preview of the all output formas
	quarto preview $(NAME).qmd

clean:  ## Delete compilation artifact files
	make -C images clean
	rm -f $(NAME).tex
	rm -f $(NAME).html
	rm -rf $(NAME)_files

distclean: clean  ## Delete output files
	make -C images distclean
	rm -f $(ALL)
