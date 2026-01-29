.PHONY: help ipynb execute-ipynb codestyle jupyter-lab distclean clean
.DEFAULT_GOAL := help

define PRINT_HELP_JLSCRIPT
rx = r"^([a-z0-9A-Z_-]+):.*?##[ ]+(.*)$$"
for line in eachline()
    m = match(rx, line)
    if !isnothing(m)
        target, help = m.captures
        println("$$(rpad(target, 20)) $$help")
    end
end
endef
export PRINT_HELP_JLSCRIPT

JULIA ?= julia
EXECUTE ?=

JLFILES := $(wildcard *.jl)
NOTEBOOKFILES := $(patsubst %.jl, %.ipynb, $(JLFILES))


help:  ## show this help
	@julia -e "$$PRINT_HELP_JLSCRIPT" < $(MAKEFILE_LIST)

ipynb: $(NOTEBOOKFILES) ## Generate all `.ipynb` files (use `make EXECUTE=--execute ipynb` to execute the examples during the conversion)

execute-ipynb:  ## Equivalent to `make EXECUTE=--execute ipynb`
	make EXECUTE=--execute ipynb

codestyle: Manifest.toml ## Apply the codestyle to all `.jl` files
	$(JULIA) --project=test -e 'using JuliaFormatter; format(["."], verbose=true)'

jupyter-lab: Manifest.toml ## Run a Jupyter lab server
	JULIA_NUM_THREADS=auto jupyter lab --no-browser

distclean:  # Remove all generated files
	rm -f $(NOTEBOOKFILES) Manifest.toml
	rm -f .jupyter_ystore.db

clean:  # Clean up temporary files

%.ipynb : | %.md Manifest.toml
	NUMEXPR_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 time jupytext --to notebook $(EXECUTE) $(*).md
	jupyter trust "$(*).ipynb"

%.ipynb : | %.py
	NUMEXPR_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 time jupytext --to notebook $(EXECUTE) $(*).py
	jupyter trust "$(*).ipynb"

%.ipynb : | %.jl Manifest.toml
	JULIA_NUM_THREADS=auto NUMEXPR_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 time jupytext --to notebook $(EXECUTE) $(*).jl
	jupyter trust "$(*).ipynb"

%.ipynb : | %.qmd Manifest.toml
	JULIA_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 time jupytext --to notebook $(EXECUTE) $(*).qmd
	jupyter trust "$(*).ipynb"

Manifest.toml: Project.toml
	$(JULIA) --project=. -e 'using Pkg; Pkg.instantiate()'
