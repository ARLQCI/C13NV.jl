# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

C13NV.jl is a Julia package for modeling nitrogen-vacancy (NV) centers in diamonds, focused on quantum control applications. It implements multi-level quantum systems with ground (G), excited (E), and metastable singlet (M) states, including nuclear spin coupling, magnetic field effects, and integration with the QuantumControl.jl ecosystem.

## Development Commands

### Core Development

- `make test` - Run full test suite
- `make coverage` - Run tests with coverage reporting
- `make devrepl` - Start interactive REPL with test dependencies loaded
- `make codestyle` - Format code using JuliaFormatter
- `make docs` - Build documentation
- `make distclean` - Restore the repo to a clean checkout state. This includes removing `Manifest.toml` files.

### Testing

- Tests are organized per-file using TestItemRunner.jl framework
- The entire test suite can be run using `make test`. This mode of test-running does not rely on the environment in the `test` subfolder to be initialized
- In order to run individual tests, the environment in the `test` subfolder must be initialized first, with `make test/Manifest.toml`. If the `test/Manifest.toml` file does not exist, always run `make test/Manifest.toml`.
- Assuming `test/Manifest.toml` exists, a quick way to run the tests in an individual file, e.g., `test/test_quantum_numbers.jl` is via `julia --project=test -e 'include("test/test_quantum_numbers.jl")'`

### Dependencies

- Requires Julia 1.11+
- Main dependencies: QuantumControl.jl, QuantumPropagators.jl

### Physics Conventions

- Frequencies in MHz/GHz, time in μs/ns/ms, magnetic fields in Gauss
- Complex matrices use ComplexF64 precision
- Quantum states constructed with `ket()` and `bra()` functions following standard bra-ket notation

### Integration Points

- Built for QuantumControl.jl optimization workflows
- ComponentArrays.jl for structured parameter handling in optimization problems
- QuantumPropagators.jl for time evolution computations

## Development Environment Setup

The project uses separate environments for different purposes:

- **Test environment** (`test/Project.toml`): Additional testing dependencies
- **Documentation** (`docs/Project.toml`): Documentation building tools
- **Examples** (`examples/Project.toml`): Jupyter notebook examples, saved as `.jl` format via the Jupytext plugin

The `test` environment encompasses the `docs` environment. Run `make test/Manifest.toml` once to make sure the test environment is properly instantiated. Then, `julia --project=test -e …` can run Julia code. For example, to apply code formatting:

	julia --project=test -e 'using JuliaFormatter; format(["src", "docs", "test"])'

This code formatting MUST be run every time any `.jl` file is modified.

## Hamiltonian Notes

The `notes` subdirectory contains a detailed derivation of the system Hamiltonian, as implemented by the `make_nv_system` function. The file `notes/hamiltonian.qmd` is the canonical source, written in Quarto markdown format. The same material is also part of the package documentation, in `docs/src/hamiltonian.md` (in Documenter-markdown format). The file `docs/src/hamiltonian.md` is maintained manually and should not be modified. To help with keeping `docs/src/hamiltonian.md` in sync with `notes/hamiltonian.qmd`, running `make -C notes hamiltonian.md` produces a file `notes/hamiltonian.md` via `pandoc` and the lua-filter `notes/qmd2documenter.lua` that partially translates Quarto markdown to Documenter markdown. The file `notes/hamiltonian.md` should be as close to `docs/src/hamiltonian.md` as possible.

## General Guidelines

* Make sure to only use explicit imports in Julia code, and that there are no imported functions or constants that are not actually used.

* When adding a new dependency to any `Project.toml` file, run `make distclean`, and then `make test/Manifest.toml`, `make docs/Manifest.toml`, etc. to recreate manifest files as necessary.

* Never commit any changes or ask to commit. I will always create git commits manually.
