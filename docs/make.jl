using C13NV
using Documenter
using DocumenterCitations
import Pkg


PROJECT_TOML = Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
VERSION = PROJECT_TOML["version"]
NAME = PROJECT_TOML["name"]
AUTHORS = join(PROJECT_TOML["authors"], ", ") * " and contributors"
GITHUB = "https://github.com/ARLQCI/C13NV.jl"

# The SVG figure files in docs/src/assets/ (levels.svg, liouville_levels.svg)
# are pre-rendered from TikZ sources in notes/images/. To update them:
#
#   cd notes/images
#   make              # builds levels.pdf and liouville_levels.pdf
#   pdf2svg levels.pdf ../../docs/src/assets/levels.svg
#   pdf2svg liouville_levels.pdf ../../docs/src/assets/liouville_levels.svg

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"))

mathengine = Documenter.MathJax3(
    Dict(
        :loader => Dict("load" => ["[tex]/physics"]),
        :tex => Dict(
            :macros => Dict(
                :tr => "\\textrm{tr}",
                :diag => "\\textrm{diag}",
                :RWA => "\\textrm{RWA}",
                :lab => "\\textrm{lab}",
                :Hilbert => "\\mathcal{H}",
                :HilbertS => "\\mathcal{H}_{\\!S}",
                :HilbertI => "\\mathcal{H}_{\\!I}",
                :HilbertGE => "\\mathcal{H}_{\\!GE}",
                :HilbertO => "\\mathcal{H}_{\\!O}",
                :HilbertOS => "\\mathcal{H}_{\\!OS}",
                :HilbertM => "\\mathcal{I}_{\\!M}",
                :HilbertG => "\\mathcal{I}_{\\!G}",
                :HilbertE => "\\mathcal{I}_{\\!E}",
                :ket => ["|#1\\rangle", 1],
                :bra => ["\\langle #1|", 1],
                :ketbra => ["|#1\\rangle\\!\\langle #2|", 2],
            ),
            :inlineMath => [["\$", "\$"], ["\\(", "\\)"]],
            :tags => "ams",
            :packages => ["base", "ams", "autoload", "physics"],
        ),
    ),
)

println("Starting makedocs")

PAGES = [
    "Home" => "index.md",
    "Hamiltonian" => "hamiltonian.md",
    "API" => "api.md",
    "References" => "references.md"
]

makedocs(;
    authors = AUTHORS,
    sitename = "$NAME.jl",
    format = Documenter.HTML(;
        prettyurls = true,
        canonical = "https://arlqci.github.io/C13NV.jl",
        edit_link = "master",
        footer = "[$NAME.jl]($GITHUB) v$VERSION docs powered by [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl).",
        assets = String["assets/citations.css"],
        mathengine = mathengine,
    ),
    pages = PAGES,
    plugins = [bib],
)

println("Finished makedocs")

deploydocs(; repo = "github.com/ARLQCI/C13NV.jl", devbranch = "master", push_preview = true)
