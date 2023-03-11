push!(LOAD_PATH,"../Mstar2t")
import Documenter
using Mstar2t
using Mstar2t: Transport

Documenter.makedocs(modules  = Module[Mstar2t, Transport],
    format = Documenter.HTML(
        canonical = "https://marcofornari.github.io/etrasport/stable",
        assets = ["assets/icon.ico"],
        sidebar_sitename = false
    ),
    source = "src",
    build = "build",
    clean = true,
    doctest = true,
    sitename = "m*2T",
    highlightsig = true,
    expandfirst = [],
    authors = "Luca Bonaldo, Terry Ethen Stearns, Ilaria Siloi, Nicholas Mecholsky, Marco Fornari.",
    pages = [
        "Home" => "index.md",
        "Installation" => "install.md",
        "Usage" => "usage.md",
        "Documentation" => "docfun.md",
        "Examples" => "examples.md",
        ]
    )

Documenter.deploydocs(
    repo = "github.com/marcofornari/etrasport.git",
)
