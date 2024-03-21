using BifrostTools
using Documenter

DocMeta.setdocmeta!(BifrostTools, :DocTestSetup, :(using BifrostTools); recursive=true)

makedocs(;
    modules=[BifrostTools],
    authors="meudnaes <eliasudnes@hotmail.com>, eilifso <eilif.oyre@gmail.com> and contributors",
    sitename="BifrostTools.jl",
    format=Documenter.HTML(;
        canonical="https://ITA-Solar.github.io/BifrostTools.jl",
        edit_link="develop-documentation",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Installation" => "install.md",
        "Example Usage" => "usage.md",
        "Documentation" => "documentation.md"
    ],
)

deploydocs(;
    repo="github.com/ITA-Solar/BifrostTools.jl",
    devbranch="develop-documentation",
)
