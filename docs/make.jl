push!(LOAD_PATH,"../src/")

using Documenter, BifrostTools

makedocs(
    sitename="BifrostTools.jl Documentation",
    pages = [
        "Home" => "index.md",
        "Installation" => "install.md",
        "Example Usage" => "usage.md"
    ]
)

deploydocs(
    repo = "github.com/ITA-solar/BifrostTools.jl.git",
    branch="gh-pages",
    devbranch="develop-documentation",
    versions = nothing
)