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