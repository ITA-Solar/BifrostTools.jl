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

#deploydocs(;
#    repo="github.com/ITA-solar/BifrostTools.jl",
#)

deploydocs(
    repo = "github.com/ITA-solar/BifrostTools.jl.git",
    devbranch="develop-documentation",
)