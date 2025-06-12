using Aqua
using Test
using BifrostTools
using Interpolations

Aqua.test_all(BifrostTools)

BASE_FOLDER = dirname(dirname(pathof(BifrostTools)))

expdir = joinpath(BASE_FOLDER,"test","sp.n064")
expname = "en48"

xp = BifrostExperiment(expname,expdir)
params = read_params(expname,xp.snaps[1],expdir)

@testset verbose = true "Total" begin
    include("experiment.jl")
    include("read_params_snap_aux.jl")
    include("stagger_operators.jl")
    include("unit_conversion.jl")
    include("eos_tables.jl")
    include("utils.jl")
end
