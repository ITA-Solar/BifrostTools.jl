using BifrostTools
using JLD2
using Test
using StaticArrays
using Interpolations

@testset "get_itp" begin
BASE_FOLDER = dirname(dirname(pathof(BifrostTools)))

expdir = joinpath(BASE_FOLDER,"test","sp.n064")
expname = "en48"

xp = BifrostExperiment(expname,expdir)


    # Create and save a 3-component interpolation object of a 3D snapshot with
    # a non-uniform z-axis.
    tmpdir = mktempdir()
    name = "tmp_itp"
    itp1, metadata1 = get_itp(
        xp,
        1,
        "bx by bz",
        [:linear, Gridded(Linear()), "constant"];
        itp_bc = [:periodic, Periodic(), "throw"],
        units="si",
        slicez=17:48,
        saveitp=true,
        name = "tmp_itp",
        dir = tmpdir,
    )
    # Calling it with Float32 arguments should give Float32 results 
    vec32 = itp1(0f0, 1f7, -1f6)
    vec64 = itp1(0.0, 1f7, -1f6)

    # Test that the metadata contains correct information.
    @test metadata1[:size           ] == (48, 48, 32)
    @test metadata1[:snaps          ] == 1
    @test metadata1[:itp_name       ] == name
    @test metadata1[:call_signature ] == "(x, y, z)"
    @test metadata1[:vars           ] == ["bx", "by", "bz"]
    @test metadata1[:experiment_name] == "en48"
    @test metadata1[:units          ] == "si"
    @test metadata1[:bounds] == [
        (0f0, 2.38125f7),
        (0f0, 2.38125f7),
        (-5.688667f6, -554340.2f0)
    ]
    @test metadata1[:itp_bc] == "(Periodic(), Periodic(), Throw())"
    @test metadata1[:itp_type] ==
        "(Gridded(Linear()), Gridded(Linear()), Gridded(Constant{Nearest}()))"

    # Test the interpolation type
    @test typeof(itp1) <: AbstractInterpolation{
        SVector{3, Float32},
        3,
        Tuple{
        Gridded{Linear{Throw{OnGrid}}},
        Gridded{Linear{Throw{OnGrid}}},
        Gridded{Constant{Nearest, Throw{OnGrid}}}
        }
    }

    # Test the interpolation results
    @test vec32 isa SVector{3, Float32}
    @test vec64 isa SVector{3, Float64}
    @test vec32 == [0.07956992f0, -0.0004932007f0, 0.07995612f0]
    @test_throws BoundsError itp1(0, 0, 0)

    # This Should not throw due to periodic boundaries in x and y.
    @test itp1(-1, -1, -1e6) isa SVector

    # Make sure the interpolation object is stored correctly
    fname = joinpath(tmpdir, name*".itp.jld2")
    @test isfile(fname)
    @load fname itp metadata
    @test itp == itp1
    @test metadata == metadata1

    # Should not create cubic spline interpolation on non-uniform grid
    @test_throws ArgumentError get_itp(xp, 1, "r", :cubic_spline)

    # Test 4D and cubic spline interpolation
    xp_exhaust = BifrostExperiment(
        "fanspine_exhaust64",
        joinpath(@__DIR__, "fanspine_exhaust64")
    )
    itp2, metadata2 = get_itp(xp_exhaust, 1:5, "ey", :cubic_spline)
    @test typeof(itp2) <: AbstractInterpolation{
        Float32,
        4,
        BSpline{Cubic{Line{OnGrid}}}
    }
    @test metadata2[:size] == (64, 64, 64, 5)
    @test metadata2[:snaps] == 1:5
    @test metadata2[:bounds][4] == (6.0012865f0, 6.005295f0)
end
