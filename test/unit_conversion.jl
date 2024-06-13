# Solar units
u_l =  1e8
u_t =  1e2
u_r =  1e-7
u_p =  1e5
u_u =  1e6
u_kr =  1e1
u_ee =  1e12
u_e =  1e5
u_te =  1e11
u_B =  1.121e3

@test 1.0*u_t == BifrostTools.convert_timeunits(1.0,params)

@testset "unit conversions" begin
    # Simplest test data
    data = [1.0]

    conversions = Dict(
        "e" => u_e,
        "px" => u_r * u_u,
        "py" => u_r * u_u,
        "pz" => u_r * u_u,
        "r" => u_r,
        "bx" => u_B, 
        "by" => u_B, 
        "bz" => u_B,
        "p" => u_p,
        "tg" => 1e0,
        "ex" => u_u * u_B / 2.99792458e10,
        "ey" => u_u * u_B / 2.99792458e10,
        "ez" => u_u * u_B / 2.99792458e10,
        "qvisc" => u_e / u_t,
        "qjoule" => u_e / u_t,
        "qpdv" => u_e / u_t,
        "qrdiff" => u_e / u_t,
        "qediff" => u_e / u_t,
        "qeadv" => u_e / u_t
    )
 
    @test all([
        BifrostTools.convert_units(data, key, params, "cgs")[1] == value 
        for (key, value) in conversions
    ])

end