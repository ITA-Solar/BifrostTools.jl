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

    test_cgs_to_si_conversions = Dict(
    # Pressure:    g/s^2/cm * 1f-3 kg/g * 1f2 cm/m = 1f-1 kg/s^2/m
    "p"  => 1e-1,
    # Gas density: g/cm^3 * 1f-3 kg/g * 1f6 cm^3/m^3 = 1f3 kg/m^3
    "r"  => 1e3,
    # Momentum:    g/cm^2/s * 1f-3 kg/g * 1f4 cm^2/m^2 = 1f1 kg/m^2/s
    "px" => 1e1,
    "py" => 1e1,
    "pz" => 1e1,
    # Bulk velocity:  cm/s * 1f-2 m/cm = 1f-2 m/s
    "ux" => 1e-2,
    "uy" => 1e-2,
    "uz" => 1e-2,
    # Internal energy:     erg/cm^3 * 1f-7 J/erg * 1f6 cm^3/m = 1f-1 J/m^3
    "e"  => 1e-1,
    # Dissipation coefficients/energy terms:     erg/cm^3/s * 1.e-7 J/erg * 1.e6 cm^3/m^3 = W/m^3
    "qvisc" => 1e-1,
    "qjoule" => 1e-1,
    "qpdv" => 1e-1,
    "qrdiff" => 1e-1,
    "qediff" => 1e-1,
    "qeadv" => 1e-1,
    # Magnetic field: G * 1f-4 T/G = 1f-4 T
    "bx" => 1e-4,
    "by" => 1e-4,
    "bz" => 1e-4,
    # Electric field: statV/cm * 1f-2 m/cm * 1f-4 T/G * c[cm/s] = 2.998f4 V/m
    "ex" => 2.99792458e4,
    "ey" => 2.99792458e4,
    "ez" => 2.99792458e4,
    # Temperature: K = K
    "tg" => 1.0,
    )

 
    @test all([
        BifrostTools.convert_units(data, key, params, "cgs")[1] == value 
        for (key, value) in conversions
    ])

    @test all([
        BifrostTools.convert_units(data, key, params, "code") == data
        for key in keys(conversions)
    ])

    @test all([
        BifrostTools.convert_units(data, key, params, "si")[1] == value*test_cgs_to_si_conversions[key]
        for (key, value) in conversions
    ])

end
