"""
# Overview
Script for converting from Bifrost's simulation units to cgs or SI units.
Based on https://github.com/ITA-Solar/Bifrost/blob/develop/IDL/util/br_make_fits.pro

Many auxiliary variables are missing conversion factors. Feel free to add.

## Note on the electric field
The conversion of electric field from cgs to SI is based on dimensional analysis of
Ohm's law. In cgs-units, Ohm's law reads
    ηJ = E + (u x B) / c.
In SI-units, one omits the lightspeed constant
    ηJ = E + u x B.
Hence the division/multiplicatoin by c during conversion.
"""


"""
    cgs_to_SI_conversion_factors
Factors for converting some physical quantities from cgs-units to SI-units.
"""
const cgs_to_SI_conversion_factors = Dict(
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
    # Hion densities: cm^-3 to m^-3
    "ne" => 1e6,
    "n1" => 1e6,
    "n2" => 1e6,
    "n3" => 1e6,
    "n4" => 1e6,
    "n5" => 1e6,
    "n6" => 1e6,
    "nh2" => 1e6,
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
    # Position
    "x" => 1e-2,
    "y" => 1e-2,
    "z" => 1e-2,
    )

const c_in_cgs = 2.99792458e10


"""
    convert_units(
        data    ::AbstractArray,
        variable::String,
        params  ::Dict{String,String},
        units   ::String,
        )
Convert the `data` from code `units` to cgs or SI. Conversion factor depends
on `variable` and snapshot `params`.
"""
function convert_units(
    data    ::AbstractArray,
    variable::String,
    params  ::Dict{String,String},
    units   ::String,
    )
    # Working floating precision
    # Data is a vector of snapshots and has type Vector{AbstractArray}.
    wfp = eltype(eltype(data))
    if lowercase(units) == "si"
        conversionfactor = code_to_cgs(variable, params)
        try conversionfactor *= cgs_to_SI_conversion_factors[variable]
        catch
            error(
            "Conversion to SI-units of variable $variable from CGS-units is" *
            " not implemented."
            )
        end
        return data * wfp(conversionfactor)
    elseif lowercase(units) == "cgs"
        return data * wfp(code_to_cgs(variable, params))
    elseif lowercase(units) == "code"
        # Do nothing
        return data
    else
        throw(ErrorException("Unit conversion '$units' does not exits"))
    end
end


"""
    code_to_cgs(
        variable::String,
        params  ::Dict{String,String},
    )
Conversion factor of `variable` from code units to cgs units.
"""
function code_to_cgs(
    variable::String,
    params  ::Dict{String,String},
    )
    if variable == "r"                       # Density
        return  parse(Float64, params["u_r"])
    elseif variable == "e"                   # Energy
         return parse(Float64, params["u_e"])
    elseif variable in ("tg", "hiontg")      # Gas temperature
        return 1.0 # nothing to do
    elseif variable in ("ne", "n1", "n2", "n3", "n4", "n5", "n6", "nh2")  # Hion densities
        return 1.0 # nothing to do, already in cgs
    elseif variable == "p"                   # Pressure
         return parse(Float64, params["u_p"])
    elseif variable in ("px", "py", "pz")    # Momentum
         return parse(Float64, params["u_r"])*parse(Float32,params["u_u"])
    elseif variable in ("bx", "by", "bz")    # Magnetic field
         return parse(Float64, params["u_B"])
    #elseif variable in ("ix", "iy", "iz")    # Current density
        # not implemented yet
    elseif variable in ("ex", "ey", "ez")    # Electric field
        u_u = parse(Float64, params["u_u"])
        u_B = parse(Float64, params["u_B"])
        return u_u * u_B / c_in_cgs
    elseif variable in ("qvisc", "qjoule", "qpdv", "qrdiff", "qediff", "qeadv")
        return parse(Float64, params["u_e"])/parse(Float64, params["u_t"])
    else
        throw(ErrorException(
            "Conversion to cgs-units of variable $variable is not implemented."
            ))
    end
end

"""
    convert_timeunits!(
    t     ::AbstractArray,
    params::Dict{String,String}
    )     ::Float64

Converts snapshot time to seconds
"""
function convert_timeunits(
    t     ::Union{AbstractArray, AbstractFloat},
    params::Dict{String,String}
    )

    t *= parse(Float64, params["u_t"])
end

function convert_axesunits(
        mesh::BifrostMesh,
        params::Dict{String, String}
        ;
        units::String="code"
    )
    wfp = eltype(mesh.x)
    if units == "cgs"
        u_l = parse(Float64, params["u_l"])
        conversionfactor = wfp(u_l)
    elseif units == "si"
        u_l = parse(Float64, params["u_l"])
        conversionfactor = wfp(u_l*cgs_to_SI_conversion_factors["x"])
    elseif units == "code"
        conversionfactor = wfp(1.0)
    else
        throw(ErrorException("Unit conversion '$units' is not implemented"))
    end
    return mesh.x*conversionfactor,
        mesh.y*conversionfactor,
        mesh.z*conversionfactor
end
