"""
Script for converting from Bifrost's simulation units to cgs or si units. 
Based on https://github.com/ITA-Solar/Bifrost/blob/develop/IDL/util/br_make_fits.pro
"""


"""
    cgs_to_SI_conversion_factors
Factors for converting some physical quantities from cgs-units to SI-units.
"""
const cgs_to_SI_conversion_factors = Dict(
    # Pressure:    g/s^2/cm * 1f-3 kg/g * 1f2 cm/m = 1f-1 kg/s^2/m
    "p"  => 1f-1,
    # Gas density: g/cm^3 * 1f-3 kg/g * 1f6 cm^3/m^3 = 1f3 kg/m^3
    "r"  => 1f3,
    # Momentum:    g/cm^2/s * 1f-3 kg/g * 1f4 cm^2/m^2 = 1f1 kg/m^2/s
    "px" => 1f1,
    "py" => 1f1,
    "pz" => 1f1,
    # Internal energy:     erg/cm^3 * 1f-7 J/erg * 1f6 cm^3/m = 1f-1 J/m^3
    "e"  => 1f-1,
    # Dissipation coefficients:     erg/cm^3/s * 1.e-7 J/erg * 1.e6 cm^3/m^3 = W/m^3
    "qvisc" => 1f-1,
    "qjoule" => 1f-1,
    # Magnetic field: G * 1f-4 T/G = 1f-4 T
    "bx" => 1f-4,
    "by" => 1f-4,
    "bz" => 1f-4,
    # Electric field: g*cm/s^2/Fr * 1f-3 kg/g * 1f-2 m/cm * 1f-1 Fr/C 
    #                   = 1f-6 kg*m/s^2/C
    "ex" => 1f-6,
    "ey" => 1f-6,
    "ez" => 1f-6,
    )


"""
    convert_units(
        data    ::AbstractArray,
        variable::String,
        params  ::Dict{String,String},
        units   ::String,
        )
Convert the `data` from code `units` to someting else.
"""
function convert_units(
    data    ::AbstractArray,
    variable::String,
    params  ::Dict{String,String},
    units   ::String,
    )
    if lowercase(units) == "si"
        return code_to_SI(data, variable, params)
    elseif lowercase(units) == "cgs"
        return code_to_cgs(data, variable, params)
    elseif lowercase(units) == "code"
        # Do nothing
        nothing
    else
        throw(ErrorException("Unit conversion '$units' does not exits"))
    end
end

"""
    code_to_SI(
        data    ::AbstractArray,
        variable::String,
        params  ::Dict{String,String},
    )
Convert the `data` from code units to SI units.
"""
function code_to_SI(
    data    ::AbstractArray,
    variable::String,
    params  ::Dict{String,String},
    )
    tmp = code_to_cgs(data, variable, params)
    return cgs_to_SI(tmp, variable)
end

"""
    code_to_cgs(
        data    ::AbstractArray,
        variable::String,
        params  ::Dict{String,String},
    )
Convert the `data` from code-units to cgs-units.
"""
function code_to_cgs(
    data    ::AbstractArray,
    variable::String,
    params  ::Dict{String,String},
    )
    if variable == "r"                       # Density
        return  data * parse(Float32, params["u_r"])
    elseif variable == "e"                   # Energy
         return data * parse(Float32, params["u_e"])
    elseif variable == "tg"                  # Gas temperature
        return data # nothing to do
    elseif variable == "p"                   # Pressure
         return data * parse(Float32, params["u_p"])
    elseif variable in ("px", "py", "pz")    # Momentum
         return data*parse(Float32, params["u_r"])*parse(Float32,params["u_u"])
    elseif variable in ("bx", "by", "bz")    # Magnetic field
         return data * parse(Float32, params["u_B"])
    #elseif variable in ("ix", "iy", "iz")    # Current density
        # not implemented yet
    elseif variable in ("ex", "ey", "ez")    # Electric field
         return data*parse(Float32, params["u_u"])*parse(Float32,params["u_B"])
    elseif variable in ("qvisc, qjoule")
        return data*parse(Float32, params["u_e"])/parse(Float32, params["u_t"]) 
    else
        throw(ErrorException(
            "Conversion to cgs-units of variable $variable is not implemented."
            ))
    end
end

"""
    cgs_to_SI(
        data    ::AbstractArray,
        variable::String,
        params  ::Dict{String,String},
    )
Convert the `data` from cgs-units to SI-units.
"""
function cgs_to_SI(
    data    ::AbstractArray,
    variable::String,
    )
    if variable == "tg" # Temperature: Kelvin -> Kelvin
        return data # nothing to do
    elseif variable in keys(cgs_to_SI_conversion_factors)
        return data * cgs_to_SI_conversion_factors[variable]
    else
        throw(ErrorException(
            "Conversion to SI-units of variable $variable is not implemented."
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

