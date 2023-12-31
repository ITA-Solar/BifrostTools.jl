"""
Script for converting from Bifrost's simulation units to cgs or si units. 
Consult https://github.com/ITA-Solar/Bifrost/blob/develop/IDL/util/br_make_fits.pro
"""

# Simulation units
const params = Dict(
        "u_l" => 1f8,
        "u_t" => 1f2,
        "u_r" => 1f-7,
        "u_u" => 1f6,
        "u_p" => 1f5,
        "u_e" => 1f5,
        "u_B" => Float32(1121)
        )

# Converts from simulation units to cgs units
const cgs_params = Dict(
        # pressure
        "p" => params["u_p"],
        # gas density
        "r" => params["u_r"],
        # momentum
        "px" => params["u_r"]*params["u_u"],
        "py" => params["u_r"]*params["u_u"],
        "pz" => params["u_r"]*params["u_u"],
        # energy
        "e" => params["u_e"],
        # B-field
        "bx"  => params["u_B"],
        "by"  => params["u_B"],
        "bz"  => params["u_B"]
        )

# Converts from simulation units to cgs units
const si_params = Dict(
        # pressure
        "p" => params["u_p"]*1f-1,
        # gas density
        "r" => params["u_r"]*1f3,
        # momentum
        "px" => params["u_r"]*params["u_u"]*1f1,
        "py" => params["u_r"]*params["u_u"]*1f1,
        "pz" => params["u_r"]*params["u_u"]*1f1,
        # energy
        "e" => params["u_e"]*1f-1,
        # B-field
        "bx"  => params["u_B"]*1f-4,
        "by"  => params["u_B"]*1f-4,
        "bz"  => params["u_B"]*1f-4
        )

"""
    convert_units!(
        snapvariable, 
        variable::String, 
        unit_conversion::String
        )

Converts the units of `snapvariable` from simulation units to 'si' or 'cgs'
units
"""
function convert_units!(
    snapvariable::AbstractArray, 
    variable::String, 
    unit_conversion::String
    )

    if unit_conversion=="cgs"
        if variable != "tg"
            snapvariable .*= cgs_params[variable]
        end
        
    elseif unit_conversion=="si"
        if variable != "tg"
            snapvariable .*= si_params[variable]
        end
    
    else
        throw(ErrorException("Unit conversion '$unit_conversion' does not exits"))
    end
end

function convert_units(
    snapvariable::AbstractArray, 
    variable::String, 
    unit_conversion::String
    )::Array{Float32, 3}

    if unit_conversion=="cgs"
        if variable != "tg"
            snapvariable = snapvariable .* cgs_params[variable]
        end
        
    elseif unit_conversion=="si"
        if variable != "tg"
            snapvariable = snapvariable .* si_params[variable]
        end
    
    else
        throw(ErrorException("Unit conversion '$unit_conversion' does not exits"))
        snapvariable = snapvariable .* 0f0
    end
end

"""
    convert_time(t::AbstractFloat)

Converts snapshot time to seconds
"""
function convert_snaptime(t::AbstractFloat)::Float64

    t*params["u_t"]

end