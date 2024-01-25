

"""
    read_params(file_name::String)

Reads and returns parameters `params` of a Bifrost simulation snapshot given 
the path `file_name` to the simulation snapshot. The input file should have the
format 'name_xxx.idl' where 'name' is the simulation name and 'xxx' is the 
snapshot number
"""
function read_params(file_name::String)
    
    params = Dict{String,String}()
  
    open(file_name, "r") do file
        for line in eachline(file)
            line = strip(line)
            if !isempty(line) && line[1] ≠ ';'
                line = replace(line, "\"" => "")
                line = replace(line, "'" => "")
                key, val = split(strip(line), '=')
                params[strip(key)] = strip(val)
            end
        end
    end

    # special case for the ixy1 variable, lives in a separate file
    if "aux" in keys(params)
        params["aux"] = replace(params["aux"], " ixy1" => "")
    end
    
    return params
end

"""
    read_params(expname::String, snap::Integer, expdir::String)
"""
function read_params(
    expname::String,
    snap::Integer,
    expdir::String
    )

    isnap = lpad(snap,3,"0")
    idl_file = string(joinpath(expdir,expname),"_",isnap,".idl")
    
    read_params(idl_file)
end

"""
    get_snap(
        expname::String,
        snap   ::Int,
        expdir ::String,
        )

Reads Bifrost *.snap binary file as an Array in dimension: (mx,my,mz,nvar).
Takes experiment name, experiment directory, and snap number as arguments.
Returns `snapdata` (the data) and `params` (the snap parameters).
Assumes single floating point precision by default.

Variables of `snapdata`:

        snapdata[:,:,:,1] : r, density
        snapdata[:,:,:,2] : px, x-component of momentum 
        snapdata[:,:,:,3] : py, y-component of momentum 
        snapdata[:,:,:,4] : pz, z-component of momentum 
        snapdata[:,:,:,5] : e, energy
    
        if params["do_mhd"] == 1 # numvars = 8, otherwise numvars = 5
            snapdata[:,:,:,6] : bx, x-component of magnetic field 
            snapdata[:,:,:,7] : by, y-component of magnetic field 
            snapdata[:,:,:,8] : bz, z-component of magnetic field 

Warning:
    variables in code units.ts
"""
function get_snap(
    expname  ::String,
    snap     ::Int,
    expdir   ::String,
    precision::DataType=Float32
    )
    # Parse filenames
    basename = string(joinpath(expdir, expname),"_$(lpad(snap,3,"0"))")
    idl_filename = string(basename, ".idl")
    snap_filename = string(basename, ".snap")
    params = read_params(idl_filename)
    return get_snap(snap_filename, params, precision), params
end


"""
    get_snap(
        file_name::String,
        params   ::Dict{String,String}
        )
"""
function get_snap(
    file_name::String,
    params   ::Dict{String,String},
    precision::DataType=Float32
    )
    datadims = 4 # 3 spatial dimensions and 1 variable dimension
    
    snapsize, numvars, _ = get_snapsize_and_numvars(params)

    file = open(file_name)
    # Use Julia standard-library memory-mapping to extract file values
    snapdata = mmap(file, Array{precision, datadims}, (snapsize..., numvars))
    close(file)
    return snapdata
end

"""
    get_aux(
        file_name::String,
        params   ::Dict{String,String}
    )
Reads Bifrost *.aux binary file using memory-mapping. The returned
`auxdata` array will have dimensions (mx,my,mz,nvars) where nvars is the
number of aux-variables. Assumes single floating point precision by default.
"""
function get_aux(
    file_name::String,
    params   ::Dict{String,String},
    precision::DataType=Float32
    )
    datadims = 4
    snapsize, _, numauxvars = get_snapsize_and_numvars(params)
    if numauxvars == 0
        return
    else 
        file = open(file_name)
        # Use Julia standard-library memory-mapping to extract file values
        auxdata = mmap(file, Array{precision, datadims},
                       (snapsize..., numauxvars)
                       )
        close(file)
        return auxdata
    end
end

"""
    get_var(
        xp::BifrostExperiment,
        snap::Union{<:Integer, AbstractVector{<:Integer}},
        variable::String,
        args...
        ;
        kwargs...
        )
Load a `variable` from one or multiple snapshots of `xp`.

# Available variables

The primary variables:
- "r":  density
- "px": x-component of momentum
- "py": y-component of momentum
- "pz": z-component of momentum
- "e":  energy

`if params["do_mhd"] == true`
- "bx": x-component of magnetic field
- "by": y-component of magnetic field
- "bz": z-component of magnetic field

auxilliary variables (variables in params["aux"]):
- "p": pressure
- "tg": gas temperature
    ...

# Optional keyword-arguments
Converts variables to "si" or "cgs" units: `units="si"` or `units="cgs"`.

To load a slice of the variable, give e.g. `slicex=[32, 410]` or `slicey=40:90`

# Example usage: 

```{julia}
exp_name = "cb24oi"
exp_dir = "/mn/stornext/d21/RoCS/matsc/3d/run/cb24oi"
snap = 700

xp = BifrostExperiment(expname, expdir)

# Load pressude for the full cube in si units
pressure = get_var(xp, snap, "p"; units="si")

# Load gas density in a slize along the xy-plane in cgs units
rho = get_var(xp, snap, "r"; units="cgs", slicez=[100])
```
"""
function get_var(
    xp      ::BifrostExperiment,
    snap    ::Union{<:Integer, AbstractVector{<:Integer}},
    variable::String,
    args...
    ;
    kwargs...
    )
    
    get_var(xp.expname, snap, xp.expdir, variable, args...; kwargs...)
end

"""
    get_var(
        expname ::String,
        snap    ::Union{<:Integer, AbstractVector{<:Integer}},
        expdir  ::String,
        variable::String,
        args...
        ;
        kwargs...
    )

Load a `variable` from one or multiple snapshots of a Bifrost experiment with 
experiment directory `expdir` and experiment name `expname`.
"""
function get_var(
    expname ::String,
    snaps   ::Union{<:Integer, AbstractVector{<:Integer}},
    expdir  ::String,
    variable::String,
    args...
    ;
    kwargs...
    )

    if variable == "t"
        # The special case of getting the snapshot time
        return get_time(expname,snaps,expdir;kwargs...)
    end

    kwarg_keys = keys(kwargs)
    kwarg_values = values(kwargs)

    # Allocate space for variable
    var = Vector{Array{precision, 3}}(undef, length(snaps))

    # decide which get_var function to use
    if :destagger in kwarg_keys && kwarg_values.destagger
        
        if variable ∉ ("px","py","pz","bx","by","bz","ex","ey","ez","ix","iy","iz")
            @warn "Variable $variable is not usually staggered. "*
                "Destagger will take place anyways."
            get_function = get_and_destagger_var
        elseif variable ∈ ("ex","ey","ez","ix","iy","iz")
            @warn "Destagger of $variable is not implemented. "*
                "Defaulting to loading variable without destaggering"
            get_function = get_var
        else
            get_function = get_and_destagger_var
        end

        # If destagger direction is not passed fall back to default direction
        if ! :direction in kwarg_keys
            kwargs[:direction] = destagger_direction(variable)
        end

    else
        get_function = get_var
    end
    
    # Loop over snapshots
    Threads.@threads for (i,snap) in collect(enumerate(snaps))
        
        params_local = read_params(expname,snap,expdir)
        varnr, file_ext = get_varnr_and_file_extension(params, variable)
        tmp_file = string(basename, "_", isnap_local, file_ext)

        var[i] = get_function(
            tmp_file,
            params_local,
            varnr,
            ;
            precision=precision,
            kwargs...)

    end

    # -------------------------------------------------------------------------
    #  Below is where you extend the functionality of get_var by handling 
    #  arbitrary keyword arguments. Please put your new implementation into a 
    #  new function.
    # -------------------------------------------------------------------------

    # UNITS: Scale from code units to something else
    # 
    #   If multiple snapshots: Assumes the same conversion factor for all
    #
    if :units in kwarg_keys
        data = convert_units(data, variable, params, kwarg_values.units)
    end

    #
    # ORIENTATION: Rotate coordinate system
    #
    if :rotate_about_x in kwarg_keys && kwarg_values.rotate_about_x
        data = rotate(data, variable, "x")
    end

    # CONCATENATION: Concatenate vector to 4D array
    if :stack in kwarg_keys && kwarg_values.stack
       data = stack(data)
    end

    # SQUEEZE: Drop empty dimensions
    if length(snaps) == 1 && :squeeze in kwarg_keys && kwarg_values.squeeze
        data = data[1]
    end

    # -------------------------------------------------------------------------
    # Additional kwargs and functionality go under here
    # -------------------------------------------------------------------------
    return data
end

"""
    get_var(
        filename       ::String,
        params         ::Dict{String,String},
        varnr          ::Integer,
        precision      ::DataType=Float32;
        slicex         ::AbstractVector{<:Integer}=Int[],
        slicey         ::AbstractVector{<:Integer}=Int[],
        slicez         ::AbstractVector{<:Integer}=Int[]
        )
Load variable nr. `varnr` from `filename`. The variable could be either 
primary or auxiliary. Slicing the snapshot is optional. Assumes single 
precision snapshot by default.
"""
function get_var(
    filename       ::String,
    params         ::Dict{String,String},
    varnr          ::Integer
    ;
    precision      ::DataType=Float32,
    slicex         ::AbstractVector{<:Integer}=Int[],
    slicey         ::AbstractVector{<:Integer}=Int[],
    slicez         ::AbstractVector{<:Integer}=Int[],
    kwargs...
    )
    datadims = 3 # 3 spatial dimensions and 1 variable dimension
    snapsize = get_snapsize(params)
    # Do slicing or not (returns the mmap)
    if isempty(slicex) && isempty(slicey) && isempty(slicez)
        slicing = false
    else
        slicing = true
    end
    # Calculate offset in file
    offset = get_variable_offset_in_file(precision, snapsize, varnr)
    file = open(filename)
    # Use Julia standard-library memory-mapping to extract file values
    if slicing
        isempty(slicex) && ( slicex = 1:snapsize[1] )
        isempty(slicey) && ( slicey = 1:snapsize[2] )
        isempty(slicez) && ( slicez = 1:snapsize[3] )
        # slice the variable
        data = mmap(file,
            Array{precision, datadims}, 
            snapsize, 
            offset
            )[slicex,slicey,slicez]
    else
        # do not slice the variable
        data = mmap(file,
            Array{precision, datadims}, 
            snapsize, 
            offset
            )
    end
    close(file)
    return data
end

function get_time(
    expname::String,
    snap   ::Union{<:Integer, AbstractVector{<:Integer}},
    expdir ::String 
    ;
    kwargs...
    )
    nsnaps = length(snap)
    if nsnaps == 1
        params = read_params(expname,snap,expdir)
        var = [parse(Float64, params["t"])]
    else 
        var = Vector{Float64}(undef, nsnaps)
        # Load the variable directly from params
        for (i,snap) in enumerate(snap)
            params = read_params(expname,snap,expdir)
            time = parse(Float64, params["t"])
            var[i] = time
        end
    end

    if :units in keys(kwargs)
        convert_timeunits!(var, params)
    end
    
    return var
end

"""
    function get_and_destagger_var(expname::String,
        snap::Integer,
        expdir::String,
        variable::String
        ;
        precision::DataType=Float32,
        units::String="none",
        direction::String="zup",
        periodic::Bool=false,
        order::Int=6,
        slicex::AbstractVector{<:Integer}=Int[],
        slicey::AbstractVector{<:Integer}=Int[],
        slicez::AbstractVector{<:Integer}=Int[]
        )

Function to load a staggered variable and interpolate it to cell center.
The staggered variables that typically need to be interpolated are the velocity
and magnetic field components. Normally you need to use `direction="zup"` for 
vz and bz with `periodic=false` (these are the default arguments), and 
`direction="xup"` for vx and bx with `periodic=true` (same for y direction).
"""
function get_and_destagger_var(
    filename::String,
    params::Dict{String,String},
    varnr::Integer,
    ;
    direction::String="zup",
    periodic::Bool=false,
    order::Int=6,
    slicex::AbstractVector{<:Integer}=Int[],
    slicey::AbstractVector{<:Integer}=Int[],
    slicez::AbstractVector{<:Integer}=Int[],
    kwargs...
    )

    shift_functions = Dict(
        "xdn" => xdn, 
        "xup" => xup,
        "ydn" => ydn, 
        "yup" => yup, 
        "zdn" => zdn, 
        "zup" => zup
    )

    shift = shift_functions[direction]

    slicing = true
    ( isempty(slicex) && isempty(slicey) && isempty(slicez) ) && ( slicing = false )
    
    if slicing
        if ( direction == "xup" ) || ( direction == "xdn" )
            # Load var
            var = get_var(filename,params,varnr,slicey=slicey,slicez=slicez)
            if isempty(slicex)
                # All indices in 'x' are loaded, don't worry about slicing
                var = shift(var,periodic,order)
            else
                # Call to the function that slices in x
                var = shift(var,slicex,periodic,order)
            end
        elseif ( direction == "yup" ) || ( direction == "ydn" )
            var = get_var(filename,params,varnr,slicex=slicex,slicez=slicez)
            if isempty(slicey)
                # All indices in 'y' are loaded, don't worry about slicing
                var = shift(var,periodic,order)
            else
                # Call to the function that slices in y
                var = shift(var,slicey,periodic,order)
            end
        else
            var = get_var(filename,params,varnr,slicex=slicex,slicey=slicey)
            if isempty(slicez)
                # All indices in 'z' are loaded, don't worry about slicing
                var = shift(var,periodic,order)
            else
                # Call to the function that slices in z
                var = shift(var,slicez,periodic,order)
            end
        end
    else
        # load the entire variable and shift it in the desired direction
        var = get_var(filename,params,varnr,)
        var = shift(var,periodic,order)
    end

    return var

end


"""
    destagger(data::AbstractArray, variable::String)
De-stagger the Bifrost `data` of type `variable` to cell centre. The input may 
be read-only, so return a copy
"""
function destagger(
    data    ::AbstractArray{<:Real, 3},
    variable::String,
    )
    if variable in ("r", "e", "tg", "p")
        return data # nothing to do, already cell centred
    elseif variable in ("px", "bx")
        return xup(data)
    elseif variable in ("py", "by")
        return yup(data)
    elseif variable in ("pz", "bz")
        return zup(data)
    elseif variable in ("ex", "ix")
        return yup(zup(data)) # not 100% sure about this operation
    elseif variable in ("ey", "iy")
        return zup(xup(data)) # not 100% sure about this operation
    elseif variable in ("ez", "iz")
        return xup(yup(data)) # not 100% sure about this operation
    else
        error("Destaggering of variable $variable is not implemented.")
    end
end

"""
    destagger_direction(variable::String)

Takes in a variable, and determines the default direction to destagger the variable
"""
function destagger_direction(variable::String)
    if variable in ("px", "bx")
        return "xup"
    elseif variable in ("py", "by")
        return "yup"
    elseif variable in ("pz", "bz")
        return "zup"
    else
        error("Give direction manually to destagger variable $variable.")
    end
end


"""
    rotate(
        data         ::AbstractArray,
        variable     ::String,
        rotation_axis::String,
        )
Rotate the data about an `rotation_axis`.
"""
function rotate(
    data         ::AbstractArray,
    variable     ::String,
    rotation_axis::String,
    )
    xcomponents = ("px", "bx", "ex", "ix")
    ycomponents = ("py", "by", "ey", "iy")
    zcomponents = ("pz", "bz", "ez", "iz")

    if variable in ("r", "p", "tg")
        return data # Scalar fields, do nothing
    else
        if rotation_axis == "x"
            if variable in xcomponents
                return data
            elseif variable in ycomponents || variable in zcomponents
                return -data
            end
        else
            error("Rotation about $rotation_axis-axis is not implemented")
        end
    end
end


"""
    function get_electron_density(
        xp::BifrostExperiment,
        snap::Integer,
        kwargs...)

Function to calculate the electron density from a snapshot `snap`. Supports
slicing. Gas density `rho` and internal energy `e` are optional arguments and
can be passed (but they MUST be in cgs units). If these quantities already 
exist, passing them will speed up the calculation of electron density.

`kwargs`:
    units::String="si",
    slicex::AbstractVector{<:Integer}=Int[],
    slicey::AbstractVector{<:Integer}=Int[],
    slicez::AbstractVector{<:Integer}=Int[],
    rho::Array{AbstractFloat,3}=Float32[;;;],
    e::Array{AbstractFloat,3}=Float32[;;;],
    tabfile::String="tabparam.in"
"""
function get_electron_density(
    xp::BifrostExperiment,
    snaps::Union{<:Integer, AbstractVector{<:Integer}};
    slicex::AbstractVector{<:Integer}=Int[],
    slicey::AbstractVector{<:Integer}=Int[],
    slicez::AbstractVector{<:Integer}=Int[],
    stack::Bool=true,
    kwargs...)

    if typeof(snaps) <: Integer
        var = get_electron_density(xp.expname,snaps,xp.expdir;
		    slicex=slicex,slicey=slicey,slicez=slicez,kwargs...)
        return var

    elseif typeof(snaps) <: AbstractVector{<:Integer}
        var = Vector{Array{Float32,3}}(undef, length(snaps))
        Threads.@threads for (i,snap) in collect(enumerate(snaps))
            var[i] = get_electron_density(xp.expname,snap,xp.expdir;
                        slicex=slicex,slicey=slicey,slicez=slicez,kwargs...)
        end

        if stack
            return stack(var)
        else
            return var
        end
    end

end

"""
    function get_electron_density(
        expname::String,
        snap::Integer,
        expdir::String;
        units::String="si",
        slicex::AbstractVector{<:Integer}=Int[],
        slicey::AbstractVector{<:Integer}=Int[],
        slicez::AbstractVector{<:Integer}=Int[],
        rho::Array{T,3}=Float32[;;;],
        e::Array{T,3}=Float32[;;;],
        tabfile::String="tabparam.in"
        ) where {T<:AbstractFloat}
"""
function get_electron_density(
    expname::String,
    snap::Integer,
    expdir::String;
    units::String="si",
    slicex::AbstractVector{<:Integer}=Int[],
    slicey::AbstractVector{<:Integer}=Int[],
    slicez::AbstractVector{<:Integer}=Int[],
    rho::Array{T,3}=Float32[;;;],
    e::Array{T,3}=Float32[;;;],
    tabfile::String="tabparam.in"
    ) where {T<:AbstractFloat}

    basename, basename_isnap = filename(expname, snap, expdir)
    params = read_params(string(basename_isnap, ".idl"))
    
    # rho in g/cm^3
    if isempty(rho)
        
        varnr, file_ext = get_varnr_and_file_extension(params, "r")
        tmp_file = string(basename, "_", lpad(snap,3,"0"), file_ext)

        rho = get_var(
            tmp_file,
            params,
            varnr,
            slicex=slicex,
            slicey=slicey,
            slicez=slicez
        )
        rho = convert_units(rho, "r", params, "cgs")

    end
    
    # internal energy in ergs
    if isempty(e)
        varnr, file_ext = get_varnr_and_file_extension(params, "e")
        tmp_file = string(basename, "_", lpad(snap,3,"0"), file_ext)

        e = get_var(
            tmp_file,
            params,
            varnr,
            slicex=slicex,
            slicey=slicey,
            slicez=slicez
        )
        e = convert_units(e, "e", params, "cgs")

    end

    # Calculate internal energy per mass
    ee = e ./ rho
    
    # construct the EOS tables for interpolation of electron density
    tabfile = joinpath(expdir,tabfile)
    eos = EOSTables(tabfile)

    if maximum(rho) > parse(Float64,eos.params["RhoMax"])
        @printf """tab_interp: density outside table bounds.
        Table rho max=%.3e, requested rho max=%.3e""" eos.params["RhoMax"] maximum(rho)
    end        
    if minimum(rho) <parse(Float64,eos.params["RhoMin"])
        @printf """tab_interp: density outside table bounds.
        Table rho min=%.3e, requested rho min=%.3e""" eos.params["RhoMin"] minimum(rho)
    end
    
    if maximum(ee) > parse(Float64,eos.params["EiMax"])
        @printf """tab_interp: density outside table bounds.
        Table rho max=%.3e, requested rho max=%.3e""" eos.params["EiMax"] maximum(Ei)
    end        
    if minimum(ee) < parse(Float64,eos.params["EiMin"])
        @printf """tab_interp: density outside table bounds.
        Table rho min=%.3e, requested rho min=%.3e""" eos.params["EiMin"] minimum(Ei)
    end
    
    # Create interpolation table, takes the log of coordinates
    itp_table = eos_interpolate(eos,3)

    x = log.(ee)
    y = log.(rho)

    ne = itp_table.(x, y)

    # take exp to remove log
    ne = exp.(ne)

    # Convert to si on request (cm^-3 --> m^-3)
    if lowercase(units) == "si"
        ne .*= 1f6
    end

    return ne
end
