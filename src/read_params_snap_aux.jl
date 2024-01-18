

"""
    br_read_params(file_name::String)

Reads and returns parameters `params` of a Bifrost simulation snapshot given 
the path `file_name` to the simulation snapshot. The input file should have the
format 'name_xxx.idl' where 'name' is the simulation name and 'xxx' is the 
snapshot number
"""
function br_read_params(file_name::String)
    
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
    basename = string(expdir, "/", expname, "_$(lpad(snap,3,"0"))")
    idl_filename = string(basename, ".idl")
    snap_filename = string(basename, ".snap")
    params = br_read_params(idl_filename)
    get_snap(snap_filename, params, precision)
end

"""
    get_snap(
        file_name::String,
        params   ::Dict{String,String}
        )

Reads Bifrost *.snap binary file as an Array in dimension: (mx,my,mz,nvar).
Takes the snap-filename and its paramters as arguments. Returns `snapdata`
(the data). Assumes single floating point precision by default.

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
    variables in code units.
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
end # function br_load_snapdata

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
end # function br_load_auxdata


"""
    get_var(
        xp       ::BifrostExperiment,
        snap     ::Union{<:Integer, AbstractVector{<:Integer}},
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
    snap    ::Union{<:Integer, AbstractVector{<:Integer}},
    expdir  ::String,
    variable::String,
    args...
    ;
    kwargs...
    )
    nsnaps = length(snap)
    basename, basename_isnap = get_basename(expname, snap, expdir)
    params = br_read_params(string(basename_isnap, ".idl"))

    if variable == "t"
        # The special case of getting the snapshot time
        return get_time(basename, snap, params; kwargs...)
    else
        varnr, file_ext = get_varnr_and_file_extension(params, variable)
    end
    
    # Check if we're dealing with multiple snapshots or not, then fetch the 
    # data
    if nsnaps == 1
        data = get_var(
            string(basename_isnap, file_ext), 
            params, 
            varnr, 
            args...
            )
    else
        data = get_var(
            basename, 
            snap,
            params,
            varnr, 
            file_ext, 
            args...
            )
    end
    # -------------------------------------------------------------------------
    # KEYWORD ARGUMENTS
    #
    #  Below is where you extend the functionality of get_var by handling 
    #  arbitrary keyword arguments. Please put your new implementation into a 
    #  new function.
    # -------------------------------------------------------------------------
    kwarg_keys = keys(kwargs)
    kwarg_values = values(kwargs)
    # -------------------------------------------------------------------------
    #
    # Add more kwargs here
    #
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
    varnr          ::Integer,
    precision      ::DataType=Float32,
    slicex         ::AbstractVector{<:Integer}=Int[],
    slicey         ::AbstractVector{<:Integer}=Int[],
    slicez         ::AbstractVector{<:Integer}=Int[]
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

"""
    get_var(
        basename ::String,
        snaps    ::AbstractVector{<:Integer},
        params   ::Dict{String, String},
        varnr    ::Integer,
        file_ext ::String,
        args...
        )
Load variable nr `varnr` from multiple snapshots. Requires experiment basename 
(expdir + expname) and the file extension (.aux/.snap). Slicing the snapshot
is optional. Assumes single precision snapshot by default.
"""
function get_var(
    basename ::String,
    snaps    ::AbstractVector{<:Integer},
    params   ::Dict{String, String},
    varnr    ::Integer,
    file_ext ::String,
    precision::DataType=Float32;
    slicex   ::AbstractVector{<:Integer}=Int[],
    slicey   ::AbstractVector{<:Integer}=Int[],
    slicez   ::AbstractVector{<:Integer}=Int[]
    )
    # Get spatial size
    mx, my, mz = get_snapsize(params, slicex, slicey, slicez)
    # Allocate space for variable
    var = Array{precision}(undef, mx, my, mz, length(snaps))
    # Loop over snapshots
    Threads.@threads for (i,snap) in collect(enumerate(snaps))
        # !! Create thread-local variables to avoid race conditions !!
        isnap_local = lpad(snap,3,"0")
        idl_file_local = string(basename, "_", isnap_local, ".idl")
        params_local = br_read_params(idl_file_local)        
        tmp_file = string(basename, "_", isnap_local, file_ext)

        var[:,:,:,i] .= get_var(
            tmp_file,
            params_local,
            varnr,
            precision,
            slicex,
            slicey,
            slicez
            )
        
        # Need manual call to run garbage collector within threads
        GC.safepoint()
    end
    
    return var
end

function get_time(
    filename_prefix::String,
    snap           ::Union{<:Integer, AbstractVector{<:Integer}},
    params         ::Dict{String,String} 
    ;
    kwargs...
    )
    nsnaps = length(snap)
    if nsnaps == 1
    var = parse(Float64, params["t"])
        return var
    else 
        var = Vector{Float64}(undef, nsnaps)
        file_ext = ".idl"
        # Load the variable directly from params
        for (i,snap) in enumerate(snap)
            isnap = lpad(snap,3,"0")
            idl_file = string(filename_prefix, isnap, file_ext)
            params = br_read_params(idl_file) 
            time = parse(Float64, params["t"])
            var[i] = time
        end
        if :units in keys(kwargs)
            convert_timeunits!(var, params)
        end
        return var
    end
end

"""
    function get_staggered_var(
        xp::BifrostExperiment,
        snap::Integer,
        variable::String
        ;
        slicex::AbstractVector{<:Integer}=Int[],
        slicey::AbstractVector{<:Integer}=Int[],
        slicez::AbstractVector{<:Integer}=Int[]
        kwargs...)

Function to load a staggered variable and interpolate it to cell center.
The staggered variables that typically need to be interpolated are the velocity
and magnetic field components. Normally you need to use `direction="zup"` for 
vz and bz with `periodic=false` (these are the default arguments), and 
`direction="xup"` for vx and bx with `periodic=true` (same for y direction).

`kwargs`:
    precision::DataType=Float32,
    units::String="none",
    direction::String="zup",
    periodic::Bool=false,
    order::Int=6,
    slicex::AbstractVector{<:Integer}=Int[],
    slicey::AbstractVector{<:Integer}=Int[],
    slicez::AbstractVector{<:Integer}=Int[]
"""
function get_staggered_var(
    xp::BifrostExperiment,
    snaps::Union{<:Integer, AbstractVector{<:Integer}},
    variable::String;
    slicex::AbstractVector{<:Integer}=Int[],
    slicey::AbstractVector{<:Integer}=Int[],
    slicez::AbstractVector{<:Integer}=Int[],
    kwargs...)

    if typeof(snaps) <: Integer
        return get_staggered_var(xp.expname,snaps,xp.expdir,variable;
                slicex=slicex,slicey=slicey,slicez=slicez,kwargs...)
    elseif typeof(snaps) <: AbstractVector{<:Integer}
        mx, my, mz = get_snapsize(xp.mesh, slicex, slicey, slicez)
        var = Array{Float32}(undef, mx, my, mz, length(snaps))
        Threads.@threads for (i,snap) in collect(enumerate(snaps))
            var[:,:,:,i] .= get_staggered_var(xp.expname,snap,xp.expdir,variable;
                            slicex=slicex,slicey=slicey,slicez=slicez,kwargs...)
            
        end
        return var
    end
end

"""
    function get_staggered_var(expname::String,
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
"""
function get_staggered_var(
    expname::String,
    snap::Integer,
    expdir::String,
    variable::String,
    precision::DataType=Float32
    ;
    units::String="none",
    direction::String="zup",
    periodic::Bool=false,
    order::Int=6,
    slicex::AbstractVector{<:Integer}=Int[],
    slicey::AbstractVector{<:Integer}=Int[],
    slicez::AbstractVector{<:Integer}=Int[]
    )

    shift_functions = Dict(
        "xdn" => br_xdn, 
        "xup" => br_xup,
        "ydn" => br_ydn, 
        "yup" => br_yup, 
        "zdn" => br_zdn, 
        "zup" => br_zup
    )

    allowed_directions = collect(keys(shift_functions))
    if ! ( direction in allowed_directions )
        throw(ErrorException("$direction not a valid keyword argument, use "*
        "one of the following: $((allowed_directions[1:end-1] .* ", ")...)"*
        "$(allowed_directions[end])."))
    end

    shift = shift_functions[direction]

    if ! ( variable in keys(primary_vars) )
        throw(ErrorException("Variable $variable is not a primary variable"))
    end

    # io stuff
    isnap = lpad(snap,3,"0")
    idl_file = string(expname,"_",isnap,".idl")
    params = br_read_params(joinpath(expdir,idl_file))
    filename = string(expname,"_",isnap,".snap")
    filename = joinpath(expdir,filename)

    slicing = true
    ( isempty(slicex) && isempty(slicey) && isempty(slicez) ) && ( slicing = false )
    
    if slicing
        if ( direction == "xup" ) || ( direction == "xdn" )
            if isempty(slicex)
                # All indices in 'x' are loaded, don't worry about slicing
                var = get_var(filename,params,variable,precision,
                    units=units,slicey=slicey,slicez=slicez)
                var = shift(var,periodic,order)
            else
                var = get_var(filename,params,variable,precision,
                    units="none",slicey=slicey,slicez=slicez)
                var = shift(var,slicex,periodic,order)
                if units ≠ "none"
                    convert_units!(var,variable,units)
                end
            end
        elseif ( direction == "yup" ) || ( direction == "ydn" )
            if isempty(slicey)
                # All indices in 'y' are loaded, don't worry about slicing
                var = get_var(filename,params,variable,precision,
                    units=units,slicex=slicex,slicez=slicez)
                var = shift(var,periodic,order)
            else
                var = get_var(filename,params,variable,precision,
                    units="none",slicex=slicex,slicez=slicez)
                var = shift(var,slicey,periodic,order)
                if units ≠ "none"
                    convert_units!(var,variable,units)
                end
            end
        else
            if isempty(slicez)
                # All indices in 'z' are loaded, don't worry about slicing
                var = get_var(filename,params,variable,precision,
                    units=units,slicex=slicex,slicey=slicey)
                var = shift(var,periodic,order)
            else
                var = get_var(filename,params,variable,precision,
                    units="none",slicex=slicex,slicey=slicey)
                var = shift(var,slicez,periodic,order)
                if units ≠ "none"
                    convert_units!(var,variable,units)
                end
            end
        end
    else
        # load the entire variable and shift it in the desired direction
        var = get_var(expname,snap,expdir,variable,precision,
            units=units)
        var = shift(var,periodic,order)
    end

    GC.safepoint()
    return var

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
    kwargs...)

    if typeof(snaps) <: Integer
        return get_electron_density(xp.expname,snaps,xp.expdir;
		slicex=slicex,slicey=slicey,slicez=slicez,kwargs...)
    elseif typeof(snaps) <: AbstractVector{<:Integer}
        mx, my, mz = get_dims(slicex, slicey, slicez, xp.mesh)
        var = Array{Float32}(undef, mx, my, mz, length(snaps))
        Threads.@threads for (i,snap) in collect(enumerate(snaps))
            var[:,:,:,i] .= get_electron_density(xp.expname,snap,xp.expdir;
                            slicex=slicex,slicey=slicey,slicez=slicez,kwargs...)
            
            # Need manual call to run garbage collector within threads
            GC.safepoint()
        end
        return var
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

    # rho in g/cm^3
    ( isempty(rho) ) && ( rho = get_var(expname,snap,expdir,"r",
    units="cgs",slicex=slicex,slicey=slicey,
    slicez=slicez) )
    
    # internal energy in ergs 
    ( isempty(e) ) && ( e = get_var(expname,snap,expdir,"e",
    units="cgs",slicex=slicex,slicey=slicey,
    slicez=slicez) )

    # Calculate internal energy per mass
    ee = e ./ rho
    
    # construct the EOS tables for interpolation of electron density
    tabfile = joinpath(expdir,tabfile)
    eos = EOS_tables(tabfile)

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
    itp_table = br_eos_interpolate(eos,3)

    x = log.(ee)
    y = log.(rho)

    ne = itp_table.(x, y)

    # take exp to remove log
    ne = exp.(ne)

    # Convert to si on request (cm^-3 --> m^-3)
    if units == "si"
        ne .*= 1f6
    end

    return ne
end
