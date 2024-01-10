const primary_vars = Dict(
    "r" => 1,
    "px" => 2, 
    "py" => 3, 
    "pz" => 4, 
    "e" => 5, 
    "bx" => 6, 
    "by" => 7, 
    "bz" => 8,  
)

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
    br_load_snapdata(
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
function br_load_snapdata(
    expname  ::String,
    snap     ::Int,
    expdir   ::String,
    precision::DataType=Float32
    )
    datadims = 4 # 3 spatial dimensions and 1 variable dimension
    # Parse filenames
    basename = string(expdir, "/", expname, "_$(lpad(snap,3,"0"))")
    idl_filename = string(basename, ".idl")
    snap_filename = string(basename, ".snap")
    params = br_read_params(idl_filename)
    
    snapsize, numvars, _ = get_snapsize_and_numvars(params)

    file = open(snap_filename)
    # Use Julia standard-library memory-mapping to extract file values
    snapdata = mmap(file, Array{precision, 4}, (snapsize..., numvars))
    close(file)
    return snapdata, params
end # function br_load_snapdata

"""
    br_load_snapdata(
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
function br_load_snapdata(
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
    br_load_auxdata(
        file_name::String,
        params   ::Dict{String,String}
    )
Reads Bifrost *.aux binary file using memory-mapping. The returned
`auxdata` array will have dimensions (mx,my,mz,nvars) where nvars is the
number of aux-variables. Assumes single floating point precision by default.
"""
function br_load_auxdata(
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
    get_variable_offset_in_file(
        precision::DataType,
        snapsize::Tuple{Integer, Integer, Integer},
        varnr   ::Integer
        )
Given the precision and size of a snapshot, find the offset for reading the
variable with index `varnr` directly from file. Offset given in number of bytes.
"""
function get_variable_offset_in_file(
    precision::DataType,
    snapsize::Tuple{Integer, Integer, Integer},
    varnr   ::Integer
    )
    if precision == Float32
        bytes_per_value = 4
    elseif precision == Float64
        bytes_per_value = 8
    end 
    values_per_variable = snapsize[1]*snapsize[2]*snapsize[3]
    offset = bytes_per_value*values_per_variable*(varnr - 1)
    return offset
end

"""
    get_snapvarnr(
        variable::String,
        upperlimit::Integer=8
        )
Given a primary `variable` (in string format), return its index in a Bifrost
snapshot.
"""
function get_snapvarnr(
    variable::String,
    upperlimit::Integer=8
    )
    if variable in keys(primary_vars)
        varnr = primary_vars[variable]
    else
        error("Variable name not known.")
    end
    # In case one tries to get magnetic field with do_mhd = false
    if varnr > upperlimit
        error("Variable number exceeding number of available variables.")
    end 
    return varnr
end # function find_snapvarnr

"""
    get_auxvarnr(
        params::Dict{String,String},
        auxvar::String
        )
Given the snapshot `params` and an auxiliary variable (in string format), return
its index in the ".aux"-array.
"""
function get_auxvarnr(
    params::Dict{String,String},
    auxvar::String
    )
    indices = findall(x -> x == auxvar, split(params["aux"]))
    if length(indices) > 1
        error("Multiple matches for given aux variable name.")
    elseif length(indices) == 0
        error("Auxiliary variable not found in file.")
    end
    return indices[1]
end

"""
    br_load_snapvariable(
        file_name::String,
        params   ::Dict{String,String},
        variable ::String,
        precision::DataType=Float32;
        units::String="none",
        slicex::AbstractVector{<:Integer}=Int[],
        slicey::AbstractVector{<:Integer}=Int[],
        slicez::AbstractVector{<:Integer}=Int[]
        )
Reads a single `variable` of the Bifrost snapshot `file_name`.
Assumes single floating point precision as default.
The available variables are: 
    "r":  density
    "px": x-component of momentum
    "py": y-component of momentum
    "pz": z-component of momentum
    "e":  energy

   if params["do_mhd"] == true
    "bx": x-component of magnetic field
    "by": y-component of magnetic field
    "bz": z-component of magnetic field

Can convert variable to si or cgs units by passing `units="si"` or
`units="cgs"`.
"""
function br_load_snapvariable(
    file_name::String,
    params   ::Dict{String,String},
    variable ::String,
    precision::DataType=Float32;
    units::String="none",
    slicex::AbstractVector{<:Integer}=Int[],
    slicey::AbstractVector{<:Integer}=Int[],
    slicez::AbstractVector{<:Integer}=Int[]
    )
    datadims = 3 # 3 spatial dimensions and 1 variable dimension

    snapsize, _, _ = get_snapsize_and_numvars(params)

    # Do slicing or not (returns the mmap)
    slicing = true
    ( isempty(slicex) && isempty(slicey) && isempty(slicez) ) && ( slicing = false )

    # Get variable number/position in snapdata from argument.
    varnr = get_snapvarnr(variable)

    # Calculate offset in file
    offset = get_variable_offset_in_file(precision, snapsize, varnr)

    file = open(file_name)

    # Use Julia standard-library memory-mapping to extract file values
    if slicing
        isempty(slicex) && ( slicex = 1:snapsize[1] )
        isempty(slicey) && ( slicey = 1:snapsize[2] )
        isempty(slicez) && ( slicez = 1:snapsize[3] )
        # slice the variable
        snapvariable = mmap(file,
                            Array{precision, datadims}, 
                            snapsize, 
                            offset)[slicex,slicey,slicez]
        
        if units ≠ "none"
            convert_units!(snapvariable, variable, units)
        end
    else
        # do not slice the variable
        snapvariable = mmap(file,
                            Array{precision, datadims}, 
                            snapsize, 
                            offset)
        
        if units ≠ "none"
            # Allocate the variable before changing units
            snapvariable = convert_units(snapvariable, variable, units)
        end
    end

    close(file)

    return snapvariable
end # function br_load_snapvariable

"""
    br_load_snapvariable(
        expname ::String,
        snap    ::Vector{T} where {T<:Integer},
        expdir  ::String,
        variable::String,
        precision::DataType=Float32;
        units::String="none"
        )
Reads a single primary `variable` of one or multiple Bifrost snapshots.
Takes the snapshot-numbers in the vector `snap`.
Assumes single floating point precision as default.
The available variables are: 
    "r":  density
    "px": x-component of momentum
    "py": y-component of momentum
    "pz": z-component of momentum
    "e":  energy

   if params["do_mhd"] == true
    "bx": x-component of magnetic field
    "by": y-component of magnetic field
    "bz": z-component of magnetic field

Can convert variable to si or cgs units by passing  `units="si"` or
`units="cgs"`.
"""
function br_load_snapvariable(
    expname ::String,
    snap    ::Vector{T} where {T<:Integer},
    expdir  ::String,
    variable::String,
    precision::DataType=Float32;
    units::String="none"
    )
    datadims = 3 # 3 spatial dimensions and 1 variable dimension

    # Parse filenames
    basename = string(expdir, "/", expname, "_$(lpad(snap[1],3,"0"))")
    idl_filename  = string(basename, ".idl")
    params = br_read_params(idl_filename)

    snapsize, numvars = get_snapsize_and_numvars(params)

    # Get variable number/position in snapdata from argument.
    varnr = get_snapvarnr(variable)

    # Calculate offset in file
    offset = get_variable_offset_in_file(precision, snapsize, varnr)

    # Loop through all snapshots and load variable
    numsnaps = length(snap)
    snapvariable = zeros(precision, snapsize..., numsnaps)
    for i = 1:numsnaps
        isnap = "_"*lpad(snap[i],3,"0")
        basename = string(expdir, "/", expname, isnap)
        snap_filename = string(basename, ".snap")
        file = open(snap_filename)
        # Use Julia standard-library memory-mapping to extract file values
        snapvariable[:,:,:,i] .= mmap(file,
                                     Array{precision, datadims},
                                     snapsize,
                                     offset)
        close(file)
    end

    if units ≠ "none"
        convert_units!(snapvariable, variable, units)
    end

    return snapvariable
end # function br_load_snapvariable

"""
    br_load_auxvariable(
        file_name::String,
        params   ::Dict{String,String},
        auxvar   ::String,
        precision::DataType=Float32;
        units::String="none",
        slicex::AbstractVector{<:Integer}=Int[],
        slicey::AbstractVector{<:Integer}=Int[],
        slicez::AbstractVector{<:Integer}=Int[]
        )
Reads a single auxiliary variable from a Bifrost ".aux"-file `file_name`. The
snapshot parameters must be given together with a string for the aux-variable,
`auxvar`. Assumes single floating point precision as default. Can convert 
variable to si or cgs units by passing  `units="si"` or 
`units="cgs"`.
"""
function br_load_auxvariable(
    file_name::String,
    params::Dict{String,String},
    auxvar::String,
    precision::DataType=Float32;
    units::String="none",
    slicex::AbstractVector{<:Integer}=Int[],
    slicey::AbstractVector{<:Integer}=Int[],
    slicez::AbstractVector{<:Integer}=Int[]
    )

    datadims = 3
    snapsize, _, _ = get_snapsize_and_numvars(params)

    # Do slicing or not (returns the mmap)
    slicing = true
    ( isempty(slicex) && isempty(slicey) && isempty(slicez) ) && ( slicing = false )

    auxvarnr = get_auxvarnr(params, auxvar)
    # Load file and auxvariable
    file = open(file_name)
    offset = get_variable_offset_in_file(precision, snapsize, auxvarnr)

    # Use Julia standard-library memory-mapping to extract file values
    if slicing
        isempty(slicex) && ( slicex = 1:snapsize[1] )
        isempty(slicey) && ( slicey = 1:snapsize[2] )
        isempty(slicez) && ( slicez = 1:snapsize[3] )
        # slice the variable
        auxvariable = mmap(file,
                           Array{precision, datadims}, 
                           snapsize, 
                           offset)[slicex,slicey,slicez]
        
        if units ≠ "none"
            convert_units!(auxvariable, auxvar, units)
        end
    else
        # do not slice the variable
        auxvariable = mmap(file,
                           Array{precision, datadims}, 
                           snapsize, 
                           offset)
        
        if units ≠ "none"
            # Allocate the variable to change units
            auxvariable = convert_units(auxvariable, auxvar, units)
        end
    end

    close(file)

    return auxvariable
end # function br_load_auxdata

"""
    br_load_auxvariable(
        expname ::String,
        snap    ::Vector{T} where {T<:Integer},
        expdir  ::String,
        auxvar  ::String,
        precision::DataType=Float32;
        units::String="none"
        )
Reads a single axiliary variable (`auxvar`) of one or multiple Bifrost
snapshots. Takes the snapshot-numbers in the vector `snap`.
Assumes single floating point precision as default. Can convert variable to si 
or cgs units by passing  `units="si"` or `units="cgs"`.
"""
function br_load_auxvariable(
    expname ::String,
    snap    ::Vector{T} where {T<:Integer},
    expdir  ::String,
    auxvar  ::String,
    precision::DataType=Float32;
    units::String="none"
    )
    datadims = 3

    # Parse filenames
    basename = string(expdir, "/", expname, "_$(lpad(snap[1],3,"0"))")
    idl_filename  = string(basename, ".idl")
    params = br_read_params(idl_filename)

    snapsize, _, numauxvars = get_snapsize_and_numvars(params)
    auxvarnr = get_auxvarnr(params, auxvar)
    # Calculate offset in file
    offset = get_variable_offset_in_file(precision, snapsize, auxvarnr)

    # Loop through all snapshots and load variable
    numsnaps = length(snap)
    auxvariable = zeros(precision, snapsize..., numsnaps)
    for i = 1:numsnaps
        isnap = "_"*lpad(snap[i],3,"0")
        basename = string(expdir, "/", expname, isnap)
        aux_filename = string(basename, ".aux")
        file = open(aux_filename)
        # Use Julia standard-library memory-mapping to extract file values
        auxvariable[:,:,:,i] = mmap(file,
                                    Array{precision, datadims},
                                    snapsize,
                                    offset)
        close(file)
    end

    if units ≠ "none"
        convert_units!(auxvariable,auxvar,units)
    end

    return auxvariable
end

"""
    get_var(
        xp::BifrostExperiment,
        snap::Union{<:Integer, AbstractVector{<:Integer}},
        variable::String
        )
 Loads a variable from a simulation snapshot. `snap` can either be an integer 
snapshot or an integer list of snapshots. 

Available variables

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
___
OPTIONAL KEYWORD-ARGUMENTS
Converts variables to "si" or "cgs" units: `units="si"` or `units="cgs"`.

To load a slice of the variable, give e.g. `slicex=[32, 410]` or `slicey=40:90`

Example usage: 

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
    xp::BifrostExperiment,
    snap::Union{<:Integer, AbstractVector{<:Integer}},
    variable::String;
    kwargs...
    )
    
    get_var(xp.expname, snap, xp.expdir, variable; kwargs...)
end

function get_var(
    expname::String,
    snap::Integer,
    expdir::String,
    variable::String,
    precision::DataType=Float32
    ;
    units::String="none",
    slicex::AbstractVector{<:Integer}=Int[],
    slicey::AbstractVector{<:Integer}=Int[],
    slicez::AbstractVector{<:Integer}=Int[]
    )


    isnap = lpad(snap,3,"0")
    idl_file = string(expname,"_",isnap,".idl")
    params = br_read_params(joinpath(expdir,idl_file))

    if variable in keys(primary_vars)
        filename = string(expname,"_",isnap,".snap")
        filename = joinpath(expdir,filename)
        return br_load_snapvariable(filename,params,variable,precision,
            units=units,slicex=slicex,slicey=slicey,slicez=slicez)
            
    elseif variable in split(params["aux"])
        filename = string(expname,"_",isnap,".aux")
        filename = joinpath(expdir,filename)
        return br_load_auxvariable(filename,params,variable,precision,
            units=units,slicex=slicex,slicey=slicey,slicez=slicez)
    
    elseif variable == "t"
        var = parse(Float64, params["t"])
        if units == "si" || units == "cgs"
            var = convert_snaptime(var)
        end
    
        return var

    else
        throw(ErrorException("Variable $variable does not exist"))
    end

end

function get_var(
    expname::String,
    snaps::AbstractVector{<:Integer},
    expdir::String,
    variable::String,
    precision::DataType=Float32
    ;
    units::String="none",
    slicex::AbstractVector{<:Integer}=Int[],
    slicey::AbstractVector{<:Integer}=Int[],
    slicez::AbstractVector{<:Integer}=Int[]
    )

    isnap = lpad(snaps[1],3,"0")
    idl_file = string(expname,"_",isnap,".idl")
    params = br_read_params(joinpath(expdir,idl_file))

    filename = joinpath(expdir,expname*"_")
    
    if variable in keys(primary_vars)
        load_var = br_load_snapvariable
        file_ext = ".snap"
    elseif variable in split(params["aux"])
        load_var = br_load_auxvariable
        file_ext = ".aux"
    elseif variable == "t"
        var = Vector{Float64}(undef, length(snaps))
        file_ext = ".idl"
        
        # Load the variable directly from params
        for (i,snap) in collect(enumerate(snaps))
            
            isnap = lpad(snap,3,"0")
            idl_file = string(filename,isnap,file_ext)
            params = br_read_params(idl_file) 
            
            time = parse(Float64, params["t"])
            
            if ( units == "si" ) || ( units == "cgs" )
                time = convert_snaptime(time)
            end
            var[i] = time
        end

        return var
    else
        throw(ErrorException("Variable $variable does not exist"))
    end
    
    # Get spatial size
    mx, my, mz = get_dims(slicex, slicey, slicez, params)
    
    # Allocate space for variable
    var = Array{precision}(undef, mx, my, mz, length(snaps))
    
    # Loop over snapshots
    Threads.@threads for (i,snap) in collect(enumerate(snaps))
        # !! Create thread-local variables to avoid race conditions !!
        isnap_local = lpad(snap,3,"0")
        idl_file_local = string(filename,isnap_local,".idl")
        params_local = br_read_params(idl_file_local)        
        tmp_file = string(filename,isnap_local,file_ext)

        var[:,:,:,i] .= load_var(tmp_file,params_local,variable,
            units=units,slicex=slicex,slicey=slicey,slicez=slicez)
        
        # Need manual call to run garbage collector within threads
        GC.safepoint()
    end
    
    return var
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
        mx, my, mz = get_dims(slicex, slicey, slicez, xp.mesh)
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
                var = br_load_snapvariable(filename,params,variable,precision,
                    units=units,slicey=slicey,slicez=slicez)
                var = shift(var,periodic,order)
            else
                var = br_load_snapvariable(filename,params,variable,precision,
                    units="none",slicey=slicey,slicez=slicez)
                var = shift(var,slicex,periodic,order)
                if units ≠ "none"
                    convert_units!(var,variable,units)
                end
            end
        elseif ( direction == "yup" ) || ( direction == "ydn" )
            if isempty(slicey)
                # All indices in 'y' are loaded, don't worry about slicing
                var = br_load_snapvariable(filename,params,variable,precision,
                    units=units,slicex=slicex,slicez=slicez)
                var = shift(var,periodic,order)
            else
                var = br_load_snapvariable(filename,params,variable,precision,
                    units="none",slicex=slicex,slicez=slicez)
                var = shift(var,slicey,periodic,order)
                if units ≠ "none"
                    convert_units!(var,variable,units)
                end
            end
        else
            if isempty(slicez)
                # All indices in 'z' are loaded, don't worry about slicing
                var = br_load_snapvariable(filename,params,variable,precision,
                    units=units,slicex=slicex,slicey=slicey)
                var = shift(var,periodic,order)
            else
                var = br_load_snapvariable(filename,params,variable,precision,
                    units="none",slicex=slicex,slicey=slicey)
                var = shift(var,slicez,periodic,order)
                if units ≠ "none"
                    convert_units!(var,variable,units)
                end
            end
        end
    else
        # load the entire variable and shift it in the desired direction
        var = br_load_snapvariable(filename,params,variable,precision,
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
