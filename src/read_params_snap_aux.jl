const destaggeroperation = Dict(
    "px" => xup,
    "py" => yup,
    "pz" => zup,
    "bx" => xup,
    "by" => yup,
    "bz" => zup,
    "ex" => yupzup,
    "ey" => zupxup,
    "ez" => xupyup,
    "ix" => yupzup,
    "iy" => zupxup,
    "iz" => xupyup
    )


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
    read_params(xp::BifrostExperiment, snap::Integer)
"""
function read_params(
    xp::BifrostExperiment,
    snap::Integer,
    )
    read_params(xp.expname, snap, xp.expdir)
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
    return get_snap(snap_filename, params, precision)
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

    snapsize, (numvars, _) = get_snapsize_and_numvars(params)

    file = open(file_name)
    # Use Julia standard-library memory-mapping to extract file values
    snapdata = mmap(file, Array{precision, datadims}, (snapsize..., numvars))
    close(file)
    return snapdata
end

"""
    get_snap(
        xp::BifrostExperiment,
        snap::Int
    )
"""
function get_snap(
    xp::BifrostExperiment,
    snap::Int
)
    get_snap(xp.expname, snap, xp.expdir)
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
    snapsize, (_, numauxvars) = get_snapsize_and_numvars(params)
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


function get_aux(
    expname  ::String,
    snap     ::Int,
    expdir   ::String,
    precision::DataType=Float32
    )
    # Parse filenames
    basename = string(joinpath(expdir, expname),"_$(lpad(snap,3,"0"))")
    idl_filename = string(basename, ".idl")
    aux_filename = string(basename, ".aux")
    params = read_params(idl_filename)
    return get_aux(aux_filename, params, precision)
end


function get_aux(
    xp::BifrostExperiment,
    snap::Int
)
    get_aux(xp.expname, snap, xp.expdir)
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
    precision::DataType=Float32,
    kwargs...
    )

    if variable == "t"
        # The special case of getting the snapshot time
        return get_time(expname,snaps,expdir;kwargs...)
    end

    # Allocate space for variable
    data = Vector{Array{precision, 3}}(undef, length(snaps))

    # Check if user wants data to be destaggered. If so we have to
    # call get_and_destagger_var. If not, we may call get_var
    if get(kwargs, :destagger, false)
        # Check if destagger-operation is passed as a keyword-argument.
        # If not, use default operation corresponding to the requested
        # variable.

        # get correct boundary condition if it is not given
        # only for b-field and velocity vector
        if !haskey(kwargs,:periodic)
            if variable in ["bx","by","px","py"]
                kwargs = addtokwargs(
                    ;periodic=true,
                    kwargs...
                )
            elseif variable in ["bz","pz"]
                kwargs = addtokwargs(
                    ;periodic=false,
                    kwargs...
                )
            end
        end

        if haskey(kwargs, :destaggeroperation)
            get_function = get_and_destagger_var
        elseif variable in keys(destaggeroperation)
            get_function = get_and_destagger_var
            kwargs = addtokwargs(
                ;destaggeroperation=destaggeroperation[variable],
                kwargs...
                )
        else
            error("Destaggering of $variable is not implemented. "*
                "Set the keyword-argument `destaggeroperation`"
                )
        end
    else
        get_function = get_var
    end

    # What we need in params should be constant for experiment
    params = read_params(expname,snaps[1],expdir)
    varnr, file_suff = get_varnr_and_file_suffix(params, variable)

    # Loop over snapshots
    for (i,snap) in enumerate(snaps)

        tmp_file = string(joinpath(expdir,expname),
                          Printf.format(file_suff, lpad(snap,3,"0")))

        data[i] = get_function(
            tmp_file,
            params,
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
    #   If multiple snapshots: Assumes the same conversion factor for all
    if haskey(kwargs,:units)
        data = convert_units(data, variable, params, kwargs[:units])
    end

    # ORIENTATION: Rotate coordinate system
    if haskey(kwargs,:rotate_about)
        data = rotate(data, variable, kwargs[:rotate_about])
    end

    # SQUEEZE: Drop empty dimensions
    #   Allocates new data with fewer dims, and copies this into data
    if get(kwargs,:squeeze,false)
        # number of non-empty dimension
        dims = count( size(data[1]) .≠ 1 )

        if dims == 0
            new_data = Vector{precision}(undef,length(snaps))
        else
            new_data = Vector{Array{precision,dims}}(undef,length(snaps))
        end
        for i in eachindex(data)
            new_data[i] = squeeze(data[i])
        end
        data = new_data
    end

    # CONCATENATION: Concatenate vector to 3D array if single snapshot
    if length(snaps) == 1
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
    # Calculate offset in file
    offset = get_variable_offset_in_file(precision, snapsize, varnr)
    file = open(filename)
    
    # Do slicing or not (returns the mmap)
    # Use Julia standard-library memory-mapping to extract file values
    if isempty(slicex) && isempty(slicey) && isempty(slicez)
        # do not slice the variable
        data = mmap(file,
            Array{precision, datadims},
            snapsize,
            offset
            )
    else
        isempty(slicex) && ( slicex = 1:snapsize[1] )
        isempty(slicey) && ( slicey = 1:snapsize[2] )
        isempty(slicez) && ( slicez = 1:snapsize[3] )
        # slice the variable
        data = mmap(file,
            Array{precision, datadims},
            snapsize,
            offset
            )[slicex,slicey,slicez]
    end
    close(file)
    return data
end

function get_time(
    expname::String,
    snap   ::Union{<:Integer, AbstractVector{<:Integer}},
    expdir ::String
    ;
    units::String="code",
    kwargs...
    )
    nsnaps = length(snap)
    if typeof(snap) <: Integer
        params = read_params(expname,snap,expdir)
        data = parse(Float64, params["t"])
    else
        data = Vector{Float64}(undef, nsnaps)
        # Load the variable directly from params
        for (i,snap) in enumerate(snap)
            params = read_params(expname,snap,expdir)
            time = parse(Float64, params["t"])
            data[i] = time
        end
    end

    if units != "code"
        data = convert_timeunits(data, params)
    end

    return data
end

"""
    function get_and_destagger_var(expname::String,
        filename::String,
        params::Dict{String,String},
        varnr::Integer,
        ;
        destaggeroperation::Function,
        units::String="none",
        periodic::Bool=false,
        order::Int=6,
        slicex::AbstractVector{<:Integer}=Int[],
        slicey::AbstractVector{<:Integer}=Int[],
        slicez::AbstractVector{<:Integer}=Int[],
        kwargs...
        )

Function to load a staggered variable and interpolate it to cell center.
The staggered variables that typically need to be interpolated are the velocity
and magnetic field components. Normally you need to use `destaggeroperation=zup`
for vz and bz with `periodic=false`, and `destaggeroperation=xup` for vx and bx
with `periodic=true` (same for y direction).
"""
function get_and_destagger_var(
    filename::String,
    params::Dict{String,String},
    varnr::Integer,
    ;
    destaggeroperation::Function,
    periodic::Bool=false,
    order::Int=6,
    slicex::AbstractVector{<:Integer}=Int[],
    slicey::AbstractVector{<:Integer}=Int[],
    slicez::AbstractVector{<:Integer}=Int[],
    kwargs...
    )
    if isempty(slicex) && isempty(slicey) && isempty(slicez)
        # load the entire variable and destaggeroperation it in the desired direction
        data = get_var(filename,params,varnr,)
        data = destaggeroperation(data,periodic,order)
    else
        if destaggeroperation in (xup, xdn)
            # Load var
            data = get_var(filename,params,varnr;
                slicey=slicey,slicez=slicez,kwargs...
                )
            if isempty(slicex)
                # All indices in 'x' are loaded, don't worry about slicing
                data = destaggeroperation(data,periodic,order)
            else
                # Call to the function that slices in x
                data = destaggeroperation(data,slicex,periodic,order)
            end
        elseif destaggeroperation in (yup, ydn)
            data = get_var(filename,params,varnr;
                slicex=slicex,slicez=slicez,kwargs...
                )
            if isempty(slicey)
                # All indices in 'y' are loaded, don't worry about slicing
                data = destaggeroperation(data,periodic,order)
            else
                # Call to the function that slices in y
                data = destaggeroperation(data,slicey,periodic,order)
            end
        elseif destaggeroperation in (zup, zdn)
            data = get_var(filename, params, varnr;
                slicex=slicex, slicey=slicey, kwargs...
                )
            if isempty(slicez)
                # All indices in 'z' are loaded, don't worry about slicing
                data = destaggeroperation(data,periodic,order)
            else
                # Call to the function that slices in z
                data = destaggeroperation(data,slicez,periodic,order)
            end
        #
        # POSSIBLE TO SIMPLIFY THIS?
        # Always passing slice to the operation and handling it there?
        #
        elseif destaggeroperation == yupzup
            data = get_var(filename, params, varnr; slicex=slicex, kwargs...)
            if isempty(slicez) && isempty(slicey)
                # All indices in 'z' and 'y' are loaded, don't worry about slicing
                data = destaggeroperation(data,periodic,order)
            elseif isempty(slicey)
                # Call to the function that slices in z
                data = zup(data, slicez, periodic, order)
                data = yup(data, periodic, order)
            elseif isempty(slicez)
                # All indices in 'z' are loaded, don't worry about slicing
                data = zup(data, periodic, order)
                data = yup(data, slicey, periodic, order)
            else
                # Call to the function that slices in z and y
                data = zup(data, slicez, periodic, order)
                data = yup(data, slicey, periodic, order)
            end
        elseif destaggeroperation == zupxup
            data = get_var(filename, params, varnr; slicey=slicey, kwargs...)
            if isempty(slicez) && isempty(slicex)
                # All indices in 'z' and 'x' are loaded, don't worry about slicing
                data = destaggeroperation(data,periodic,order)
            elseif isempty(slicez)
                # Call to the function that slices in x
                data = xup(data, slicex, periodic, order)
                data = zup(data, periodic, order)
            elseif isempty(slicex)
                # All indices in 'y' are loaded, don't worry about slicing
                data = xup(data, periodic, order)
                data = zup(data, slicez, periodic, order)
            else
                # Call to the function that slices in x and z
                data = xup(data, slicex, periodic, order)
                data = zup(data, slicez, periodic, order)
            end
        elseif destaggeroperation == xupyup
            data = get_var(filename, params, varnr; slicez=slicez, kwargs...)
            if isempty(slicex) && isempty(slicey)
                # All indices in 'x' and 'y' are loaded, don't worry about slicing
                data = destaggeroperation(data,periodic,order)
            elseif isempty(slicex)
                # Call to the function that slices in y
                data = yup(data, slicey, periodic, order)
                data = xup(data, periodic, order)
            elseif isempty(slicey)
                # All indices in 'y' are loaded, don't worry about slicing
                data = yup(data, periodic, order)
                data = xup(data, slicex, periodic, order)
            else
                # Call to the function that slices in y and x.
                data = yup(data, slicey, periodic, order)
                data = xup(data, slicex, periodic, order)
            end
        else
            error("Destaggering of varnnr $varnr is not implemented. "*
                "Set the keyword-argument `destaggeroperation`"
                )
        end

    end
    return data
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

function rotate(
    data         ::AbstractVector,
    variable     ::String,
    rotation_axis::String,
    )

    return [rotate(data[i], variable, rotation_axis) for i in eachindex(data)]
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
    kwargs...)

    if typeof(snaps) <: Integer
        var = get_electron_density(xp.expname,snaps,xp.expdir;
		    kwargs...)
        return var

    elseif typeof(snaps) <: AbstractVector{<:Integer}
        var = Vector{Array{Float32,3}}(undef, length(snaps))
        for (i,snap) in enumerate(snaps)
            var[i] = get_electron_density(xp.expname,snap,xp.expdir;
                        kwargs...)
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
    expdir::String
    ;
    units::String="si",
    verbose::Bool=true,
    kwargs...
    )

    params = read_params(expname,snap,expdir)

    # If simulation is Hion, it has electron density in .hion file
    # If simulation is not Hion, we interpolate electron density from EOS
    try
        if parse(Int, params["do_hion"]) == 1
            if verbose
                @info "--- Snap $snap: reading electron density from hion file"
            end
            # Precision is hard-coded to be Float32
            return get_var(expname,snap,expdir,"ne";
                units=units,precision=Float32,kwargs...)
        else
            throw(KeyError("Condition not met"))
        end
    catch err # ´e´ is the internal energy variable, use ´err´ instead
        if isa(err, KeyError)
            if verbose
                @info "--- Snap $snap: interpolating electron density from EOS tables"
            end
            return interpolate_electron_density(expname,snap,expdir,params;
                units=units,kwargs...)
        else
            rethrow(err)
        end
    end

end
