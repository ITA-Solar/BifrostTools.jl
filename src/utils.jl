"""
    get_snapsize_and_numvars(
        params::Dict{String,String},
    )
Returns snapsize (mx, my, mz), number of primary variables and number of
auxiliary variables, given the snapshot-parameters.
"""
function get_snapsize_and_numvars(
    params::Dict{String,String},
    )
    return get_snapsize(params), get_numvars(params)
end


"""
    get_snapsize(
        params::Dict{String,String},
    )
Returns snapsize (mx, my, mz) given the snapshot-parameters.
"""
function get_snapsize(
    params::Dict{String,String},
    )
    mx = parse(Int, params["mx"])
    my = parse(Int, params["my"])
    mz = parse(Int, params["mz"])
    snapsize::Tuple{Int64, Int64, Int64} = mx, my, mz
    return snapsize 
end
"""
    get_snapsize(
        params::Dict{String,String},
        slicex::AbstractVector{<:Integer},
        slicey::AbstractVector{<:Integer},
        slicez::AbstractVector{<:Integer}
    )
Returns snapsize (mx, my, mz) given the snapshot-parameters.
"""
function get_snapsize(
    params::Dict{String,String},
    slicex::AbstractVector{<:Integer},
    slicey::AbstractVector{<:Integer},
    slicez::AbstractVector{<:Integer}
    )
    if isempty(slicex)
        mx = parse(Int, params["mx"])
    else
        mx = length(slicex)
    end
    if isempty(slicey)
        my = parse(Int, params["my"])
    else
        my = length(slicey)
    end
    if isempty(slicez)
        mz = parse(Int, params["mz"])
    else
        mz = length(slicez)
    end
    return mx, my, mz
end
"""
    get_snapsize(
        mesh::BifrostMesh,
    )
Returns snapsize (mx, my, mz) given a Bifrost-mesh.
"""
function get_snapsize(
    mesh::BifrostMesh,
    )
    return (mesh.mx, mesh.my, mesh.mz)
end
"""
    get_snapsize(
        mesh::BifrostMesh,
        slicex::AbstractVector{<:Integer},
        slicey::AbstractVector{<:Integer},
        slicez::AbstractVector{<:Integer}
    )
Returns snapsize (mx, my, mz) given a Bifrost-mesh.
"""
function get_snapsize(
    mesh::BifrostMesh,
    slicex::AbstractVector{<:Integer},
    slicey::AbstractVector{<:Integer},
    slicez::AbstractVector{<:Integer}
    )    
    if isempty(slicex)
        mx = mesh.mx
    else
        mx = length(slicex)
    end
    if isempty(slicey)
        my = mesh.my
    else
        my = length(slicey)
    end
    if isempty(slicez)
        mz = mesh.mz
    else
        mz = length(slicez)
    end
    return mx, my, mz
end


"""
    get_numvars(
        params::Dict{String,String},
    )
Returns number of primary variables and number of
auxiliary variables, given the snapshot-parameters.
"""
function get_numvars(
    params::Dict{String,String},
    )
    if parse(Int, params["do_mhd"]) == 1
        numvars = 8
    else
        numvars = 5
    end
    numauxvars::Int64 = length(split(params["aux"]))
    return numvars, numauxvars
end

"""
    get_snap_numbers(
        expdir::String, 
        expname::String="none"
        ;
        findall::Bool=false, 
        filenames::Vector{String}=String[]
        )

Finds all files in the format 'expname_XXX.snap' in the experiment directory
`exp_dir`, and returns a vector of the snapshots XXX. If `expname` is not
given, is is assumed that the directory of the simulation matches the
experiment name.
"""
function get_snap_numbers(
    expdir::String, 
    expname::String="none"
    ;
    findall::Bool=false, 
    filenames::Vector{String}=String[]
    )

    if expname=="none"
        expname = splitpath(expdir)[end]
    end 

    if isempty(filenames)
        filenames = readdir(expdir)
    end

    if ! findall
        # Regex magic to match the filenames with 'expname' and '.snap'
        pattern = r"^" * expname * r"_(\d+)\.snap$"
    else
        # wildcard that finds all files on format 'abXYcd_xyz.snap' 
        pattern = r"^.*_(\d+)\.snap$"
    end

    # Initialize an empty list to store the XXX numbers
    snaps = Vector{Int}()

    # Loop through the filenames and extract XXX numbers
    for filename in filenames
        match_result = match(pattern, filename)
        if match_result ≠ nothing
            isnap = Meta.parse(match_result.captures[1])
            push!(snaps, isnap)
        end
    end

    return sort(snaps)

end


"""
    get_varnr_and_file_ext(
        params::Dict{String,String},
        variable::String
        )
Given the snapshot `params` and desired `variable`, return
its index in the binary file, as well as the extension of this file.
(Either ".aux" or ".snap").
"""
function get_varnr_and_file_extension(
    params  ::Dict{String,String},
    variable::String,
    )
    if variable in keys(primary_vars)
        file_ext = ".snap"
        varnr = primary_vars[variable]
    elseif variable in split(params["aux"])
        file_ext = ".aux"
        indices = findall(x -> x == variable, split(params["aux"]))
        if length(indices) > 1  
            error("Multiple matches for given aux-variable name.")
        elseif length(indices) == 0
            throw(ErrorException("Auxiliary variable not found in file."))
        end
        varnr =  indices[1]
    else
        throw(ErrorException("Variable $variable does not exist"))
    end
    return varnr, file_ext
end



"""
    get_basename(
        expname ::String,
        snap    ::Union{<:Integer, AbstractVector{<:Integer}},
        expdir  ::String,
        )
Return the basename of snapshots in the experiment `expname`, located in the 
directory `expdir`. Also return the filename (withou file extension) of the 
first snapshot of the experiment.
"""
function get_basename(
    expname ::String,
    snap    ::Union{<:Integer, AbstractVector{<:Integer}},
    expdir  ::String,
    )
    isnap = lpad(snap[1],3,"0")
    basename =  joinpath(expdir, expname)
    return  basename, string(basename, "_", isnap)
end
