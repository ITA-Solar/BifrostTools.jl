"""
    get_snapsize_and_numvars(
        params::Dict{String,Any},
    )
Returns snapsize (nx, ny, nz), number of primary variables and number of
auxiliary variables, given the snapshot-parameters.
"""
function get_snapsize_and_numvars(
    params::Dict{String,Any},
    )
    if params["do_mhd"] == 1
        numvars = 8
    else
        numvars = 5
    end
    numauxvars::Int64 = length(split(params["aux"]))
    snapsize::Tuple{Int64, Int64, Int64} = params["mx"], params["my"], params["mz"]
    return snapsize, numvars, numauxvars
end


"""
    br_find_snap_numbers

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


function get_dims(slicex::AbstractVector{<:Integer},
                  slicey::AbstractVector{<:Integer},
                  slicez::AbstractVector{<:Integer},
                  params::Dict{String, Any})

    if isempty(slicex)
        mx = params["mx"]
    else
        mx = length(slicex)
    end

    if isempty(slicey)
        my = params["my"]
    else
        my = length(slicey)
    end

    if isempty(slicez)
        mz = params["mz"]
    else
        mz = length(slicez)
    end

    return mx, my, mz

end

function get_dims(slicex::AbstractVector{<:Integer},
    slicey::AbstractVector{<:Integer},
    slicez::AbstractVector{<:Integer},
    mesh::BifrostMesh)

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
