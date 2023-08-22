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