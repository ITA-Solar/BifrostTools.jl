# --- generic stagger operations for vectors

function br_up(vec::Vector{T}) where {T<:AbstractFloat}
    n = size(vec)
    if (n[1] == 1)
        return vec
    else
        out = zeros(T, (n[1]+10))
        out[5:n[1]+4] .= vec[:]
        for i = 1:4
            @inbounds out[5-i]            = 2.0 .* out[5] .- out[5+i]
            @inbounds out[n[1]+4+i] = 2.0 .* out[n[1]+4] .- out[n[1]+4-i]
        end

        c =     3.0/256.0
        b = -25.0/256.0
        a =     0.5-b-c

        out = (a .* (circshift(out,-1) .+ circshift(out,0)) .+
                     b .* (circshift(out,-2) .+ circshift(out,1)) .+
                     c .* (circshift(out,-3) .+ circshift(out,2))    )

        return out[5:n[1]+4]
    end
end

function br_dup(vec::Vector{T}) where {T<:AbstractFloat}
    n = size(vec)
    if (n[1] == 1)
        return vec
    else
        out = zeros(T, (n[1]+10))
        out[5:n[1]+4] .= vec[:]
        for i = 1:4
            @inbounds out[5-i]            = 2.0 .* out[5] .- out[5+i]
            @inbounds out[n[1]+4+i] = 2.0 .* out[n[1]+4] .- out[n[1]+4-i]
        end

        d =     -75.0 / 107520.0
        c =    1029.0 / 107520.0
        b = -8575.0 / 107520.0
        a = 1. - 3 * b - 5 * c - 7 * d

        out = (a .* (circshift(out,-1) .- circshift(out,0)) .+
                     b .* (circshift(out,-2) .- circshift(out,1)) .+
                     c .* (circshift(out,-3) .- circshift(out,2)) .+
                     d .* (circshift(out,-4) .- circshift(out,3))    )

        return out[5:n[1]+4]
    end
end

function br_dn(vec::Vector{T}) where {T<:AbstractFloat}
    n = size(vec)
    if (n[1] == 1)
        return vec
    else
        out = zeros(T, (n[1]+10))
        out[5:n[1]+4] .= vec[:]
        for i = 1:4
            @inbounds out[5-i]            = 2.0 .* out[5] .- out[5+i]
            @inbounds out[n[1]+4+i] = 2.0 .* out[n[1]+4] .- out[n[1]+4-i]
        end

        c =     3.0/256.0
        b = -25.0/256.0
        a =     0.5-b-c

        out = (a .* (circshift(out, 0) .+ circshift(out,1)) .+
                     b .* (circshift(out,-1) .+ circshift(out,2)) .+
                     c .* (circshift(out,-2) .+ circshift(out,3))    )

        return out[5:n[1]+4]
    end
end

function br_ddn(vec::Vector{T}) where {T<:AbstractFloat}
    n = size(vec)
    if (n[1] == 1)
        return vec
    else
        out = zeros(T, (n[1]+10))
        out[5:n[1]+4] .= vec[:]
        for i = 1:4
            @inbounds out[5-i]            = 2.0 .* out[5] .- out[5+i]
            @inbounds out[n[1]+4+i] = 2.0 .* out[n[1]+4] .- out[n[1]+4-i]
        end

        d =     -75.0 / 107520.0
        c =    1029.0 / 107520.0
        b = -8575.0 / 107520.0
        a = 1. - 3 * b - 5 * c - 7 * d

        out = (a .* (circshift(out, 0) .- circshift(out,1)) .+
                     b .* (circshift(out,-1) .- circshift(out,2)) .+
                     c .* (circshift(out,-2) .- circshift(out,3)) .+
                     d .* (circshift(out,-3) .- circshift(out,4))    )

        return out[5:n[1]+4]
    end
end

# --- fast stagger operations

function br_cdxdn(arr::Array{T,3}) where {T<:AbstractFloat}
    return arr .- circshift(arr,[1,0,0])
end

function br_cdydn(arr::Array{T,3}) where {T<:AbstractFloat}
    return arr .- circshift(arr,[0,1,0])
end

function br_cdzdn(arr::Array{T,3}) where {T<:AbstractFloat}
    return arr .- circshift(arr,[0,0,1])
end

function br_cdxup(arr::Array{T,3}) where {T<:AbstractFloat}
    return -(arr .- circshift(arr,[-1,0,0]))
end

function br_cdyup(arr::Array{T,3}) where {T<:AbstractFloat}
    return -(arr .- circshift(arr,[0,-1,0]))
end

function br_cdzup(arr::Array{T,3}) where {T<:AbstractFloat}
    return -(arr .- circshift(arr,[0,0,-1]))
end

function br_cdivup(varr::Vector{Array{T,3}}) where {T<:AbstractFloat}
    return br_cdxup(varr[1]) .+ br_cdyup(varr[2]) .+ br_cdzup(varr[3])
end

function br_cdivdn(varr::Vector{Array{T,3}}) where {T<:AbstractFloat}
    return br_cdxdn(varr[1]) .+ br_cdydn(varr[2]) .+ br_cdzdn(varr[3])
end

function br_cgrad_dn(arr::Array{T,3}) where {T<:AbstractFloat}
    return [br_cdxdn(arr),br_cdydn(arr),br_cdzdn(arr)]
end

function br_cgrad_up(arr::Array{T,3}) where {T<:AbstractFloat}
    return [br_cdxup(arr),br_cdyup(arr),br_cdzup(arr)]
end

function br_claplace_du(arr::Array{T,3}) where {T<:AbstractFloat}
    return br_cdivdn(br_cgrad_up(arr))
end

function br_claplace_ud(arr::Array{T,3}) where {T<:AbstractFloat}
    return br_cdivup(br_cgrad_dn(arr))
end

# --- stagger operations

"""
    function br_xup(
        arr::Array{T,3},
        periodic::Bool=true,
        order::Int=6
    ) where {T<:AbstractFloat}

Stagger operation on `arr` by a 5th order polynomial interpolation, 
shifting the variable half a grid point upwards in the x-direction
"""
function br_xup(
    arr::Array{T,3},
    periodic::Bool=true,
    order::Int=6
) where {T<:AbstractFloat}

    n = size(arr)
    @assert length(n) == 3
    if (n[1] == 1)
        return arr[:, :, :]
    else
        if order == 6

            c = 3.0 / 256.0
            b = -25.0 / 256.0
            a = 0.5 - b - c

            # inflate matrix in x direction

            out = zeros(T, (n[1] + 5, n[2], n[3]))
            tmp = zeros(T, (n[1] + 5, n[2], n[3]))
            tmp[3:n[1]+2, :, :] .= arr[:, :, :]

            if periodic
                tmp[1:2, :, :] .= arr[end-1:end, :, :]
                tmp[n[1]+3:end, :, :] .= arr[1:3, :, :]
            else
                # extrapolate bottom
                for i = 1:2
                    @inbounds tmp[3-i, :, :] = 2.0 .* tmp[3, :, :] .- tmp[3+i, :, :]
                end
                # extrapolate top
                for i = 1:3
                    @inbounds tmp[n[1]+2+i, :, :] = 2.0 .* tmp[n[1]+2, :, :] .- tmp[n[1]+2-i, :, :]
                end
            end

            Threads.@threads for k = 1:n[3]
                for j = 1:n[2]
                    @simd for i = 3:n[1]+2
                        @inbounds out[i, j, k] =
                            a * (tmp[i, j, k] + tmp[i+1, j, k]) +
                            b * (tmp[i-1, j, k] + tmp[i+2, j, k]) +
                            c * (tmp[i-2, j, k] + tmp[i+3, j, k])
                    end
                end
            end
            return out[3:end-3, :, :]
        else # order 2
            out = zeros(T, (n[1], n[2], n[3]))
            Threads.@threads for k = 1:n[3]
                for j = 1:n[2]
                    @simd for i = 1:n[1]-1
                        @inbounds out[i, j, k] = 0.5f0 * (arr[i, j, k] + arr[i+1, j, k])
                    end
                end
            end
            if periodic
                out[:, :, end] .= 0.5f0 .* (arr[:, :, 1] .+ arr[:, :, end])
            else
                out[:, :, end] .= arr[:, :, end]
            end
            return out
        end
    end
end

"""
	br_dzup(
		arr::Array{T,3},
		dz::Vector{T}, 
		periodic::Bool=false, 
		order::Int=6
		)

Computes the spatial derivative in the x-direction of every entry in `arr` 
shifted a half grid point upwards. Defaults to the 6th order accurate Bifrost 
derivative with `order=6`, optional 2nd order accurate derivative with keyword 
`order=2`
"""
function br_dxup(
    arr::Array{T,3},
    dx::Vector{T},
    periodic::Bool=true,
    order::Int=6
) where {T<:AbstractFloat}

    n = size(arr)
    @assert length(n) == 3
    if (n[1] == 1)
        return zeros(T, (n[1], n[2], n[3]))
    else
        if order == 6

            c = (-1.0 + (3.0^5 - 3.0) / (3.0^3 - 3.0)) / (5.0^5 - 5.0 - 5.0 * (3.0^5 - 3))
            b = (-1.0 - 120.0 * c) / 24.0
            a = (1.0 - 3.0 * b - 5.0 * c)

            # inflate matrix in x direction

            out = zeros(T, (n[1] + 5, n[2], n[3]))
            tmp = zeros(T, (n[1] + 5, n[2], n[3]))
            tmp[3:n[1]+2, :, :] .= arr[:, :, :]

            if periodic
                tmp[1:2, :, :] .= arr[end-1:end, :, :]
                tmp[n[1]+3:end, :, :] .= arr[1:3, :, :]
            else
                # extrapolate bottom
                for i = 1:2
                    @inbounds tmp[3-i, :, :] = 2.0 .* tmp[3, :, :] .- tmp[3+i, :, :]
                end
                # extrapolate top
                for i = 1:3
                    @inbounds tmp[n[1]+2+i, :, :] = 2.0 .* tmp[n[1]+2, :, :] .- tmp[n[1]+2-i, :, :]
                end
            end

            Threads.@threads for k = 1:n[3]
                for j = 1:n[2]
                    @simd for i in 3:n[1]+2
                        @inbounds out[i, j, k] =
                            dx[i-2] * (
                                a * (tmp[i+1, j, k] - tmp[i, j, k]) +
                                b * (tmp[i+2, j, k] - tmp[i-1, j, k]) +
                                c * (tmp[i+3, j, k] - tmp[i-2, j, k])
                            )
                    end
                end
            end
            return out[3:end-3, :, :]
        else # oder 2
            out = zeros(T, (n[1], n[2], n[3]))
            Threads.@threads for k = 1:n[3]
                for j = 1:n[2]
                    @simd for i = 1:n[1]-1
                        @inbounds out[i, j, k] = dx[i] * (arr[i+1, j, k] - arr[i, j, k])
                    end
                end
            end
            if periodic
                out[:, :, end] .= dx[end] .* (arr[:, :, 1] .- arr[:, :, end])
            end
            return out
        end
    end
end

"""
    function br_xdn(
        arr::Array{T,3},
        periodic::Bool=true,
        order::Int=6
    ) where {T<:AbstractFloat}

Stagger operation on `arr` by a 5th order polynomial interpolation, 
shifting the variable half a grid point downwards in the x-direction
"""
function br_xdn(
    arr::Array{T,3},
    periodic::Bool=true,
    order::Int=6
) where {T<:AbstractFloat}

    n = size(arr)
    @assert length(n) == 3
    if (n[1] == 1)
        return arr[:, :, :]
    else

        if order == 6
            c = 3.0 / 256.0
            b = -25.0 / 256.0
            a = 0.5 - b - c

            # inflate matrix in x direction

            out = zeros(T, (n[1] + 5, n[2], n[3]))
            tmp = zeros(T, (n[1] + 5, n[2], n[3]))
            tmp[4:n[1]+3, :, :] .= arr[:, :, :]

            if periodic
                tmp[1:3, :, :] .= arr[end-2:end, :, :]
                tmp[n[1]+4:end, :, :] .= arr[1:2, :, :]
            else
                # extrapolate bottom
                for i = 1:3
                    @inbounds tmp[4-i, :, :] = 2.0 .* tmp[4, :, :] .- tmp[4+i, :, :]
                end
                # extrapolate top
                for i = 1:2
                    @inbounds tmp[n[1]+3+i, :, :] = 2.0 .* tmp[n[1]+3, :, :] .- tmp[n[1]+3-i, :, :]
                end
            end

            Threads.@threads for k = 1:n[3]
                for j = 1:n[2]
                    @simd for i = 4:n[1]+3
                        @inbounds out[i, j, k] =
                            a * (tmp[i-1, j, k] + tmp[i, j, k]) +
                            b * (tmp[i-2, j, k] + tmp[i+1, j, k]) +
                            c * (tmp[i-3, j, k] + tmp[i+2, j, k])
                    end
                end
            end
            return out[4:end-2, :, :]
        else # order 2
            out = zeros(T, (n[1], n[2], n[3]))
            Threads.@threads for k = 1:n[3]
                for j = 1:n[2]
                    @simd for i = 2:n[1]
                        @inbounds out[i, j, k] = 0.5f0 * (arr[i-1, j, k] + arr[i, j, k])
                    end
                end
            end
            if periodic
                out[:, :, 1] .= 0.5f0 .* (arr[:, :, end] .+ arr[:, :, 1])
            else
                out[:, :, 1] .= arr[:, :, 1]
            end
            return out
        end
    end
end

"""
	br_dxdn(
		arr::Array{T,3},
		dz::Vector{T}, 
		periodic::Bool=false, 
		order::Int=6
		)

Computes the spatial derivative in the x-direction of every entry in `arr` 
shifted a half grid point downwards. Defaults to the 6th order accurate Bifrost 
derivative with `order=6`, optional 2nd order accurate derivative with keyword 
`order=2`
"""
function br_dxdn(
    arr::Array{T,3},
    dx::Vector{T},
    periodic::Bool=true,
    order::Int=6
) where {T<:AbstractFloat}

    n = size(arr)
    @assert length(n) == 3
    if (n[1] == 1)
        return zeros(T, (n[1], n[2], n[3]))
    else
        if order == 6

            c = (-1.0 + (3.0^5 - 3.0) / (3.0^3 - 3.0)) / (5.0^5 - 5.0 - 5.0 * (3.0^5 - 3))
            b = (-1.0 - 120.0 * c) / 24.0
            a = (1.0 - 3.0 * b - 5.0 * c)

            # inflate matrix in x direction

            out = zeros(T, (n[1] + 5, n[2], n[3]))
            tmp = zeros(T, (n[1] + 5, n[2], n[3]))
            tmp[4:n[1]+3, :, :] .= arr[:, :, :]

            if periodic
                tmp[1:3, :, :] .= arr[end-2:end, :, :]
                tmp[n[1]+4:end, :, :] .= arr[1:2, :, :]
            else
                # extrapolate bottom
                for i = 1:3
                    @inbounds tmp[4-i, :, :] = 2.0 .* tmp[4, :, :] .- tmp[4+i, :, :]
                end
                # extrapolate top
                for i = 1:2
                    @inbounds tmp[n[1]+3+i, :, :] = 2.0 .* tmp[n[1]+3, :, :] .- tmp[n[1]+3-i, :, :]
                end
            end

            Threads.@threads for k = 1:n[3]
                for j = 1:n[2]
                    @simd for i = 4:n[1]+3
                        @inbounds out[i, j, k] =
                            dx[i-3] * (
                                a * (tmp[i, j, k] - tmp[i-1, j, k]) +
                                b * (tmp[i+1, j, k] - tmp[i-2, j, k]) +
                                c * (tmp[i+2, j, k] - tmp[i-3, j, k])
                            )
                    end
                end
            end
            return out[4:end-2, :, :]
        else # oder 2
            out = zeros(T, (n[1], n[2], n[3]))
            Threads.@threads for k = 1:n[3]
                for j = 1:n[2]
                    @simd for i = 2:n[1]
                        @inbounds out[i, j, k] = dx[i] * (arr[i, j, k] - arr[i-1, j, k])
                    end
                end
            end
            if periodic
                out[1, :, :] .= dx[1] .* (arr[1, :, :] .- arr[end, :, :])
            end
            return out
        end
    end
end

"""
    function br_yup(
        arr::Array{T,3},
        periodic::Bool=true,
        order::Int=6
    ) where {T<:AbstractFloat}

Stagger operation on `arr` by a 5th order polynomial interpolation, 
shifting the variable half a grid point upwards in the y-direction
"""
function br_yup(
    arr::Array{T,3},
    periodic::Bool=true,
    order::Int=6
) where {T<:AbstractFloat}

    n = size(arr)
    @assert length(n) == 3
    if (n[2] == 1)
        return arr[:, :, :]
    else

        if order == 6

            c = 3.0 / 256.0
            b = -25.0 / 256.0
            a = 0.5 - b - c

            # inflate matrix in x direction

            out = zeros(T, (n[1], n[2] + 5, n[3]))
            tmp = zeros(T, (n[1], n[2] + 5, n[3]))
            tmp[:, 3:n[2]+2, :] .= arr[:, :, :]

            if periodic
                tmp[:, 1:2, :] .= arr[:, end-1:end, :]
                tmp[:, n[2]+3:end, :] .= arr[:, 1:3, :]
            else
                # extrapolate bottom
                for j = 1:2
                    @inbounds tmp[:, 3-j, :] = 2.0 .* tmp[:, 3, :] .- tmp[:, 3+j, :]
                end
                # extrapolate top
                for j = 1:3
                    @inbounds tmp[:, n[2]+2+j, :] = 2.0 .* tmp[:, n[2]+2, :] .- tmp[:, n[2]+2-j, :]
                end
            end

            Threads.@threads for k = 1:n[3]
                for j = 3:n[2]+2
                    @simd for i = 1:n[1]
                        @inbounds out[i, j, k] =
                            a * (tmp[i, j, k] + tmp[i, j+1, k]) +
                            b * (tmp[i, j-1, k] + tmp[i, j+2, k]) +
                            c * (tmp[i, j-2, k] + tmp[i, j+3, k])
                    end
                end
            end
            return out[:, 3:end-3, :]
        else # order 2
            out = zeros(T, (n[1], n[2], n[3]))
            Threads.@threads for k = 1:n[3]
                for j = 1:n[2]-1
                    @simd for i = 1:n[1]
                        @inbounds out[i, j, k] = 0.5f0 * (arr[i, j, k] + arr[i, j+1, k])
                    end
                end
            end
            if periodic
                out[:, end, :] .= 0.5f0 .* (arr[:, end, :] .+ arr[:, 1, :])
            else
                out[:, end, :] .= arr[:, end, :]
            end
            return out
        end
    end
end

"""
	br_dyup(
		arr::Array{T,3},
		dz::Vector{T}, 
		periodic::Bool=false, 
		order::Int=6
		)

Computes the spatial derivative in the y-direction of every entry in `arr` 
shifted a half grid point upwards. Defaults to the 6th order accurate Bifrost 
derivative with `order=6`, optional 2nd order accurate derivative with keyword 
`order=2`
"""
function br_dyup(
    arr::Array{T,3},
    dy::Vector{T},
    periodic::Bool=true,
    order::Int=6
) where {T<:AbstractFloat}

    n = size(arr)
    @assert length(n) == 3
    if (n[2] == 1)
        return zeros(T, (n[1], n[2], n[3]))
    else
        if order == 6

            c = (-1.0 + (3.0^5 - 3.0) / (3.0^3 - 3.0)) / (5.0^5 - 5.0 - 5.0 * (3.0^5 - 3))
            b = (-1.0 - 120.0 * c) / 24.0
            a = (1.0 - 3.0 * b - 5.0 * c)

            # inflate matrix in x direction

            out = zeros(T, (n[1], n[2] + 5, n[3]))
            tmp = zeros(T, (n[1], n[2] + 5, n[3]))
            tmp[:, 3:n[2]+2, :] = arr[:, :, :]

            if periodic
                tmp[:, 1:2, :] .= arr[:, end-1:end, :]
                tmp[:, n[2]+3:end, :] .= arr[:, 1:3, :]
            else
                # extrapolate bottom
                for j = 1:2
                    @inbounds tmp[:, 3-j, :] .= 2.0 .* tmp[:, 3, :] .- tmp[:, 3+j, :]
                end
                # extrapolate top
                for j = 1:3
                    @inbounds tmp[:, n[2]+2+j, :] .= 2.0 .* tmp[:, n[2]+2, :] .- tmp[:, n[2]+2-j, :]
                end
            end

            Threads.@threads for k = 1:n[3]
                for j = 3:n[2]+2
                    @simd for i = 1:n[1]
                        @inbounds out[i, j, k] =
                            dy[j-2] * (
                                a * (tmp[i, j+1, k] - tmp[i, j, k]) +
                                b * (tmp[i, j+2, k] - tmp[i, j-1, k]) +
                                c * (tmp[i, j+3, k] - tmp[i, j-2, k])
                            )
                    end
                end
            end
            return out[:, 3:end-3, :]
        else # order 2
            out = zeros(T, (n[1], n[2], n[3]))
            Threads.@threads for k = 1:n[3]
                for j = 1:n[2]-1
                    @simd for i = 1:n[1]
                        @inbounds out[i, j, k] = dy[j] * (arr[i, j+1, k] - arr[i, j, k])
                    end
                end
            end
            if periodic
                out[:, end, :] .= dy[end] .* (arr[:, 1, :] - arr[:, end, :])
            end
            return out
        end
    end
end

"""
    function br_ydn(
        arr::Array{T,3},
        periodic::Bool=true,
        order::Int=6
    ) where {T<:AbstractFloat}

Stagger operation on `arr` by a 5th order polynomial interpolation, 
shifting the variable half a grid point downwards in the y-direction
"""
function br_ydn(
	arr::Array{T,3}, 
	periodic::Bool=true, 
	order::Int=6
) where {T<:AbstractFloat}

    n = size(arr)
    @assert length(n) == 3
    if (n[2] == 1)
        return arr[:, :, :]
    else
        if order == 6

            c = 3.0 / 256.0
            b = -25.0 / 256.0
            a = 0.5 - b - c

            # inflate matrix in x direction

            out = zeros(T, (n[1], n[2] + 5, n[3]))
            tmp = zeros(T, (n[1], n[2] + 5, n[3]))
            tmp[:, 4:n[2]+3, :] .= arr[:, :, :]

            if periodic
                tmp[:, 1:3, :] .= arr[:, end-2:end, :]
                tmp[:, n[2]+4:end, :] .= arr[:, 1:2, :]
            else
                # extrapolate bottom
                for j = 1:3
                    @inbounds tmp[:, 4-j, :] .= 2.0 .* tmp[:, 4, :] .- tmp[:, 4+j, :]
                end
                # extrapolate top
                for j = 1:2
                    @inbounds tmp[:, n[2]+3+j, :] .= 2.0 .* tmp[:, n[2]+3, :] .- tmp[:, n[2]+3-j, :]
                end
            end

            Threads.@threads for k = 1:n[3]
                for j = 4:n[2]+3
                    @simd for i = 1:n[1]
                        @inbounds out[i, j, k] =
                            a * (tmp[i, j-1, k] + tmp[i, j, k]) +
                            b * (tmp[i, j-2, k] + tmp[i, j+1, k]) +
                            c * (tmp[i, j-3, k] + tmp[i, j+2, k])
                    end
                end
            end
            return out[:, 4:end-2, :]
        else # order 2
            out = zeros(T, (n[1], n[2], n[3]))
            Threads.@threads for k = 1:n[3]
                for j = 2:n[2]
                    @simd for i = 1:n[1]
                        @inbounds out[i, j, k] = 0.5f0 * (arr[i, j-1, k] + arr[i, j, k])
                    end
                end
            end
            if periodic
                out[:, 1, :] .= 0.5f0 .* (arr[:, end, :] .+ arr[:, 1, :])
            else
                out[:, 1, :] .= arr[:, 1, :]
            end
            return out
        end
    end
end

"""
	br_dydn(
		arr::Array{T,3},
		dz::Vector{T}, 
		periodic::Bool=false, 
		order::Int=6
		)

Computes the spatial derivative in the y-direction of every entry in `arr` 
shifted a half grid point downwards. Defaults to the 6th order accurate Bifrost 
derivative with `order=6`, optional 2nd order accurate derivative with keyword 
`order=2`
"""
function br_dydn(
    arr::Array{T,3},
    dy::Vector{T},
    periodic::Bool=true,
    order::Int=6
) where {T<:AbstractFloat}

    n = size(arr)
    @assert length(n) == 3
    if (n[2] == 1)
        return zeros(T, (n[1], n[2], n[3]))
    else
        if order == 6

            c = (-1.0 + (3.0^5 - 3.0) / (3.0^3 - 3.0)) / (5.0^5 - 5.0 - 5.0 * (3.0^5 - 3))
            b = (-1.0 - 120.0 * c) / 24.0
            a = (1.0 - 3.0 * b - 5.0 * c)

            # inflate matrix in x direction

            out = zeros(T, (n[1], n[2] + 5, n[3]))
            tmp = zeros(T, (n[1], n[2] + 5, n[3]))
            tmp[:, 4:n[2]+3, :] .= arr[:, :, :]

            if periodic
                tmp[:, 1:3, :] .= arr[:, end-2:end, :]
                tmp[:, n[2]+4:end, :] .= arr[:, 1:2, :]
            else
                # extrapolate bottom
                for j = 1:3
                    @inbounds tmp[:, 4-j, :] .= 2.0 .* tmp[:, 4, :] .- tmp[:, 4+j, :]
                end
                # extrapolate top
                for j = 1:2
                    @inbounds tmp[:, n[2]+3+j, :] .= 2.0 .* tmp[:, n[2]+3, :] .- tmp[:, n[2]+3-j, :]
                end
            end

            Threads.@threads for k = 1:n[3]
                for j = 4:n[2]+3
                    @simd for i = 1:n[1]
                        @inbounds out[i, j, k] =
                            dy[j-3] * (
                                a * (tmp[i, j, k] - tmp[i, j-1, k]) +
                                b * (tmp[i, j+1, k] - tmp[i, j-2, k]) +
                                c * (tmp[i, j+2, k] - tmp[i, j-3, k])
                            )
                    end
                end
            end
            return out[:, 4:end-2, :]
        else # order 2
            out = zeros(T, (n[1], n[2], n[3]))
            Threads.@threads for k = 1:n[3]
                for j = 2:n[2]
                    @simd for i = 1:n[1]
                        @inbounds out[i, j, k] = dy[j] * (arr[i, j, k] - arr[i, j-1, k])
                    end
                end
            end
            if periodic
                out[:, 1, :] .= dy[1] .* (arr[:, 1, :] .- arr[:, end, :])
            end
            return out
        end
    end
end

"""
    function br_zup(
        arr::Array{T,3},
        periodic::Bool=true,
        order::Int=6
    ) where {T<:AbstractFloat}

Stagger operation on `arr` by a 5th order polynomial interpolation, 
shifting the variable half a grid point upwards in the z-direction
"""
function br_zup(
	arr::Array{T,3}, 
	periodic::Bool=false, 
	order::Int=6
) where {T<:AbstractFloat}

    n = size(arr)
    @assert length(n) == 3
    if (n[3] == 1)
        return arr[:, :, :]
    else

        if order == 6

            c = 3.0 / 256.0
            b = -25.0 / 256.0
            a = 0.5 - b - c

            # inflate matrix in x direction

            out = zeros(T, (n[1], n[2], n[3] + 5))
            tmp = zeros(T, (n[1], n[2], n[3] + 5))
            tmp[:, :, 3:n[3]+2] .= arr[:, :, :]

            if periodic
                tmp[:, :, 1:2] .= arr[:, :, end-1:end]
                tmp[:, :, n[3]+3:end] .= arr[:, :, 1:3]
            else
                # extrapolate bottom
                for k = 1:2
                    @inbounds tmp[:, :, 3-k] = 2.0 .* tmp[:, :, 3] .- tmp[:, :, 3+k]
                end
                # extrapolate top
                for k = 1:3
                    @inbounds tmp[:, :, n[3]+2+k] = 2.0 .* tmp[:, :, n[3]+2] .- tmp[:, :, n[3]+2-k]
                end
            end

            Threads.@threads for k = 3:n[3]+2
                for j = 1:n[2]
                    @simd for i = 1:n[1]
                        @inbounds out[i, j, k] =
                            a * (tmp[i, j, k] + tmp[i, j, k+1]) +
                            b * (tmp[i, j, k-1] + tmp[i, j, k+2]) +
                            c * (tmp[i, j, k-2] + tmp[i, j, k+3])
                    end
                end
            end
            return out[:, :, 3:end-3]
        else # order 2
            out = zeros(T, (n[1], n[2], n[3]))
            Threads.@threads for k = 1:n[3]-1
                for j = 1:n[2]
                    @simd for i = 1:n[1]
                        @inbounds out[i, j, k] = 0.5f0 * (arr[i, j, k] + arr[i, j, k+1])
                    end
                end
            end
            if periodic
                out[:, :, end] .= 0.5f0 .* (arr[:, :, end] .+ arr[:, :, 1])
            else
                out[:, :, end] .= arr[:, :, end]
            end
            return out
        end
    end
end

function br_zup(
	arr::Array{T,3}, 
    slicez::AbstractVector{<:Integer},
	periodic::Bool=false, 
	order::Int=6
) where {T<:AbstractFloat}

    n = size(arr)
    @assert length(n) == 3
    if (n[3] == 1)
        return arr[:, :, :]
    else

        if order == 6

            c = 3.0 / 256.0
            b = -25.0 / 256.0
            a = 0.5 - b - c

            out = similar(arr, n[1], n[2], length(slicez))
            tmp = zeros(T, (n[1], n[2], n[3] + 5))
            tmp[:, :, 3:n[3]+2] .= arr[:, :, :]

            if periodic
                tmp[:, :, 1:2] .= arr[:, :, end-1:end]
                tmp[:, :, n[3]+3:end] .= arr[:, :, 1:3]
            else
                # extrapolate bottom
                for k = 1:2
                    @inbounds tmp[:, :, 3-k] = 2.0 .* tmp[:, :, 3] .- tmp[:, :, 3+k]
                end
                # extrapolate top
                for k = 1:3
                    @inbounds tmp[:, :, n[3]+2+k] = 2.0 .* tmp[:, :, n[3]+2] .- tmp[:, :, n[3]+2-k]
                end
            end

            slicez = slicez .+ 2
            Threads.@threads for (k,slice) in collect(enumerate(slicez))
                for j = 1:n[2]
                    @simd for i = 1:n[1]
                        @inbounds out[i, j, k] =
                            a * (tmp[i, j, slice] + tmp[i, j, slice+1]) +
                            b * (tmp[i, j, slice-1] + tmp[i, j, slice+2]) +
                            c * (tmp[i, j, slice-2] + tmp[i, j, slice+3])
                    end
                end
            end
            return out
        else # order 2
            out = similar(arr, (n[1], n[2], length(slicez)))
            
            if n[3] in slicez
                index = findfirst(x->x==n[3], slicez)
                if periodic
                    out[:, :, index] .= 0.5f0 .* (arr[:, :, end] .+ arr[:, :, 1])
                else
                    out[:, :, index] .= arr[:, :, end]
                end
                slicez = slicez[∉(index).(1:end)]
            end
            
            Threads.@threads for (k,slice) in collect(enumerate(slicez))
                for j = 1:n[2]
                    @simd for i = 1:n[1]
                        @inbounds out[i, j, k] = 0.5f0 * (arr[i, j, slice] + arr[i, j, slice+1])
                    end
                end
            end
            return out
        end
    end
end

"""
	br_dzup(
		arr::Array{T,3},
		dz::Vector{T}, 
		periodic::Bool=false, 
		order::Int=6
		)

Computes the spatial derivative in the z-direction of every entry in `arr` 
shifted a half grid point upwards. Defaults to the 6th order accurate Bifrost
derivative with `order=6`, optional 2nd order accurate derivative with keyword 
`order=2`
"""
function br_dzup(
    arr::Array{T,3},
    dz::Vector{T},
    periodic::Bool=false,
    order::Int=6
) where {T<:AbstractFloat}

    n = size(arr)
    @assert length(n) == 3
    if (n[3] == 1)
        return zeros(T, (n[1], n[2], n[3]))
    else
        if order == 6

            c = (-1.0 + (3.0^5 - 3.0) / (3.0^3 - 3.0)) / (5.0^5 - 5.0 - 5.0 * (3.0^5 - 3))
            b = (-1.0 - 120.0 * c) / 24.0
            a = (1.0 - 3.0 * b - 5.0 * c)

            # inflate matrix in x direction

            out = zeros(T, (n[1], n[2], n[3] + 5))
            tmp = zeros(T, (n[1], n[2], n[3] + 5))
            tmp[:, :, 3:n[3]+2] .= arr[:, :, :]

            if periodic
                tmp[:, :, 1:2] .= arr[:, :, end-1:end]
                tmp[:, :, n[3]+3:end] .= arr[:, :, 1:3]
            else
                # extrapolate bottom
                for k = 1:2
                    @inbounds tmp[:, :, 3-k] .= 2.0 .* tmp[:, :, 3] .- tmp[:, :, 3+k]
                end
                # extrapolate top
                for k = 1:3
                    @inbounds tmp[:, :, n[3]+2+k] .= 2.0 .* tmp[:, :, n[3]+2] .- tmp[:, :, n[3]+2-k]
                end
            end

            Threads.@threads for k = 3:n[3]+2
                for j = 1:n[2]
                    @simd for i = 1:n[1]
                        @inbounds out[i, j, k] =
                            dz[k-2] * (
                                a * (tmp[i, j, k+1] - tmp[i, j, k]) +
                                b * (tmp[i, j, k+2] - tmp[i, j, k-1]) +
                                c * (tmp[i, j, k+3] - tmp[i, j, k-2])
                            )
                    end
                end
            end
            return out[:, :, 3:end-3]
        else # order 2
            out = zeros(T, (n[1], n[2], n[3]))
            Threads.@threads for k = 1:n[3]-1
                for j = 1:n[2]
                    @simd for i = 1:n[1]
                        @inbounds out[i, j, k] = dz[k] * (arr[i, j, k+1] - arr[i, j, k])
                    end
                end
            end
            if periodic
                out[:, :, end] .= dz[end] .* (arr[:, :, 1] .- arr[:, :, end])
            end
            return out
        end
    end
end

"""
    function br_zdn(
        arr::Array{T,3},
        periodic::Bool=true,
        order::Int=6
    ) where {T<:AbstractFloat}

Stagger operation on `arr` by a 5th order polynomial interpolation, 
shifting the variable half a grid point downwards in the z-direction
"""
function br_zdn(
	arr::Array{T,3}, 
	periodic::Bool=false, 
	order::Int=6
) where {T<:AbstractFloat}
    
	n = size(arr)
    @assert length(n) == 3
    if (n[3] == 1)
        return arr[:, :, :]
    else
        if order == 6

            c = 3.0 / 256.0
            b = -25.0 / 256.0
            a = 0.5 - b - c

            # inflate matrix in x direction

            out = zeros(T, (n[1], n[2], n[3] + 5))
            tmp = zeros(T, (n[1], n[2], n[3] + 5))
            tmp[:, :, 4:n[3]+3] .= arr[:, :, :]

            if periodic
                tmp[:, :, 1:3] .= arr[:, :, end-2:end]
                tmp[:, :, n[3]+4:end] .= arr[:, :, 1:2]
            else
                # extrapolate bottom
                for k = 1:3
                    @inbounds tmp[:, :, 4-k] .= 2.0 .* tmp[:, :, 4] .- tmp[:, :, 4+k]
                end
                # extrapolate top
                for k = 1:2
                    @inbounds tmp[:, :, n[3]+3+k] .= 2.0 .* tmp[:, :, n[3]+3] .- tmp[:, :, n[3]+3-k]
                end
            end

            Threads.@threads for k = 4:n[3]+3
                for j = 1:n[2]
                    @simd for i = 1:n[1]
                        @inbounds out[i, j, k] =
                            a * (tmp[i, j, k-1] + tmp[i, j, k]) +
                            b * (tmp[i, j, k-2] + tmp[i, j, k+1]) +
                            c * (tmp[i, j, k-3] + tmp[i, j, k+2])
                    end
                end
            end
            return out[:, :, 4:end-2]
        else # order 2
            out = zeros(T, (n[1], n[2], n[3]))
            Threads.@threads for k = 2:n[3]
                for j = 1:n[2]
                    @simd for i = 1:n[1]
                        @inbounds out[i, j, k] = 0.5f0 * (arr[i, j, k-1] + arr[i, j, k])
                    end
                end
            end
            if periodic
                out[:, :, 1] .= 0.5f0 .* (arr[:, :, end] .+ arr[:, :, 1])
            else
                out[:, :, 1] .= arr[:, :, 1]
            end
            return out
        end
    end
end

function br_zdn(
	arr::Array{T,3},
    slicez::AbstractVector{<:Integer},
	periodic::Bool=false, 
	order::Int=6
) where {T<:AbstractFloat}
    
	n = size(arr)
    @assert length(n) == 3
    if (n[3] == 1)
        return arr[:, :, :]
    else
        if order == 6
            
            c = 3.0 / 256.0
            b = -25.0 / 256.0
            a = 0.5 - b - c

            out = similar(arr, n[1], n[2], length(slicez))
            tmp = zeros(T, (n[1], n[2], n[3] + 5))
            tmp[:, :, 4:n[3]+3] .= arr[:, :, :]

            if periodic
                tmp[:, :, 1:3] .= arr[:, :, end-2:end]
                tmp[:, :, n[3]+4:end] .= arr[:, :, 1:2]
            else
                # extrapolate bottom
                for k = 1:3
                    @inbounds tmp[:, :, 4-k] .= 2.0 .* tmp[:, :, 4] .- tmp[:, :, 4+k]
                end
                # extrapolate top
                for k = 1:2
                    @inbounds tmp[:, :, n[3]+3+k] .= 2.0 .* tmp[:, :, n[3]+3] .- tmp[:, :, n[3]+3-k]
                end
            end

            slicez = slicez .+ 3
            Threads.@threads for (k,slice) in collect(enumerate(slicez))
                for j = 1:n[2]
                    @simd for i = 1:n[1]
                        @inbounds out[i, j, k] =
                            a * (tmp[i, j, slice-1] + tmp[i, j, slice]) +
                            b * (tmp[i, j, slice-2] + tmp[i, j, slice+1]) +
                            c * (tmp[i, j, slice-3] + tmp[i, j, slice+2])
                    end
                end
            end
            return out
        else # order 2
            out = similar(arr, (n[1], n[2], length(slicez)))
            
            if 1 in slicez
                index = findfirst(x->x==n[3], slicez)
                if periodic
                    out[:, :, index] .= 0.5f0 .* (arr[:, :, end] .+ arr[:, :, 1])
                else
                    out[:, :, index] .= arr[:, :, 1]
                end
                slicez = slicez[∉(index).(1:end)]
            end
            
            Threads.@threads for (k,slice) in collect(enumerate(slicez))
                for j = 1:n[2]
                    @simd for i = 1:n[1]
                        @inbounds out[i, j, k] = 0.5f0 * (arr[i, j, slice-1] + arr[i, j, slice])
                    end
                end
            end
            return out
        end
    end
end

"""
	br_dzdn(
		arr::Array{T,3},
		dz::Vector{T}, 
		periodic::Bool=false, 
		order::Int=6
		)

Computes the spatial derivative in the z-direction of every entry in `arr` 
shifted a half grid point downwards. Defaults to the 6th order accurate Bifrost 
derivative with `order=6`, optional 2nd order accurate derivative with keyword 
`order=2`
"""
function br_dzdn(
    arr::Array{T,3},
    dz::Vector{T},
    periodic::Bool=false,
    order::Int=6
) where {T<:AbstractFloat}
    n = size(arr)
    @assert length(n) == 3
    if (n[3] == 1)
        return zeros(T, (n[1], n[2], n[3]))
    else
        if order == 6

            c = (-1.0 + (3.0^5 - 3.0) / (3.0^3 - 3.0)) / (5.0^5 - 5.0 - 5.0 * (3.0^5 - 3))
            b = (-1.0 - 120.0 * c) / 24.0
            a = (1.0 - 3.0 * b - 5.0 * c)

            # inflate matrix in x direction
            out = zeros(T, (n[1], n[2], n[3] + 5))
            tmp = zeros(T, (n[1], n[2], n[3] + 5))
            tmp[:, :, 4:n[3]+3] .= arr[:, :, :]

            if periodic
                tmp[:, :, 1:3] .= arr[:, :, end-2:end]
                tmp[:, :, n[3]+4:end] .= arr[:, :, 1:2]
            else
                # extrapolate bottom
                for k = 1:3
                    @inbounds tmp[:, :, 4-k] .= 2.0 .* tmp[:, :, 4] .- tmp[:, :, 4+k]
                end
                # extrapolate top
                for k = 1:2
                    @inbounds tmp[:, :, n[3]+3+k] .= 2.0 .* tmp[:, :, n[3]+3] .- tmp[:, :, n[3]+3-k]
                end
            end

            Threads.@threads for k = 4:n[3]+3
                for j = 1:n[2]
                    @simd for i = 1:n[1]
                        @inbounds out[i, j, k] =
                            dz[k-3] * (
                                a * (tmp[i, j, k] - tmp[i, j, k-1]) +
                                b * (tmp[i, j, k+1] - tmp[i, j, k-2]) +
                                c * (tmp[i, j, k+2] - tmp[i, j, k-3])
                            )
                    end
                end
            end
            return out[:, :, 4:end-2]
        else # order 2
            out = zeros(T, (n[1], n[2], n[3]))
            Threads.@threads for k = 2:n[3]
                for j = 1:n[2]
                    @simd for i = 1:n[1]
                        @inbounds out[i, j, k] = dz[k] * (arr[i, j, k] - arr[i, j, k-1])
                    end
                end
            end
            if periodic
                out[:, :, 1] .= dz[1] .* (arr[:, :, 1] .- arr[:, :, end])
            end
            return out
        end
    end
end


