
function cdivB_clean_du(
    m::BifrostMesh,
    bx::Array{T,3},
    by::Array{T,3},
    bz::Array{T,3};
    order::Int=6,
    nit::Int=256,
    damp::Float32=12.0f0) where {T<:AbstractFloat} 
    
    divb = divdn(m, bx, by, bz, order=6);
    @printf "[DivB] min: %e :: max %e\n" minimum(divb) maximum(divb)
    res = divb;

    # Poisson solver

    phi = divb ./ 12.0f0
    for i in 1:1024

            res = claplace_du(phi) .- divb
            phi .= phi .+ (res ./ 12.0f0)
            if i % 32 == 0
                        @printf "It: %03d, min = %e / max %e \n" i minimum(res) maximum(res)
            end

    end

    b_fix = cgrad_up(phi);
    fbx = bx .- b_fix[1]
    fby = by .- b_fix[2]
    fbz = bz .- b_fix[3]

    fdivb = cdivdn([fbx, fby, fbz]);

    @printf "(old) divB [min/max] %e / %e \n" minimum(divb) maximum(divb)
    @printf "(new) divB [min/max] %e / %e" minimum(fdivb) maximum(fdivb)
end

function cdivB_clean_ud(
    m::BifrostMesh,
    bx::Array{T,3},
    by::Array{T,3},
    bz::Array{T,3};
    order::Int=6,
    nit::Int=256,
    damp::Float32=12.0f0) where {T<:AbstractFloat} 
    
    divb = divup(m, bx, by, bz, order=6);
    @printf "[DivB] min: %e :: max %e\n" minimum(divb) maximum(divb)
    res = divb;

    # Poisson solver

    phi = divb ./ 12.0f0
    for i in 1:1024

            res = claplace_ud(phi) .- divb
            phi .= phi .+ (res ./ 12.0f0)
            if i % 32 == 0
                        @printf "It: %03d, min = %e / max %e \n" i minimum(res) maximum(res)
            end

    end

    b_fix = cgrad_dn(phi);
    fbx = bx .- b_fix[1]
    fby = by .- b_fix[2]
    fbz = bz .- b_fix[3]

    fdivb = cdivup([fbx, fby, fbz]);

    @printf "(old) divB [min/max] %e / %e \n" minimum(divb) maximum(divb)
    @printf "(new) divB [min/max] %e / %e" minimum(fdivb) maximum(fdivb)
end

# --- aditinal operators

function divup(
    m::BifrostMesh,
    x::Array{T,3},
    y::Array{T,3},
    z::Array{T,3};
    order::Int=6
) where {T<:AbstractFloat}
    return dxup(x, m.dxidxup, true, order) .+ dyup(y, m.dyidyup, true, order) .+ dzup(z, m.dzidzup, false, order)
end

function divdn(
    m::BifrostMesh,
    x::Array{T,3},
    y::Array{T,3},
    z::Array{T,3};
    order::Int=6
) where {T<:AbstractFloat}
    return dxdn(x, m.dxidxdn, true, order) .+ dydn(y, m.dyidydn, true, order) .+ dzdn(z, m.dzidzdn, false, order)
end

function gradup(m::BifrostMesh, x::Array{T,3}; order::Int=6) where {T<:AbstractFloat}
    return [dxup(x, m.dxidxup, true, order), dyup(x, m.dyidyup, true, order), dzup(x, m.dzidzup, false, order)]
end

function graddn(m::BifrostMesh, x::Array{T,3}; order::Int=6) where {T<:AbstractFloat}
    return [dxdn(x, m.dxidxdn, true, order), dydn(x, m.dyidydn, true, order), dzdn(x, m.dzidzdn, false, order)]
end

function laplacedu(m::BifrostMesh, x::Array{T,3}; order::Int=6) where {T<:AbstractFloat}
    tmp = gradup(m, x, order=order)
    return divdn(m, tmp[1], tmp[2], tmp[3], order=order)
end

function laplaceud(m::BifrostMesh, x::Array{T,3}; order::Int=6) where {T<:AbstractFloat}
    tmp = graddn(m, x, order=order)
    return divup(m, tmp[1], tmp[2], tmp[3], order=order)
end

function poissondu(m::BifrostMesh, x::Array{T,3}; order::Int=6, nit::Int=256, damp::Float32=12.0f0) where {T<:AbstractFloat}
    phi = x ./ damp
    for i in 1:nit
        res = laplacedu(m, phi, order=order) .- x
        phi .= phi .+ (res ./ damp)
        @printf "\t[laplaceup] It: %d, min = %e / max %e \n" i minimum(res) maximum(res)
    end
    return phi
end

function poissonud(m::BifrostMesh, x::Array{T,3}; order::Int=6, nit::Int=256, damp::Float32=12.0f0) where {T<:AbstractFloat}
    phi = x ./ damp
    for i in 1:nit
        res = laplaceud(m, phi, order=order) .- x
        phi = phi .+ (res ./ damp)
        @printf "\t[laplacedn] It: %d, min = %e / max %e \n" i minimum(res) maximum(res)
    end
    return phi
end

function divB_clean_du(
    m::BifrostMesh,
    bx::Array{T,3},
    by::Array{T,3},
    bz::Array{T,3};
    order::Int=6,
    nit::Int=256,
    damp::Float32=12.0f0) where {T<:AbstractFloat} 

    divb = divup(m, bx, by, bz, order=order);
    @printf "(before) divB [min/max] %e / %e \n" minimum(divb) maximum(divb)
    
    phi    = poissondu(m,divb,order=order,nit=nit,damp=damp);

    b_fix = graddn(m, phi, order=order);
    fbx = bx .- b_fix[1]
    fby = by .- b_fix[2]
    fbz = bz .- b_fix[3]
    
    fdivb = divup(m, fbx, fby, fbz, order=order);
    @printf "(after) divB [min/max] %e / %e" minimum(fdivb) maximum(fdivb)
    return [fbx, fby, fbz]
end

function divB_clean_ud(
    m::BifrostMesh,
    bx::Array{T,3},
    by::Array{T,3},
    bz::Array{T,3};
    order::Int=6,
    nit::Int=256,
    damp::Float32=12.0f0) where {T<:AbstractFloat} 

    divb = divdn(m, bx, by, bz, order=order);
    @printf "(before) divB [min/max] %e / %e \n" minimum(divb) maximum(divb)
    
    phi    = poissonud(m,divb,order=order,nit=nit,damp=damp);

    b_fix = gradup(m, phi, order=order);
    fbx = bx .- b_fix[1]
    fby = by .- b_fix[2]
    fbz = bz .- b_fix[3]
    
    fdivb = divdn(m, fbx, fby, fbz, order=order);
    @printf "(after) divB [min/max] %e / %e" minimum(fdivb) maximum(fdivb)
    return [fbx, fby, fbz]
end


