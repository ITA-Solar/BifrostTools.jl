
function br_cdivB_clean_du(
    m::BifrostMesh,
    bx::Array{T,3},
    by::Array{T,3},
    bz::Array{T,3};
    order::Int=6,
    nit::Int=256,
    damp::Float32=12.0f0) where {T<:AbstractFloat} 
    
    divb = br_divdn(m, bx, by, bz, order=6);
    @printf "[DivB] min: %e :: max %e\n" minimum(divb) maximum(divb)
    res = divb;

    # Poisson solver

    phi = divb ./ 12.0f0
    for i in 1:1024

            res = br_claplace_du(phi) .- divb
            phi .= phi .+ (res ./ 12.0f0)
            if i % 32 == 0
                        @printf "It: %03d, min = %e / max %e \n" i minimum(res) maximum(res)
            end

    end

    b_fix = br_cgrad_up(phi);
    fbx = bx .- b_fix[1]
    fby = by .- b_fix[2]
    fbz = bz .- b_fix[3]

    fdivb = br_cdivdn([fbx, fby, fbz]);

    @printf "(old) divB [min/max] %e / %e \n" minimum(divb) maximum(divb)
    @printf "(new) divB [min/max] %e / %e" minimum(fdivb) maximum(fdivb)
end

function br_cdivB_clean_ud(
    m::BifrostMesh,
    bx::Array{T,3},
    by::Array{T,3},
    bz::Array{T,3};
    order::Int=6,
    nit::Int=256,
    damp::Float32=12.0f0) where {T<:AbstractFloat} 
    
    divb = br_divup(m, bx, by, bz, order=6);
    @printf "[DivB] min: %e :: max %e\n" minimum(divb) maximum(divb)
    res = divb;

    # Poisson solver

    phi = divb ./ 12.0f0
    for i in 1:1024

            res = br_claplace_ud(phi) .- divb
            phi .= phi .+ (res ./ 12.0f0)
            if i % 32 == 0
                        @printf "It: %03d, min = %e / max %e \n" i minimum(res) maximum(res)
            end

    end

    b_fix = br_cgrad_dn(phi);
    fbx = bx .- b_fix[1]
    fby = by .- b_fix[2]
    fbz = bz .- b_fix[3]

    fdivb = br_cdivup([fbx, fby, fbz]);

    @printf "(old) divB [min/max] %e / %e \n" minimum(divb) maximum(divb)
    @printf "(new) divB [min/max] %e / %e" minimum(fdivb) maximum(fdivb)
end

# --- aditinal operators

function br_divup(
    m::BifrostMesh,
    x::Array{T,3},
    y::Array{T,3},
    z::Array{T,3};
    order::Int=6
) where {T<:AbstractFloat}
    return br_dxup(x, m.dxidxup, true, order) .+ br_dyup(y, m.dyidyup, true, order) .+ br_dzup(z, m.dzidzup, false, order)
end

function br_divdn(
    m::BifrostMesh,
    x::Array{T,3},
    y::Array{T,3},
    z::Array{T,3};
    order::Int=6
) where {T<:AbstractFloat}
    return br_dxdn(x, m.dxidxdn, true, order) .+ br_dydn(y, m.dyidydn, true, order) .+ br_dzdn(z, m.dzidzdn, false, order)
end

function br_gradup(m::BifrostMesh, x::Array{T,3}; order::Int=6) where {T<:AbstractFloat}
    return [br_dxup(x, m.dxidxup, true, order), br_dyup(x, m.dyidyup, true, order), br_dzup(x, m.dzidzup, false, order)]
end

function br_graddn(m::BifrostMesh, x::Array{T,3}; order::Int=6) where {T<:AbstractFloat}
    return [br_dxdn(x, m.dxidxdn, true, order), br_dydn(x, m.dyidydn, true, order), br_dzdn(x, m.dzidzdn, false, order)]
end

function br_laplacedu(m::BifrostMesh, x::Array{T,3}; order::Int=6) where {T<:AbstractFloat}
    tmp = br_gradup(m, x, order=order)
    return br_divdn(m, tmp[1], tmp[2], tmp[3], order=order)
end

function br_laplaceud(m::BifrostMesh, x::Array{T,3}; order::Int=6) where {T<:AbstractFloat}
    tmp = br_graddn(m, x, order=order)
    return br_divup(m, tmp[1], tmp[2], tmp[3], order=order)
end

function br_poissondu(m::BifrostMesh, x::Array{T,3}; order::Int=6, nit::Int=256, damp::Float32=12.0f0) where {T<:AbstractFloat}
    phi = x ./ damp
    for i in 1:nit
        res = br_laplacedu(m, phi, order=order) .- x
        phi .= phi .+ (res ./ damp)
        @printf "\t[br_laplaceup] It: %d, min = %e / max %e \n" i minimum(res) maximum(res)
    end
    return phi
end

function br_poissonud(m::BifrostMesh, x::Array{T,3}; order::Int=6, nit::Int=256, damp::Float32=12.0f0) where {T<:AbstractFloat}
    phi = x ./ damp
    for i in 1:nit
        res = br_laplaceud(m, phi, order=order) .- x
        phi = phi .+ (res ./ damp)
        @printf "\t[br_laplacedn] It: %d, min = %e / max %e \n" i minimum(res) maximum(res)
    end
    return phi
end

function br_divB_clean_du(
    m::BifrostMesh,
    bx::Array{T,3},
    by::Array{T,3},
    bz::Array{T,3};
    order::Int=6,
    nit::Int=256,
    damp::Float32=12.0f0) where {T<:AbstractFloat} 

    divb = br_divup(m, bx, by, bz, order=order);
    @printf "(before) divB [min/max] %e / %e \n" minimum(divb) maximum(divb)
    
    phi    = br_poissondu(m,divb,order=order,nit=nit,damp=damp);

    b_fix = br_graddn(m, phi, order=order);
    fbx = bx .- b_fix[1]
    fby = by .- b_fix[2]
    fbz = bz .- b_fix[3]
    
    fdivb = br_divup(m, fbx, fby, fbz, order=order);
    @printf "(after) divB [min/max] %e / %e" minimum(fdivb) maximum(fdivb)
    return [fbx, fby, fbz]
end

function br_divB_clean_ud(
    m::BifrostMesh,
    bx::Array{T,3},
    by::Array{T,3},
    bz::Array{T,3};
    order::Int=6,
    nit::Int=256,
    damp::Float32=12.0f0) where {T<:AbstractFloat} 

    divb = br_divdn(m, bx, by, bz, order=order);
    @printf "(before) divB [min/max] %e / %e \n" minimum(divb) maximum(divb)
    
    phi    = br_poissonud(m,divb,order=order,nit=nit,damp=damp);

    b_fix = br_gradup(m, phi, order=order);
    fbx = bx .- b_fix[1]
    fby = by .- b_fix[2]
    fbz = bz .- b_fix[3]
    
    fdivb = br_divdn(m, fbx, fby, fbz, order=order);
    @printf "(after) divB [min/max] %e / %e" minimum(fdivb) maximum(fdivb)
    return [fbx, fby, fbz]
end


