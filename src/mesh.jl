"""
        BifrostMesh

Stores Bifrost grid information in struct
"""
struct BifrostMesh
    mx::Int64
    x::Vector{Float32}
    xmdn::Vector{Float32}
    dxidxup::Vector{Float32}
    dxidxdn::Vector{Float32}
    my::Int64
    y::Vector{Float32}
    ymdn::Vector{Float32}
    dyidyup::Vector{Float32}
    dyidydn::Vector{Float32}
    mz::Int64
    z::Vector{Float32}
    zmdn::Vector{Float32}
    dzidzup::Vector{Float32}
    dzidzdn::Vector{Float32}
    n::Int64
    
    function BifrostMesh(expdir::String)

        # Check if `expdir` is path to mesh_file or folder of experiment
        expname = splitpath(expdir)[end]
        if split(expname, ".")[end] == "mesh"
            mesh_file = expdir
        else
            mesh_file = joinpath(expdir, expname*".mesh")
        end

        f = open(mesh_file, "r")
        l = readlines(f)
        # -- x direction
        mx = parse.(Int64, l[1])
        x = parse.(Float32, split(l[2]))
        xmdn = parse.(Float32, split(l[3]))
        dxidxup = parse.(Float32, split(l[4]))
        dxidxdn = parse.(Float32, split(l[5]))
        # -- y direction
        my = parse.(Int64, l[6])
        y = parse.(Float32, split(l[7]))
        ymdn = parse.(Float32, split(l[8]))
        dyidyup = parse.(Float32, split(l[9]))
        dyidydn = parse.(Float32, split(l[10]))
        # -- z direction
        mz = parse.(Int64, l[11])
        z = parse.(Float32, split(l[12]))
        zmdn = parse.(Float32, split(l[13]))
        dzidzup = parse.(Float32, split(l[14]))
        dzidzdn = parse.(Float32, split(l[15]))
        new(
            mx,
            x,
            xmdn,
            dxidxup,
            dxidxdn,
            my,
            y,
            ymdn,
            dyidyup,
            dyidydn,
            mz,
            z,
            zmdn,
            dzidzup,
            dzidzdn,
            mx * my * mz
        )
    end
end


function mesh2file(M::BifrostMesh, file_name::String ="bifrost.mesh")
    open(file_name,"w") do io
        println(io, @sprintf "%d" M.mx)
        println(io, join([@sprintf "%e" x for x in M.x], " "))
        println(io, join([@sprintf "%e" x for x in M.xmdn], " "))
        println(io, join([@sprintf "%e" x for x in M.dxidxup], " "))
        println(io, join([@sprintf "%e" x for x in M.dxidxdn], " "))
        println(io, @sprintf "%d" M.my)
        println(io, join([@sprintf "%e" x for x in M.y], " "))
        println(io, join([@sprintf "%e" x for x in M.ymdn], " "))
        println(io, join([@sprintf "%e" x for x in M.dyidyup], " "))
        println(io, join([@sprintf "%e" x for x in M.dyidydn], " "))
        println(io, @sprintf "%d" M.mz)
        println(io, join([@sprintf "%e" x for x in M.z], " "))
        println(io, join([@sprintf "%e" x for x in M.zmdn], " "))
        println(io, join([@sprintf "%e" x for x in M.dzidzup], " "))
        println(io, join([@sprintf "%e" x for x in M.dzidzdn], " "))
    end
end


function make_uniform_axes(
    mesh  ::BifrostMesh,
    new_mx::Integer,
    new_my::Integer,
    new_mz::Integer,
    )
    # Get new mesh-axes
    new_x = collect(LinRange(mesh.x[1], mesh.x[end], new_mx))
    new_y = collect(LinRange(mesh.y[1], mesh.y[end], new_my))
    new_z = collect(LinRange(mesh.z[1], mesh.z[end], new_mz))
    return new_x, new_y, new_z
end
