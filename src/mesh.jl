"""
    BifrostMesh

Stores Bifrost grid information in struct
"""
mutable struct BifrostMesh
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
  
    function BifrostMesh(mesh_file::String)
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
  
    function BifrostMesh(
      mx::Int64,
      dx::Float64,
      my::Int64,
      dy::Float64,
      mz::Int64,
      dz::Float64,
    )
      # Create an uniform Mesh
      if (dx == 0.0)
        dx = 1.0
      end
      if (dy == 0.0)
        dy = 1.0
      end
      if (dz == 0.0)
        dz = 1.0
      end
  
      # x
      x = collect(0:mx-1) .* dx
      xmdn = x .- (0.5 * dx)
      dxidxup = zeros(mx) .+ (1.0 / dx)
      dxidxdn = zeros(mx) .+ (1.0 / dx)
  
      # y
      y = collect(0:my-1) .* dy
      ymdn = y .- (0.5 * dy)
      dyidyup = zeros(my) .+ (1.0 / dy)
      dyidydn = zeros(my) .+ (1.0 / dy)
  
      # z
      z = collect(0:mz-1) .* dz
      zmdn = y .- (0.5 * dz)
      dzidzup = zeros(mz) .+ (1.0 / dz)
      dzidzdn = zeros(mz) .+ (1.0 / dz)
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
        mx * my * mz,
      )
    end
    function BifrostMesh(par::Dict{String,Any})
      if (par["meshfile"] == "nomesh.dat")
        return BifrostMesh(par["mx"], par["dx"], par["my"], par["dy"], par["mz"], par["dz"])
      else
        return BifrostMesh(par["meshfile"])
      end
    end
  
    function BifrostMesh(x::Vector{T},y::Vector{T},z::Vector{T}) where {T<:AbstractFloat}
      
      mx = size(x)[1]
      xd = br_dn(x)
      xddn = 1.0 ./ br_ddn(x)
      xdup = 1.0 ./ br_dup(x)
  
      if (mx == 1)
        xddn[:] = 1.0
        xdup[:] = 1.0
      end 
  
      my = size(y)[1]
      yd = br_dn(y)
      yddn = 1.0 ./ br_ddn(y)
      ydup = 1.0 ./ br_dup(y)
  
      if (my == 1)
        yddn[:] = 1.0
        ydup[:] = 1.0
      end
  
      mz = size(z)[1]
      zd = br_dn(z)
      zddn = 1.0 ./ br_ddn(z)
      zdup = 1.0 ./ br_dup(z)
  
      if (mz == 1)
        zddn[:] = 1.0
        zdup[:] = 1.0
      end
      new(
        mx,
        x,
        xd,
        xdup,
        xddn,
        my,
        y,
        yd,
        ydup,
        yddn,
        mz,
        z,
        zd,
        zdup,
        zddn,
        mx * my * mz,
      )
    end
  end