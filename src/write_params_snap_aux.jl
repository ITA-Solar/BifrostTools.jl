
"""
    change_snap_resolution(
        xp      ::BifrostExperiment,
        isnap   ::Integer,
        new_x   ::Vector{<:Real},
        new_y   ::Vector{<:Real},
        new_z   ::Vector{<:Real},
        itp_bc  ::Interpolations.BoundaryCondition,
        ;
        filename::String="out.snap"
        )
Writes a new snapfile with a new resolution, according to new axes,
given as arguments. To fill the values of the new snap this function uses
gridded linear interpolation (from Interpolations.jl) of the old snap, given
by the BifrostExperiment and isnap argument.

Boundary conditions could be either Throw(), Flat(), Line(), Periodic()
or Reflect().
"""
function change_snap_resolution(
    xp      ::BifrostExperiment,
    isnap   ::Integer,
    new_x   ::Vector{<:Real},
    new_y   ::Vector{<:Real},
    new_z   ::Vector{<:Real},
    itp_bc  ::Interpolations.BoundaryCondition,
    ;
    filename::String="out.snap"
    )
    primaries, params = load_snapdata(xp.expname, isnap, xp.expdir)
    snapsize, numvars, _ = get_snapsize_and_numvars(params)
    change_snap_resolution(primaries, xp.mesh, numvars, 
                                new_x, new_y, new_z, itp_bc;
                                filename
                                )
end

function change_snap_resolution(
    snapname::String,
    meshname::String,
    new_x   ::Vector{<:Real},
    new_y   ::Vector{<:Real},
    new_z   ::Vector{<:Real},
    itp_bc  ::Interpolations.BoundaryCondition,
    ;
    filename::String="out.snap"
    )
    params_filename = snapname[1:end-4]*"idl"
    params = read_params(params_filename)
    primaries = load_snapdata(snapname, params)
    snapsize, numvars, _ = get_snapsize_and_numvars(params)
    mesh = BifrostMesh(meshname)
    change_snap_resolution(primaries, mesh, numvars,
                                new_x, new_y, new_z, itp_bc;
                                filename
                                )
end 

function change_snap_resolution(
    primaries::Array{<:Real, 4},
    mesh     ::BifrostMesh,
    numvars  ::Int64,
    new_x    ::Vector{<:Real},
    new_y    ::Vector{<:Real},
    new_z    ::Vector{<:Real},
    itp_bc   ::Interpolations.BoundaryCondition,
    ;
    filename::String="out.snap"
    )
    # Interpolation through the package Interpolations.jl. Currently this
    # function only uses gridded linear interpolation because the BSpline 
    # interoplation in Interpolations.jl requires the axes to be 
    # AbstractRange's, and I have not yet found a workaround for the non-
    # uniform Bifrost messh.
    #
    # However, one may choose from the boundary conditions available 
    # through Interpolations.jl e.g. Periodic() or Flat().
    itp_type = Gridded(Linear())

    # Interpolations don't like meshes where axes have length 1.
    # So we check if we are dealing with a 2D snap or not.
    if length(mesh.x) == 1
        dims = "yz"
    elseif length(mesh.y) == 1
        dims = "xz"
    elseif length(mesh.z) == 1
        dims = "xy"
    else
        dims = "3D"
    end

    # Construct an interpolation object for each primary.
    itp = Array{AbstractInterpolation}(undef, numvars)
    get_primary_interpolators!(itp, mesh.x, mesh.y, mesh.z,
                                numvars, primaries, itp_type, itp_bc,
                                dims
                                )
                        
    # Construct array to hold interpolated primaries
    nx, ny, nz = (length(new_x), length(new_y), length(new_z))
    working_precision = typeof(primaries[1])
    new_primaries = Array{working_precision}(undef, nx, ny, nz, numvars)

    # fill in interpolated values
    for var = 1:numvars
        if dims == "3D"
            for i = 1:nx
                for j = 1:ny
                    for k = 1:nz
                        x, y, z = new_x[i], new_y[j], new_z[k]
                        new_primaries[i,j,k,var] = itp[var](x,y,z)
                    end
                end
            end
        elseif dims == "yz"
            for j = 1:ny
                for k = 1:nz
                    y, z = new_y[j], new_z[k]
                    new_primaries[1,j,k,var] = itp[var](y, z)
                end
            end
        elseif dims == "xz"
            for i = 1:nx
                for k = 1:nz
                    x, z = new_x[i], new_z[k]
                    new_primaries[i,1,k,var] = itp[var](x, z)
                end
            end
        elseif dims == "xy"
            for i = 1:nx
                for j = 1:nj
                    x, y = new_x[i], new_y[j]
                    new_primaries[i,j,1,var] = itp[var](x, y)
                end
            end
        end
    end

    println("Writing snap to $filename...")
    file = open(filename, "w+")
    write(file, new_primaries)
    close(file)
    return new_primaries, new_x, new_y, new_z
end 

function get_primary_interpolators!(
    itp      ::Array{AbstractInterpolation},
    x        ::Vector{<:Real},
    y        ::Vector{<:Real},
    z        ::Vector{<:Real},
    numvars  ::Integer,
    primaries::Array{<:Real, 4},
    itp_type ::Interpolations.InterpolationType,
    itp_bc   ::Interpolations.BoundaryCondition,
    dims     ::String
    )
    for var = 1:numvars
        if dims == "3D"
            itp[var] = interpolate((x, y, z),
                             primaries[:,:,:,var],
                             itp_type
                             )
        elseif dims == "xz"
            itp[var] = interpolate((x, z),
                                 primaries[:,1,:,var],
                                 itp_type
                                 )
        elseif dims == "xy"
            itp[var] = interpolate((x, y),
                                 primaries[:,:,1,var],
                                 itp_type
                                 )
        elseif dims == "yz"
            itp[var] = interpolate((y, z),
                                 primaries[1,:,:,var],
                                 itp_type
                                 )
        end
        itp[var] = extrapolate(itp[var], itp_bc)
    end
end


"""
    duplicate_xz_plane(
        xp      ::BifrostExperiment,
        isnap   ::Integer,
        ny      ::Integer,
        ;
        filename::String="out.snap"
        )
Duplicates the xz-plane into `ny` gridpoints in the y-axis. Used for creating
3D snapshots from a 2D xy-snapshot. 
"""
function duplicate_xz_plane(
    xp      ::BifrostExperiment,
    isnap   ::Integer,
    ny      ::Int64,
    ;
    filename::String="out.snap"
    )
    primaries, params = load_snapdata(xp.expname, isnap, xp.expdir)
    duplicate_xz_plane(primaries, ny; filename=filename)
end

function duplicate_xz_plane(
    snapname::String,
    ny      ::Int64,
    ;
    filename::String="out.snap"
    )
    params_filename = snapname[1:end-4]*"idl"
    params = read_params(params_filename)
    primaries = load_snapdata(snapname, params)
    duplicate_xz_plane(primaries, ny; filename=filename)
end 

function duplicate_xz_plane(
    primaries::Array{<:Real, 4},
    ny       ::Int64,
    ;
    filename::String="out.snap"
    )

    # Construct array to hold interpolated primaries
    nx, _, nz, numvars = size(primaries)
    working_precision = typeof(primaries[1])
    new_primaries = Array{working_precision}(undef, nx, ny, nz, numvars)

    # fill in values
    for var = 1:numvars
        for i = 1:nx
            for k = 1:nz
                new_primaries[i,:,k,var] .= primaries[i,1,k,var]
            end
        end
    end

    println("Writing snap to $filename...")
    file = open(filename, "w+")
    write(file, new_primaries)
    close(file)
    return new_primaries
end 
