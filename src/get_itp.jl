"""
    get_itp(brxp::BifrostExperiment, snaps, variables, itp_type)

Return an Interpolations.jl-object for the Bifrost-experiment. Also return
metadata for the interpolation.

The interpolation can be in both space and time, and on both scalar and
vector fields.

# Arguments
- `brxp`: The `BifrostExperiment` to interpolate.
- `snaps`: The Bifrost snapshot(s) to interpolate. An `AbstractVector` of
    multiple snapshots yields interpolation in time (as if it was an additional
    spatial dimension).
- `variables`: The variables to interpolate. Multiple variables yields vector
    field interpolation (returning static arrays). The variables can be given as
    a space-separated string or as a vector of `String`s/`Symbol`s.
- `itp_type`: The interpolation type.

    Can be any `Interpolations.InterpolationType` or any of the `Symbol`s
    `:constant`, `:linear`, `:quadratic_spline`, `:cubic_spline` (or their
    `String` equivalent).

    It can also be a `Vector`/`Tuple` for applying different interpolation
    types to different axes.

    When using `:constant` or `:linear`, the interpolation uses `BSpline`s if
    all axes are uniform. It uses `Gridded` interpolation otherwise. Non-linear
    interpolation using `:quadratic_spline` and `:cubic_spline` requires
    uniform axes.

# Keyword arguments

- `itp_bc=Throw()`: Extrapolation boundary conditions. Can be any
    `Interpolations.Boundarycondition` or the `Symbol`s (or `String`s);
    `:periodic`, `:flat`, `:throw`, `:reflect`, `:line`.
    It can also be a vector of boundary conditions as long as the length is
    larger than 1 and equal to the number of dimensions.
- `units="code"`: E.g. `"si"` or `"cgs"`.
- `slicex/y/z=[]`: Ranges for which to slice the axes of the snapshot(s).
- `saveitp=false`: Whether to save the interpolator to a JLD2-file or not.
- `name="no_name"`: of the interpolation-object (and possibly the saved file).
- `dir="./"`: The directory where the interpolation object should be saved.
- `precision=Float32`: The floating point precision of the Bifrost-snapshots.

# Examples
```julia-repl
julia> using BifrostTools, JLD2

julia> xp = BifrostExperiment("fanspine_exhaust64", "test/fanspine_exhaust64");

julia> itp, metadata = get_itp(xp, 1:5, "ey", :cubic_spline);

julia> itp(18f0, 2.5f0, -5f0, 6.002f0) # itp(x, y, z, t)
0.00053988147f0

julia> xp = BifrostExperiment("en48", "test/sp.n064");

julia> get_itp(xp, 1, "bx by bz", [:linear, :constant];
           slicex=24:24,
           itp_bc=[:periodic, :throw],
           saveitp=true,
           name="test"
       );

julia> @load "test.itp.jld2" itp metadata

julia> itp(1f1, -1f0) # itp(y, z)
3-element StaticArraysCore.SVector{3, Float32} with indices SOneTo(3):
 0.70152
 0.00021198636
 0.7067044
```
"""
function get_itp(
    brxp::BifrostExperiment,
    snaps::Union{Integer,AbstractVector},
    variables,
    itp_type;
    itp_bc=Throw(),
    units::String="code",
    slicex::AbstractVector{<:Integer}=Int[],
    slicey::AbstractVector{<:Integer}=Int[],
    slicez::AbstractVector{<:Integer}=Int[],
    saveitp=false,
    name="no_name",
    dir="./",
    precision::DataType=Float32,
)
    # Put the snapshot(s) of each variable in a vector.
    # If `variable` is a String, not a vector of strings, `vectorfield` will
    # not be a vector, just the snapshot(s) itself.
    variables = parse_variables(variables)
    fetched_vars = map(variables) do var
        data = get_var(
            brxp, snaps, var;
            units=units,
            destagger=isstaggered(var),
            slicex=slicex,
            slicey=slicey,
            slicez=slicez,
            precision=precision,
            squeezeresult=true # dimensions of length 1 is not supported
        )
        # `data` is a vector if we fetched multiple snapshots. We want time to be
        # a dimension in addition to the spatial ones, so we stack the data.
        if data isa Vector
            return stack(data)
        else
            return data
        end
    end

    # Get spatial axes
    spatial_axes = get_axes(
        brxp,
        units=units,
        slicex=slicex,
        slicey=slicey,
        slicez=slicez
    )
    # Get temporal axis
    temporal_axis = precision.(get_var(brxp, snaps, "t"; units=units))
    axes = (spatial_axes..., temporal_axis)

    # A string of the call signature. Used in the metadata.
    # The signature must be created before squeezing.
    call_signature = get_callsignature(axes)

    # Squeeze because Interpolation.jl does not allow dimensions of length 1.
    # Make uniform axes ranges because Interpolations.scale does not allow
    # vectors.
    axes = make_uniform_axes_ranges(squeeze(axes))

    # Resolve the the interpolation type and boundary conditions from the
    # arguments. May throw `ArgumentError`.
    itp_type, gridded_or_bspline = resolve_itp_type(itp_type, axes)
    itp_bc = resolve_itp_bc(itp_bc, length(axes))

    # From the fetched variables, axes and itp-type, crate the interpolation
    # object using interpolations.jl
    itp = create_interpolator(
        fetched_vars,
        itp_type,
        itp_bc,
        axes,
        gridded_or_bspline
    )

    metadata = Dict(
        :itp_name => name,
        :experiment_name => brxp.expname,
        :snaps => snaps,
        :vars => variables,
        :units => units,
        :size => Tuple(length(ax) for ax in axes),
        :bounds => [(first(ax), last(ax)) for ax in axes], 
        :call_signature => call_signature,
        :itp_type => string(itp_type),
        :itp_bc => string(itp_bc),
        :uniform_axes => isuniform.(axes),
        :julia => VERSION,
        :BifrostTools => pkgversion(@__MODULE__),
        :Interpolations => pkgversion(Interpolations),
        :StaticArrays => pkgversion(StaticArrays),
        :JLD2 => pkgversion(JLD2)
    )

    if saveitp
        jldsave(joinpath(dir, name*".itp.jld2"); itp, metadata)
    end
    return (; itp=itp, metadata=metadata)
end

"""
    isuniform(x; tolerance)
Return `true` if the axis `x` is regarded as uniform.

# Method
1. Calculate the average grid spacing of the axis, `dx`.
2. Use `dx` to create a uniform version of the axis.
3. Find the maximum deviation between the original axis and its uniform version.
4. If the maximum deviation is less than `tolerance*dx`, return `true`.
   Return `false` otherwise.
"""
function isuniform(x; tolerance=1e-3)
    n = length(x)
    x_f64 = Float64.(x)
    dx = (x_f64[end] - x_f64[1]) / (n - 1) 
    uniform_x = [x_f64[1] + (i - 1)*dx for i in 1:n]
    maxdeviation = maximum(abs(x_f64[i] - uniform_x[i]) for i in 1:n)
    return maxdeviation / abs(dx) <= tolerance 
end

"""
    make_uniform_axes_ranges(axes)
Create ranges from uniform axes so that we can use Interpolations.scale
"""
function make_uniform_axes_ranges(axes)
    axes_are_uniform = all(isuniform.(axes))
    if axes_are_uniform
        map(axes) do ax
            if !(ax isa AbstractRange)
                return range(ax[1], ax[end]; length=length(ax))
            else
                return ax
            end
        end
    else
        return axes
    end
end

"""
    resolve_itp_type(itp_type, axes::Tuple)
Given the axes on which the interpolation is acting, resolve which
`Interpolation.InterpolationType` that `itp_type` corresponds to.

`itp_type` can be a string, symbol, an actual `Interpolation.InterpolationType`,
or a Vector/Tuple of them. E.g `"linear"`, `:linear`, `Gridded(Linear())`,  or
`[:cubic_spline, BSpline(Constant())]`.

For `BSpline` interpolation to work, the axes need to be ranges, not vectors.
If they are ranges, `BSpline` will be prefered over `Gridded` if `itp_type` is
given as symbols or strings.

Passing `itp_type` as a mix of `Gridded` and `BSpline` will throw an error.
"""
function resolve_itp_type(itp_type, axes::Tuple)
    # Check whether the axes are vectors or ranges.
    # B-spline interpolation requires all axes to be ranges.
    axes_are_ranges = all(ax isa AbstractRange for ax in axes)

    if itp_type isa Union{Tuple, Vector}
        # Check that the length matches the number of dimensions.
        n = length(itp_type)
        if n != length(axes) || n == 1
            throw(ArgumentError("""itp_type is the wrong length.

The length should be equal to the number of dimensions (including time). If
you want to apply one interpolation type to all dimensions, use a scalar, not a
vector or tuple.
"""
            ))
        end

        # Resolve the interpolation type based on the arguments and whether the
        # axis is a range or not.
        itp_type = Tuple(map(
            match_itp_type,
            itp_type,
            fill(axes_are_ranges, length(itp_type))
        ))
        # Each axis can have a different interpolation order, but they must have
        # the same type, e.g. either Gridded or Bspline.
        if all(@. itp_type isa Gridded)
            return itp_type, :gridded
        elseif all(@. itp_type isa BSpline)
            return itp_type, :bspline
        else
            throw(ArgumentError("Cannot mix BSpline and Gridded interpolation"))
        end
    else
        itp_type = match_itp_type(itp_type, axes_are_ranges)  
        if itp_type isa Gridded
            return itp_type, :gridded
        else
            return itp_type, :bspline
        end
    end
end

"""
    match_itp_type(it, isrange)
Return an appropriate Interpolation.InterpolationType based on
the given `it` and whether its axis is a range or not.
"""
function match_itp_type(it, isrange::Bool)
    if it isa Union{Symbol, String}
        it = Symbol(it)
        if isrange
            if it == :constant             return BSpline(Constant())
            elseif it == :linear           return BSpline(Linear())
            elseif it == :quadratic_spline return BSpline(Quadratic())
            elseif it == :cubic_spline     return BSpline(Cubic())
            else return unknown_itp_type_error()
            end
        else
            if it == :constant             return Gridded(Constant())
            elseif it == :linear           return Gridded(Linear())
            elseif it == :quadratic_spline return nonlinear_nonuniform_error()
            elseif it == :cubic_spline     return nonlinear_nonuniform_error()
            else return unknown_itp_type_error()
            end
        end
    elseif it isa Interpolations.InterpolationType return it
    else return unknown_itp_type_error()
    end
end

function nonlinear_nonuniform_error()
        throw(ArgumentError("""Non-linear interpolation on non-uniform axis:
This is currently not implemented. However, it is possible. For example, by
transforming the coordinates to the uniform index-space before interpolation.
"""
        ))
end

function unknown_itp_type_error()
    throw(ArgumentError(
        "Uknown itp_type. Try :constant, :linear or :cubic_spline."
    ))
end

"""
    resolve_itp_bc(itp_bc, num_dims)
Resolve and validate `itp_type`.

If it is a `Vector` or a `Tuple`, check that the `num_dims` matches the number of
boundary conditions.
"""
function resolve_itp_bc(itp_bc, num_dims)
    if itp_bc isa Union{Tuple, Vector}
        # Check that the number of boundary conditions matches the number of
        # dimensions.
        n = length(itp_bc)
        if n != num_dims || n == 1
            throw(ArgumentError("""itp_bc is the wrong length.

The length should be equal to the number of dimensions (including time). If
you want to apply one boundary condition to all dimensions, use a scalar, not a
vector or tuple.
"""
            ))
        end
        return Tuple(map(match_itp_bc, itp_bc))
    else return match_itp_bc(itp_bc)
    end
end

"""
    match_itp_bc(bc)
Match `bc` to corresponding `Interpolations.BoundaryCondition` if `bc`
is a `String` or a `Symbol`. Return `bc` if it already is an
`Interpolations.BoundaryCondition`.
"""
function match_itp_bc(bc)
    # If the given itp_bc is a symbol or string, map it to the
    # corresponding boundary-condtition-function from Interpolations.jl.
    # If the argument is already given as a Interpolations.BoundaryCondition
    # type, use that.
    if bc isa Union{Symbol, String}
        bc = Symbol(bc)
        if bc == :periodic  return Periodic()
        elseif bc == :flat  return Flat()
        elseif bc == :throw return Throw()
        elseif bc == :line  return Line()
        elseif bc == :reflect return Reflect()
        else return unknown_itp_bc_error()
        end
    elseif bc isa Interpolations.BoundaryCondition return bc
    else return unknown_itp_bc_error()
    end
end

function unknown_itp_bc_error()
    throw(ArgumentError("""Uknown itp_bc or unknown formatting. Try something like:

"flat"
Flat()
["flat", "periodic"]
(Flat(), Periodic()).
"""))
end



"""
    create_interpolator(field, itp_type, itp_bc, axes, uniform_grid)

Return an callable interpolation object for structured grids using
Interpolations.jl.

# Arguments
- `field`: The field quantity to interpolate. It should have the same number
    of dimensions as `axes`.
- `itp_type`: The type of interpolation.
- `itp_bc`: Boundary conditions/extrapolation.
- `axes`: The axes of the field.
- `gridded_or_bspline`: Whether we are creating a B-spline or Gridded
    interpolation.
"""
function create_interpolator(
    fetched_vars,
    itp_type,
    itp_bc,
    axes,
    gridded_or_bspline
)
    n_components = length(fetched_vars)

    # Are we interpolating a vector field?
    vectorfield = n_components > 1

    # If there is only one component, `fetched_vars` may be a
    # single-element vector. Discard the Vector-wrapping in this case,
    # because we want fetched_vars to be an array.
    if fetched_vars isa Vector && n_components == 1
        fetched_vars = first(fetched_vars)
    end

    # Use an array of static vectors if we are interpolating a vector field
    field = vectorfield ?
        to_array_of_svectors(Val(n_components), fetched_vars) :
        fetched_vars

    # If we are creating a BSpline interpolation we have to `scale`.
    if gridded_or_bspline == :bspline
        return extrapolate(
            scale(interpolate(field, itp_type), axes...),
            itp_bc
        )
    elseif gridded_or_bspline == :gridded
        return extrapolate(interpolate(axes, field, itp_type), itp_bc)
    else
        throw(ArgumentError("Unknown `gridded_or_bspline`-value."))
    end
end

"""
    to_array_of_svectors(::Val{N}, vec) where {N}
Turn a `vec`tor of arrays into an array of `SVector`s. Assumes all elements
of `vec` have the same shape, size and type.

# Example
```julia-repl
julia> vec = [bitrand(2,2) for _ in 1:3]
3-element Vector{BitMatrix}:
 [1 1; 0 0]
 [1 1; 0 1]
 [1 1; 0 1]

julia> to_array_of_svectors(Val(3), vec)
2×2 Matrix{StaticArraysCore.SVector{3, Bool}}:
 [1, 1, 1]  [1, 1, 1]
 [0, 0, 0]  [0, 1, 1]
```
"""
function to_array_of_svectors(::Val{N}, vec) where {N}
    A = first(vec)
    out = similar(A,SVector{N,(eltype(A))})
    @inbounds for i in eachindex(A)
        out[i] = SVector{N}(ntuple(k -> vec[k][i], Val(N)))
    end
    return out
end

"""
    get_callsignature(axes)
Return the physical symbol each variable in the axes represents. Axes of
length one is squeezed away when creating an interpolation object. Assumes
`spatial_axes` is either of
- (xaxis)
- (xaxis, yaxis)
- (xaxis, yaxis, zaxis)
- (xaxis, yaxis, zaxis, time)
"""
function get_callsignature(axes)
    axis_symbols = ["x", "y", "z", "t"]
    call_signature = [
        axis_symbols[i] for i in eachindex(axes) if length(axes[i]) > 1
    ]
    return "("*join(call_signature, ", ")*")"
end
