module BifrostTools

using FortranFiles: FortranFile, read, readlines
using OffsetArrays
using DelimitedFiles
using Printf
using Interpolations
using Mmap
using LoopVectorization

include("mesh.jl")
include("utils.jl")
include("stagger_operators.jl")
include("experiment.jl")
include("read_params_snap_aux.jl")
include("write_params_snap_aux.jl")
include("eos_tables.jl")
include("unit_conversion.jl")

#-------------------------------------------------------------------------------
# Exports
#-------------------------------------------------------------------------------

# Structs
export BifrostMesh
export BifrostExperiment
export EOSTables

# mesh.jl
export make_uniform_axes
export get_axes
export mesh2file

# utils.jl
export change_snap_resolution, duplicate_xz_plane

# read_params_snap_aux.jl
export read_params
export get_var, get_snap_numbers, get_electron_density

# eos_tables.jl
export get_eostable
export eos_interpolate


# stagger_operators.jl
# Basic stagger operations with optional BC extrapolation
export up
export dup

export dn
export ddn 

export xup
export dxup

export xdn
export dxdn

export yup
export dyup

export ydn
export dydn

export zup
export dzup

export zdn
export dzdn

export destaggeroperation

end # module BifrostTools
