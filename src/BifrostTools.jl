module BifrostTools

using FortranFiles: FortranFile, read, readlines
using OffsetArrays
using DelimitedFiles
using Printf
using Interpolations
using Mmap

include("mesh.jl")
include("utils.jl")
include("stagger_operators.jl")
include("experiment.jl")
include("div_operators.jl")
include("read_params_snap_aux.jl")
include("write_params_snap_aux.jl")
include("eos_tables.jl")
include("unit_conversion.jl")

# Exports
export squeeze
export read_params
export load_snapdata
export load_auxdata
export load_snapvariable
export load_auxvariable
export get_var, get_staggered_var, get_snap_numbers, get_electron_density
export change_snap_resolution, duplicate_xz_plane

export make_uniform_axes

export arr_ffile
export heatmap_xy
export heatmap_xz
export heatmap_yz

export get_eostable
export get_expieos_err

export get_epstable
export get_temtable
export get_opatable

export get_ne_epstable
export get_ne_temtable
export get_ne_opatable

export get_lndlnT_table
export get_theta_rho_table

export eos_interpolate
export mesh2file 
export fix_mesh 

export BifrostMesh
export BifrostExperiment
export EOS_Tables

# debugging

export spitzer_debug_file

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

# fast stagger operations using circular shift

export cdxup
export cdxdn

export cdyup
export cdydn

export cdzup
export cdzdn

export cdivup
export cdivdn

export cgrad_dn
export cgrad_up

export claplace_du
export claplace_ud

export cdivB_clean_du

# special functions

export divup
export divdn

export gradup
export graddn

export laplacedu
export laplaceud

export poissondu
export poissonud

export divB_clean_ud
export divB_clean_du



end # module BifrostTools
