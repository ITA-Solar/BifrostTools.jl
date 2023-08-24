module BifrostTools

using FortranFiles: FortranFile, read, readlines
using Plots: heatmap
using OffsetArrays
using DelimitedFiles
using Printf
using Interpolations: CubicSplineInterpolation, Line
using Base
using Mmap

include("utils.jl")
include("stagger_operators.jl")
include("mesh.jl")
include("experiment.jl")
include("div_operators.jl")
include("read_params_snap_aux.jl")
include("EOS_tables.jl")
include("plotting.jl")
include("unit_conversion.jl")

# Exports
export br_squeeze
export br_read_params
export br_load_snapdata
export br_load_auxdata
export br_load_snapvariable
export br_load_auxvariable
export get_var, get_staggered_var, get_snap_numbers, get_electron_density

export br_arr_ffile
export br_heatmap_xy
export br_heatmap_xz
export br_heatmap_yz

export br_get_eostable
export br_get_expieos_err

export br_get_epstable
export br_get_temtable
export br_get_opatable

export br_get_ne_epstable
export br_get_ne_temtable
export br_get_ne_opatable

export br_get_lndlnT_table
export br_get_theta_rho_table

export br_eos_interpolate
export br_mesh2file 
export br_fix_mesh 

export BifrostMesh
export BifrostExperiment
export EOS_tables

# debugging

export spitzer_debug_file

# Basic stagger operations with optional BC extrapolation

export br_up
export br_dup

export br_dn
export br_ddn 

export br_xup
export br_dxup

export br_xdn
export br_dxdn

export br_yup
export br_dyup

export br_ydn
export br_dydn

export br_zup
export br_dzup

export br_zdn
export br_dzdn

# fast stagger operations using circular shift

export br_cdxup
export br_cdxdn

export br_cdyup
export br_cdydn

export br_cdzup
export br_cdzdn

export br_cdivup
export br_cdivdn

export br_cgrad_dn
export br_cgrad_up

export br_claplace_du
export br_claplace_ud

export br_cdivB_clean_du

# special functions

export br_divup
export br_divdn

export br_gradup
export br_graddn

export br_laplacedu
export br_laplaceud

export br_poissondu
export br_poissonud

export br_divB_clean_ud
export br_divB_clean_du



end # module BifrostTools
