# Example usage

To load the package, type the following in the REPL

```{julia}
using BifrostTools
```

## Using `BifrostExperiment`

In this example, we look at the simulation *cb24oi*.
we start with defining the part to the simulation directory and the name of the simulation.

```{julia}
expdir = "/mn/stornext/d21/RoCS/matsc/3d/run/cb24oi/"
expname = "cb24oi"
```

These variables can be passed to the `BifrostExperiment` structure, which creates an instance that lets us access the mesh file, snapshot numbers, and so on.

```{julia}
xp = BifrostExperiment(expname,expdir)

# Mesh file that holds grid info etc.
mesh = xp.mesh

x = mesh.x
y = mesh.y
z = mesh.z

# vector with snap numbers
snaps = xp.snaps
```

## Reading data from the simulation with `get_var`

Reading data from the snapshot or aux files is handled through the `get_var` function. Due to `Julia`'s multiple dispatch functionality, there are several ways to call the `get_var` function, but we recommend using the following function to simplify the calling signature:

```{julia}
get_var(
    xp::BifrostExperiment,
    snap::Union{<:Integer, AbstractVector{<:Integer}},
    variable::String,
    args...
    ;
    kwargs...
)
```

This funciton loads a *variable* from one or multiple snapshots of `xp`. If multiple snapshots are loaded, the variable is returned as a `Vector` of snapshots. The available variables are

The primary variables:

- "r":  density
- "px": x-component of momentum
- "py": y-component of momentum
- "pz": z-component of momentum
- "e":  energy

if `params["do_mhd"]` is `true` in the params file, we also have the magnetic field:

- "bx": x-component of magnetic field
- "by": y-component of magnetic field
- "bz": z-component of magnetic field

and auxilliary variables (variables in params["aux"]):

- "p": pressure
- "tg": gas temperature
- Q terms

The following are optional keyword-arguments

- `units::String`: Converts variables to "si" or "cgs" units. `units="si"` or `units="cgs"`.
- `destagger::Bool`: Performs 5th-order interpolation to center the variable, and should only be used when reading variables that are staggered, (e.g. velocity or magnetic field). The function uses the default direction (destaggeroperation) associated with the variable unless otherwise stated by the `destaggeroperation` keyword.
- `destaggeroperation::String`: Determines which direction to destagger the variable. This is by default handled automatically by the `destaggeroperation` dictionary. 
- `rotate_about_x::String`: Rotate coordinate system to the normal "right hand system" with *z*-axis pointing upwards by passing `rotate_about="x"`
- `slicex::AbstractVector{<:Integer}`: Load a slice or slices in x-axis. Give e.g. `slicex=[32, 410]` or `slicex=40:90`
- `slicey::AbstractVector{<:Integer}`: Load a slice or slices in y-axis
- `slicez::AbstractVector{<:Integer}`: Load a slice or slices in z-axis
- `squeeze::Bool`: Removes singleton dimensions (dimensions with length 1)

### Loading a single snapshot

With `xp` as defined above, we define a snapshot that we want to investigate. When loading the full cube in code units, the variables are memory mapped, making them fast to load.

```{julia}
snap = 700
# Load some quantities for the full cube in code units
pressure = get_var(xp, snap, "p")
density = get_var(xp, snap, "r")
temperature = get_var(xp, snap, "tg")
```

### Converting units

If we want *si* or *cgs* units:
```{julia}
snap = 700
# Load some quantities for the full cube in si or cgs units
pressure = get_var(xp, snap, "p", units="cgs")
rho = get_var(xp, snap, "r", units="si")
# The temperature is written in Kelvin from before
temperature = get_var(xp, snap, "tg")
```

### Reading a slice of the full cube

If we're only interested in a small part of the cube, we can use the slicing functionality of `get_var`. Use the `squeeze` keyword to drop singleton dimensions.

We can load only the surface

```{julia}
idz = argmin(abs.(mesh.z))

rho = get_var(xp, snap, "r"; units="si", slicez=[idz], squeeze=true)
temperature = get_var(xp, snap, "tg"; units="si", slicez=[idz], squeeze=true)
```

or a smaller cube around the surface

```{julia}
rho = get_var(xp, snap, "r",
     units="si", slicex=100:200, slicey=400:500, slicez=[idz-20:idz+20])
temperature = get_var(xp, snap, "tg",
     units="si", slicex=100:200, slicey=400:500, slicez=[idz-20:idz+20])
```

### Interpolating staggered variables to the grid center

Interpolating staggered variables (destaggering) can be handled through `get_var`. This is recommended because `get_var` can determine the interpolation direction, and if you want to slice a variable, it takes care of the correct ordering of interpolating and slicing. 

```{julia}
# Read and destagger vertical momentum in si units
pz = get_var(xp, snap, "pz", units="si", destagger=true)
```

### Reading multiple snapshots

If you want the time evolution 

If you want to get the time evolution of a quantity, you can simply pass a vector of snapshots. The `get_var` function uses Julia's threads functionality to read multiple snapshots in parallel. 

```{julia}
snaps = 100:150

rho = get_var(xp, snaps, "r", units="si", slicez=[idz], squeeze=true)
# Calculate vertical velocity
vz = get_var(xp, snaps, "pz", units="si", slicez=[idz], squeeze=true) ./ rho
```

### Rotating grid

Rotate about the x-axis to get vector quantities on a *z*-axis that points upwards.

```
# When we rotate about the x-axis, this is what happens to the grid
z = -mesh.z
y = -mesh.y

# Load x-component of B-field and rotate about x-axis
bx = get_var(xp, isnap, "bx", units="si", destagger=true, rotate_about="x")
```

## Loading the simulation parameters with `read_params`

The `read_params` function reads the params file. It can be called by giving the full filename, like the following

```{julia}
snap = 500
params_file = joinpath(expdir,string(expname,"_",snap,".idl"))
params = read_params(params_file)
```

or

```{julia}
read_params(expname,snap,expdir)
```