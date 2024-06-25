# BifrostTools.jl

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ITA-Solar.github.io/BifrostTools.jl/dev/)
[![Build Status](https://github.com/ITA-Solar/BifrostTools.jl/actions/workflows/CI.yml/badge.svg?branch=develop)](https://github.com/ITA-Solar/BifrostTools.jl/actions/workflows/CI.yml?query=branch%3Adevelop)
[![Coverage](https://codecov.io/gh/ITA-Solar/BifrostTools.jl/branch/develop/graph/badge.svg)](https://codecov.io/gh/ITA-Solar/BifrostTools.jl)

Tools for reading and working with simulation output from [The stellar atmosphere simulation code Bifrost](https://ui.adsabs.harvard.edu/abs/2011A%26A...531A.154G/abstract) in Julia.

This Julia package is created for working with *Bifrost* data **efficiently**. 
Load single or multiple simulation snapshots, and analyse data with Julia speed. 

This package is an extension of `Bifrost.jl`, a script written by Mikolaj Szydlarski.

## Documentation
The documentation is available at [https://ita-solar.github.io/BifrostTools.jl](https://ita-solar.github.io/BifrostTools.jl)

## Quick user guide
To load the package, type the following in the REPL

```{julia}
using BifrostTools
```

### Using `get_var`
The function `get_var` is the main function for reading data from Bifrost simulations.
It can read single or multiple snapshots, and it can read full data cubes or slices.
It can read primary variables or auxiliary variables.

The command

```{julia}
variable = get_var(expname, snap, expdir, variable)
```
loads the (primary or auxiliary) variable `variable` from snapshot `snap` in the simulation `expname` located in the directory `expdir`.

By creating a `BifrostExperiment` object 

```{julia}
brxp = BifrostExperiment(expname, expdir)
```
we can access the mesh file

```{julia}
brxp.mesh
```
snapshot numbers
```{julia}
brxp.snaps
```
and the calling signature of `get_var` can be simplified

```{julia}
variable = get_var(brxp, snap, variable)
```

Using optional keyword-arguments in `get_var` allows us to convert units, destagger variables, rotate the grid, and read slices of the full cube.

The command

```{julia}
bx = get_var(brxp, snap, "bx"; units="si", destagger=true)
```
will load the $x$-component of the magnetic field in SI units and destagger it to the cell center.

See the [documentation](https://ITA-Solar.github.io/BifrostTools.jl/dev/) for further information and more elaborate example usage.
