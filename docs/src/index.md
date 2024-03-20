# BifrostTools.jl

BifrostTools is a Julia package for working with Bifrost simulations.

The goal of this package is to analyse *Bifrost* data **faster** than if you were to use python. 

BifrostTools lets you can load variables from single or multiple simulation snapshots, destaggers data automatically, and calculates derivatives. It is made to analyse data efficiently. 


```@docs
get_var(
        xp::BifrostExperiment,
        snap::Union{<:Integer, AbstractVector{<:Integer}},
        variable::String,
        args...
        ;
        kwargs...
        )
```