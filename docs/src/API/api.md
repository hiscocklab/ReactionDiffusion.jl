# API for ReactionDiffusion.jl

ReactionDiffusion aims to be an easy-to-use *and* computationally-efficient pipeline to simulate biologically-inspired reaction-diffusion models. It is our hope that models can be built with just a few lines of code and solved without the user having any knowledge of PDE solver methods. 

This is achieved by drawing from a range of numerical routines from the SciML packages, including [Catalyst.jl](https://github.com/SciML/Catalyst.jl), [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl), [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl), and [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl). Currently, there are limited options to customize the simulation method (e.g., currently only 1D simulations are supported with reflective boundary conditions). 


## Data structures

## Functions

```@docs
model_parameters
save_turing
screen_values
get_params
get_param
returnTuringParams
simulate
```



