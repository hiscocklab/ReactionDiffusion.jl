# API for ReactionDiffusion.jl

ReactionDiffusion aims to be an easy-to-use *and* computationally-efficient pipeline to simulate biologically-inspired reaction-diffusion models. It is our hope that models can be built with just a few lines of code and solved without the user having any knowledge of PDE solver methods. 

This is achieved by drawing from a range of numerical routines from the SciML packages, including [Catalyst.jl](https://github.com/SciML/Catalyst.jl), [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl), [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl), and [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl). Currently, there are limited options to customize the simulation method (e.g., currently only 1D simulations are supported with reflective boundary conditions). 


## Data structures
We introduce two structures

- `save_turing` is an object that records parameter sets that undergo a Turing instability. It has the following fields:

    - `steady_state_values`: The computed steady state values of each variable
    - `reaction_params`: The reaction parameters
    - `diffusion_constants`: The diffusion constants
    - `initial_conditions`: The initial conditions of each variable used to compute the steady state values
    - `pattern_phase: The predicted phase of each of the variables in the final pattern. `[1 1]` would be in-phase, `[1 -1]` would be out-of-phase
    - `wavelength`: The wavelength that is maximally unstable
    - `max_real_eigval`: The maximum real eigenvalue associated with the Turing instability
    - `non_oscillatory`: If `true`, this parameter set represents a stationary Turing pattern. If `false`, the unstable mode has complex eigenvalues and thus may be oscillatory.

- `model_parameters` is an object that records parameter values for `model` simulations. It has the following fields:

    -`reaction`: reaction parameters
    -`diffusion`: diffusion constants
    -`initial_condition`: values of each of the variables used for homogeneous initial conditions
    -`initial_noise`: level of (normally distributed) noise added to the homogeneous initial conditions
    -`domain_size`: size of 1D domain
    -`random_seed`: seed associated with random number generation; only set this when you need to reproduce exact simulations each run


## Functions

```@docs
screen_values
get_params
get_param
returnTuringParams
simulate
```



