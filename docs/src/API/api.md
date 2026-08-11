# API for ReactionDiffusion.jl

ReactionDiffusion aims to be an easy-to-use *and* computationally-efficient pipeline to simulate biologically-inspired reaction-diffusion models. It is our hope that models can be built with just a few lines of code and solved without the user having any knowledge of PDE solver methods. 

This is achieved by drawing on a range of SciML packages, including [Catalyst.jl](https://github.com/SciML/Catalyst.jl), [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl), [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl), and [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl). 

For users more familiar with PDE solvers, it is possible to specify optional arguments to e.g., control the method of discretisation or the specific solver algorithm used (see the API below). However, ReactionDiffusion does not aim to be a fully customizable PDE solving package that covers all bases; our focus is to make something that is easy-to-use *and* performant but only for the specific use case of 1D reaction-diffusion systems.

If you require more customization or a more flexible PDE-solving framework, we highly recommend [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) and/or [Catalyst.jl](https://github.com/SciML/Catalyst.jl).

# Functions
```@docs
Model
@reaction_network
@diffusion_system
parameter_set
simulate
PsuedoSpectralProblem
PseudoSpectralSolution
EnsembleProblem
solve
turing_wavelength
is_turing
filter_turing
timeseries_plot
interactive_plot
```



