# Simulate pattern formation in 1D

Having now screened through parameter sets, we can now pick a single one and simulate the corresponding PDE. The simulation is performed on a 1D domain with reflective boundary conditions.  

We first select the parameter values from `turing_params`; here we choose the 1000th parameter set found. 

```julia
param1 = get_params(model, turing_params[1000])
```

Optionally, we can manually change the parameters too:

```julia
param1.reaction["n1"] = 10
```

The size of the 1D domain over which the simulation is performed is chosen according to the `turing_params.wavelength` value by default, although you can change this manually too:

```julia
param1.domain_size = 100
```

To simulate the script, you then simply run:

```julia
sol = simulate(model,param1)
```

*Note: as with all Julia functions, there is a significant compilation time associated with the first time you run this function. This overhead time will not be present if you re-run the simulation e.g., for different parameters*

The resulting PDE dynamics are represented by the solution object `sol` from [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/). To visualize the results, you can use a variety of plotting packages (e.g., Makie.jl, Plots.jl). We provide simple helper functions for the Plots.jl package.

To visualize the final pattern (once steady state has been reached), use: `endpoint()` 

```julia
using Plots
plot(endpoint(),model,sol)
```

To visualize intermediate timepoints (e.g., for making a movie of the dynamics), use `timepoint()`, e.g.,:

```julia
@gif for t in 0.0:0.01:1
    plot(timepoint(),model,sol,t)
end fps=20
```

Here, for example, `plot(timepoint(),model,sol,0.1)` plots the solution at a time that is 10% through the simulation. 