# ReactionDiffusion.jl for modelling pattern formation in biological systems

Reaction-diffusion dynamics are present across many areas of the physical and natural world, and allow complex spatiotemporal patterns to self-organize *de novo*. `ReactionDiffusion.jl` aims to be an easy-to-use and performant pipeline to simulate reaction-diffusion PDEs of arbitrary complexity, with a focus on pattern formation in biological systems. Using this package, complex, biologically-inspired reaction-diffusion models can be:

- specified using an intuitive syntax
- screened across millions of parameter combinations to identify pattern-forming networks (i.e., those that undergo a Turing instability)
- rapidly simulated to predict spatiotemporal patterns

## Quick start guide

Here we show how `ReactionDiffusion.jl` can be used to simulate a biologically-inspired reaction-diffusion system, responsible for generating evenly spaced joints along the length of your fingers and toes (from [Grall et el, 2024](https://www.pnas.org/doi/10.1073/pnas.2304470121)).

We begin by specifying the reaction-diffusion dynamics via the intuitive syntax developed in [Catalyst.jl](https://github.com/SciML/Catalyst.jl), which naturally mirrors biochemical feedbacks and interactions.

```@example quickstart
using ReactionDiffusion

model = @reaction_network begin
    # complex formation
    (k₊, k₋),               GDF5 + NOG <--> COMPLEX 
    # degradation
    δ₁,                     GDF5 --> ∅
    δ₂,                     NOG --> ∅
    δ₃,                     pSMAD --> ∅
    # transcriptional feedbacks (here: repressive hill functions)
    hillr(pSMAD,μ₁,K₁,n₁),  ∅ --> GDF5
    hillr(pSMAD,μ₂,K₂,n₂),  ∅ --> NOG
    # signalling
    μ₃*GDF5,                ∅ --> pSMAD
end  
```

We can then specify values for each parameter:

```@example quickstart
params = model_parameters()

# constant parameters
params.reaction["δ₃"] = [1.0]
params.reaction["μ₁"] = [1.0]
params.reaction["μ₃"] = [1.0]
params.reaction["n₁"] = [8]
params.reaction["n₂"] = [2]

# varying parameters
num_params = 5
params.reaction["δ₁"] = screen_values(min = 0.1,max = 10, number=num_params)
params.reaction["δ₂"] = screen_values(min = 0.1,max = 10,number=num_params)
params.reaction["μ₂"] = screen_values(min = 0.1,max = 10, number=num_params)
params.reaction["k₊"] = screen_values(min = 10.0,max = 100.0, number=num_params)
params.reaction["k₋"] = screen_values(min = 10.0,max = 100.0, number=num_params)
params.reaction["K₁"] = screen_values(min = 0.01,max = 1,number=num_params)
params.reaction["K₂"] = screen_values(min = 0.01,max = 1, number=num_params)
params.reaction

```

We must then also specify the diffusion coefficients for each variable:

```@example quickstart
params.diffusion = Dict(
    "NOG"       => [1.0],
    "GDF5"      => screen_values(min = 0.1,max = 30, number=10),
    "COMPLEX"   => screen_values(min = 0.1,max = 30, number=10)
)
params.diffusion
```

Then, with a single line of code, we can perform a Turing instability analysis across all combinations of parameters:

```@repl quickstart
turing_params = returnTuringParams(model, params);
```

**NEED TO UPDATE VERSION, TOTAL NUMBER OF PARAMETERS IS INCORRECT**

This returns all parameter combinations that can break symmetry from a homogeneous initial condition. We take advantage of the highly performant numerical solvers in [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) to be able to simulate millions of parameter sets per minute on a standard laptop. 

We may then take a single parameter set and simulate its spatiotemporal dynamics directly, using `Plots.jl` to visualize the resulting pattern:

```@example quickstart
param1 = get_params(model, turing_params[1000])
sol = simulate(model,param1)

using Plots
plot(endpoint(),model,sol)
```

We may also view the full spatiotemporal dynamics:

```@example quickstart
@gif for t in 0.0:0.01:1
    plot(timepoint(),model,sol,t)
end fps=20
```


## Support, citation and future developments

If you find `ReactionDiffusion.jl` helpful in your research, teaching, or other activities, please star the repository and consider citing this paper:

```
@article{TBD,
 doi = {TBD},
 author = {Muzatko, Daniel AND Daga, Bijoy AND Hiscock, Tom W.},
 journal = {biorXiv},
 publisher = {TBD},
 title = {TBD},
 year = {TBD},
 month = {TBD},
 volume = {TBD},
 url = {TBD},
 pages = {TBD},
 number = {TBD},
}
```

We are a small team of academic researchers who build mathematical models of developing embryos and tissues. We hope to extend the functionality of `ReactionDiffusion.jl` as funding and time allows. We also welcome suggestions and contributions from the community. 

