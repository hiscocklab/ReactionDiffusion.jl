```@raw html
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/banner-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="assets/banner.svg">
  <img alt="ReactionDiffusion.jl" src="assets/banner.svg">
</picture>
```

# A simulation and visualisation package for 1D reaction-diffusion systems.

Reaction-diffusion dynamics are present across many areas of the physical and natural world, and allow complex spatiotemporal patterns to self-organize *de novo*. `ReactionDiffusion.jl` aims to be an easy-to-use and performant pipeline to simulate reaction-diffusion PDEs of arbitrary complexity, with a focus on pattern formation in biological systems. Using this package, complex, biologically-inspired reaction-diffusion models can be:

- specified using an intuitive syntax
- screened across millions of parameter combinations to identify pattern-forming networks (i.e., those that undergo a Turing instability)
- rapidly simulated to predict spatiotemporal patterns

## Quick start guide

Here we show how `ReactionDiffusion.jl` can be used to simulate a biologically-inspired reaction-diffusion system, responsible for generating evenly spaced joints along the length of your fingers and toes (from [Grall et el, 2024](https://www.pnas.org/doi/10.1073/pnas.2304470121)).

We begin by specifying the reaction-diffusion dynamics via the intuitive syntax developed in [Catalyst.jl](https://github.com/SciML/Catalyst.jl), which naturally mirrors biochemical feedbacks and interactions.

```julia
using ReactionDiffusion

reaction = @reaction_network begin
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

diffusion = @diffusion_system L begin
    D₁, GDF5,
    D₂, NOG,
    D₃, COMPLEX
end

model = Model(reaction, diffusion)
```

We can then specify possible values for each parameter:

```julia
params = product(
    μ₁ = [1.0],
    μ₂ = range(0.1,10,5),
    k₊ = range(10,100, 5),
    k₋ = range(10,100,5),
    μ₃ = [1.0],
    δ₁ = range(0.1,10,5),
    δ₂ = range(0.1,10,5),
    δ₃ = [1.0],
    K₁ = range(0.01,1.0,5),
    K₂ = range(0.01,1.0,5),
    n₁ = [8.0],
    n₂ = [2.0],
    D₁ = [1.0],
    D₂ = range(0.1,30,10),
    D₃ = range(0.1,30,10)
)
```

Then, with a single line of code, we can perform a Turing instability analysis across all combinations of parameters:

```julia
turing_params = filter_turing(model, params);
```

This returns all parameter combinations that can break symmetry from a homogeneous initial condition.

We may then take a collection of parameter sets and simulate their spatiotemporal dynamics directly.

```julia
sol = simulate(model,turing_params)
```

We may also view the full spatiotemporal dynamics of a particular parameter set:

```julia
using WGLMakie
timeseries_plot(model,turing_params[4])
```

Or view the results of adjusting parameters of interest in real time.

```julia
param_ranges = dict(
    μ₁ = range(0.5,2.0,100),
    μ₂ = [x->exp(-λ*x) for λ in range(1.0,10.0,100)],
    D₃ = range(1.0,100.0,100)
)
interactive_plot(model, param_ranges)
```

We take advantage of psuedo-spectral methods in combination with the highly performant numerical solvers in [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) to be able to simulate millions of parameter sets per minute on a standard laptop. 


## Support, citation and future developments

If you find `ReactionDiffusion.jl` helpful in your research, teaching, or other activities, please star the repository and consider citing [this paper](https://doi.org/10.1242/dev.205067).

We are a small team of academic researchers from the [Hiscock Lab](https://twhiscock.github.io/), who build mathematical models of developing embryos and tissues. We have found these scripts helpful in our own research, and make them available in case you find them helpful in your research too. We hope to extend the functionality of `ReactionDiffusion.jl` as our future projects, funding and time allows.

This work is supported by ERC grant SELFORG-101161207, and UK Research and Innovation (Biotechnology and Biological Sciences Research Council, grant number BB/W003619/1) 

*Funded by the European Union. Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Council Executive Agency. Neither the European Union nor the granting authority can be held responsible for them*

![ERC_logo](./assets/LOGO_ERC-FLAG_FP.png)

