# Example: Nodal-Lefty signalling

Here we illustrate a slightly more complex example, motivated by the Nodal-Lefty system that has been described in zebrafish embryos (e.g., [1](https://www.science.org/doi/full/10.1126/science.1221920), [2](https://elifesciences.org/articles/66397), [3](https://elifesciences.org/articles/54894)). We consider possible interactions between Nodal (a diffusible ligand), Lefty (a diffusible inhibitor), Oep (the relevant Nodal co-receptor), and pSmad2 (the intracellular transducer of Nodal signalling). We now show how a set of assumed biochemical interactions can be directly modelled by ReactionDiffusion.jl.

```@example nodal
using ReactionDiffusion

model = @reaction_network begin

    # Lefty binds to Nodal targeting it for degradation
    k₁,                     Nodal + Lefty --> ∅ 

    # Lefty binds to Oep targeting it for degradation
    k₂,                     Oep + Lefty --> ∅

    # Degradation of free Nodal, Lefty, and pSmad2
    δ₁,                     Nodal --> ∅
    δ₂,                     Lefty --> ∅
    δ₃,                     pSmad2 --> ∅

    # pSmad2 (i.e., Nodal signalling) activates Nodal transcription 
    hill(pSmad2,μ₁,K,n),  ∅ --> Nodal

    # pSmad2 (i.e., Nodal signalling) activates Lefty transcription 
    hill(pSmad2,μ₂,K,n),  ∅ --> Lefty

    # pSmad2 (i.e., Nodal signalling) activates Oep transcription (plus a basal level of Oep transcription)
    hill(pSmad2,μ₃,K,n) + basal,  ∅ --> Oep

    # Signalling occurs when Nodal and Oep combine to activate the intracellular signal transduction machinery
    β,                      Oep + Nodal --> pSmad2
end 
```

We then specify parameter ranges:

```@example nodal
params = model_parameters()

num_params = 5

params.reaction["δ₁"] = [1]
params.reaction["δ₂"] = [1.0]
params.reaction["δ₃"] = [1.0]


params.reaction["μ₁"] = [1.0]
params.reaction["μ₂"] = screen_values(min = 0.1,max = 10.0, number=num_params, mode="log")
params.reaction["μ₃"] = screen_values(min = 0.1,max = 10.0, number=num_params, mode="log")
params.reaction["basal"] = screen_values(min = 0.1,max = 1.0, number=num_params, mode="log")


params.reaction["K"] = screen_values(min = 0.01,max = 100.0, number=num_params, mode="log")
params.reaction["n"] = [4]


params.reaction["β"] = screen_values(min = 10.0,max = 100.0, number=num_params, mode="log")

params.reaction["k₁"] = screen_values(min = 10.0,max = 100.0, number=num_params, mode="log")
params.reaction["k₂"] = screen_values(min = 10.0,max = 100.0, number=num_params, mode="log")


params.diffusion = Dict(
    "Nodal"       => [1.0],
    "Lefty"      => screen_values(min = 1.0,max = 20.0, number=num_params, mode="log")
)
```

And then screen for pattern-forming parameter sets:

```@repl nodal
turing_params = returnTuringParams(model, params)
```

We now plot an example pattern-forming parameter set for further investigation.

```@example nodal
param1 = get_params(model, turing_params[1])
sol = simulate(model,param1)

using Plots
@gif for t in 0.0:0.01:1
    plot(timepoint(),model,sol,t)
end fps=20
```




