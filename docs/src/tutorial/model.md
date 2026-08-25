# Specify the reaction-diffusion model

In the following sections, we go step-by-step through the code used in the Quick Start Guide.

We first specify the *reaction* part of the reaction-diffusion system, using the intuitive syntax developed in [Catalyst.jl](https://github.com/SciML/Catalyst.jl). This allows models to be written in a very natural way which reflects the (bio)-chemical interactions involved in the system: 

```julia
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
```

Here, reaction rates are assumed to follow mass action kinetics according to the stoichiometries of the reactants. So, for example, the `(k₊, k₋), GDF5 + NOG <--> COMPLEX` term represents the binding/unbinding of `GDF5` and `NOG` to form a `COMPLEX`. The rate of the forward reaction (i.e., the rate of complex formation) is $k_+ [\mathrm{GDF5}] [\mathrm{NOG}]$, and the rate of the reverse reaction (i.e., the rate of complex dissociation) is $k_- [\mathrm{COMPLEX}]$.

Reaction rates can also be overriden by user-specified functions, for example to denote regulatory feedbacks. In this example, the effect of `pSMAD` on `GDF5` expression is captured by the repressive hill function `hillr`: $\frac{\mu_1}{ 1 + ([\mathrm{pSMAD}]/K_1)^{n_1}}$. 

Arbitrary user-defined functions may also be used directly, for example the following code would reproduce the repressive hill function interaction:

```julia
function myOwnHillrFunction(input,μ,K,n)
    return μ/((input/K)^n + 1)
end

model = @reaction_network begin
# ...
    myOwnHillrFunction(pSMAD,μ₁,K₁,n₁),  ∅ --> GDF5
# ...
end

Second we define the *diffusion* part of the system by providing a set of transport reactions, each of which associate a diffusion rate parameter with one of the reactants. We also supply a parameter `L` which determines the size of the domain in which diffusion occurs.

diffusion = @diffusion_system L begin
    D_GDF5,     GDF5
    D_NOG,      NOG
    D_COMPLEX,  COMPLEX
end
```

Note, any reactants that are not assigned a diffusion constant are assumed to be non-diffusing (i.e., `pSMAD` in this example).

Finally we combine the two compontents to produce a `Model` object which describes the complete system.

```
model = Model(reaction, diffusion; name="BMP Signalling Dynamics")
```