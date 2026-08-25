## Import the ReactionDiffusion library.
using ReactionDiffusion, WGLMakie
WGLMakie.activate!(resize_to=:body) # Make plots fill the available space.

anterior = <(1/12)
posterior = >(11/12)

## Define a reaction network.
reaction = @reaction_network begin
    δ_bcd,                      BCD --> ∅   # Bicod degredation
    δ_nos,                      NOS --> ∅ # Nanos degredation
    hillar(BCD, NOS, μ_hb,K,n), ∅ --> HB  # Hunchback expression. Hill function with BCD activating and NOS inhibiting.
    δ_hb,                       HB --> ∅ # Hunchback degredation 
end

b0 = @reaction_network begin μ_bcd, ∅ --> BCD end # Bicoid expression at anterior.
b1 = @reaction_network begin μ_nos, ∅ --> NOS end  # Nanos expression at posterior.
boundary = (b0,b1)

## Define a system of diffusing species on a 1D domain of size `L`.
diffusion = @diffusion_system L begin
    D_bcd,  BCD
    D_nos,  NOS
    D_hb,   HB
end

## Combine `reaction` and `diffusion` into a single `Model` object.
model = Model(reaction, diffusion, boundary; name="Maternal Gradients")

## Pick some parameter sets to test.
params = dict(
    L = 2.0,
    μ_bcd = 1.0,
    μ_nos = 1.0,
    δ_bcd = 1.0,
    δ_nos = 1.0,
    μ_hb = 1.0,
    K = 1.0,
    n = 1.0,
    δ_hb = 1.0,
    D_bcd = 1.0,
    D_nos = 1.0,
    D_hb = 0.1
)

## Simulate the system with one of the "good" parameter sets and plot the results over time. 
# timeseries_plot(model, params; abstol=1e-3)

## Define some plausible ranges of parameter values to explore.
param_ranges = dict(
    L = range(1.0, 50.0, 50),
    μ_bcd = range(0.1, 2.0, 50),
    μ_nos = range(0.1, 2.0, 50),
    δ_bcd = range(0.1, 2.0, 50),
    δ_nos = range(0.1, 2.0, 50),
    μ_hb = range(0.1, 2.0, 50),
    K = range(0.1, 2.0, 50),
    n = range(1.0, 8.0, 50),
    δ_hb = range(0.1, 2.0, 50),
    D_bcd = range(0.1, 2.0, 50),
    D_nos = range(0.1, 2.0, 50),
    D_hb = range(0.1, 2.0, 50)
)




fig=interactive_plot(model, param_ranges; dt=0.01, num_verts=32)