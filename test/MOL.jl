
using ReactionDiffusion
include("../examples/Schnakenberg.jl")    
model = Schnakenberg.model
params = dict(a = 0.2, b = 2.0, γ = 1.0, Dᵤ = 1.0, Dᵥ = 50.0, L=100.0, U0=1.0,V0=1.0)
expected_periods = 1/turing_wavelength(model,params)
u,t = simulate(model, params)
u_mol,t_mol = simulate(model, params; discretisation=:mol)


u = u[:,1]
u_mol = u_mol[:,1]
approx(u, u_mol)