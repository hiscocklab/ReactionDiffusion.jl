
using ReactionDiffusion
using Makie, WGLMakie
using BenchmarkTools
include("../examples/Schnakenberg.jl")    
model = Schnakenberg.model
params = dict(a = 0.2, b = 2.0, γ = 1.0, Dᵤ = 1.0, Dᵥ = 50.0, L=100.0, U0=1.0,V0=1.0)
expected_periods = 1/turing_wavelength(model,params)

n=64
sim_ps= simulate(model; seed=100, dt=0.1, max_attempts=1, num_verts=n)
sim_mol =  simulate(model; discretisation=:mol, seed=100,abstol=1e0, reltol=1e1, num_verts=n)


##
u_ps, t_ps = sim_ps(params)
u_mol, t_mol = sim_mol(params)
u_ps = u_ps[:,1]
u_mol = u_mol[:,1]
isapprox(u_ps, u_mol)
fig,ax = lines(u_ps)
lines!(ax, u_mol)
fig

##
@benchmark sim_ps(params)
@benchmark sim_mol(params)
