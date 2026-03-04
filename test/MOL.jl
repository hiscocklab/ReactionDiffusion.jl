
using ReactionDiffusion
using Makie, WGLMakie
using BenchmarkTools
include("../examples/Schnakenberg.jl")    
model = Schnakenberg.model

n=64
sim_ps= simulate(model; seed=101, dt=0.1, max_attempts=1, num_verts=n)
sim_mol =  simulate(model; discretisation=:mol, seed=101,abstol=1e-5, reltol=1e-2, num_verts=n)

##

u_ps, t_ps = sim_ps(Schnakenberg.params)
u_mol, t_mol = sim_mol(Schnakenberg.params)
u_ps = u_ps[:,1]
u_mol = u_mol[:,1]
isapprox(u_ps, u_mol)
fig,ax = lines(u_ps)
lines!(ax, u_mol)
fig

##
params = product(a = range(0.1,2.0,2), b = range(0.1,5.0,2), γ = 1.0, Dᵤ = 1.0, Dᵥ = 50.0, L=100.0, U0=1.0,V0=1.0)
# params = Schnakenberg.params
b_ps = @benchmark sim_ps(params)
b_mol = @benchmark sim_mol(params)

##
model = Schnakenberg.model
params = dict(a = 0.2, b = 2.0, γ = 1.0, Dᵤ = 1.0, Dᵥ = 50.0, L=100.0)
u_ps,t_ps = simulate(model, params; seed = 10)
u_mol,t_mol = simulate(model, params; discretisation=:mol, seed=10,  reltol=1e-5)

u_ps = u_ps[:,1]
u_mol = u_mol[:,1]
isapprox(u_ps, u_mol)
fig,ax = lines(u_ps)
lines!(ax, u_mol)
fig
