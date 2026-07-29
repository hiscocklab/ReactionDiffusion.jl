using ReactionDiffusion
using PseudoSpectral
using SciMLBase: init
using OrdinaryDiffEqExponentialRK: ETDRK4
include("../examples/Schnakenberg.jl")
using .Schnakenberg: U,model

params = dict(a = 0.2, b = 2.0, γ = 1.0, Dᵤ = 1.0, Dᵥ = 50.0, L=100.0, U0=1.0,V0=1.0)
p = parameter_set(model,params)
sol = simulate(model, params)
prob = PseudoSpectralProblem(model, 64; noise=0.0, dt=0.1, p)
@show prob.ode_problem

A = prob.ode_problem.f.f1

init(prob.ode_problem, ETDRK4(), dt=0.1, tspan=(0.0,0.1))