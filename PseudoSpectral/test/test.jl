using PseudoSpectral
using SciMLBase: solve
using Symbolics: @variables
using OrdinaryDiffEqExponentialRK: ETDRK4
using Test


@variables U,g0,g1,d,a,b
R = [0]
D = [1/(pi)^2] # Divide by pi^2 for a domain of size pi.
n=128
dt=0.001

@variables D
d=(1:3)/pi^2 |> collect  # Divide by pi^2 for a domain of size pi.
B = [0,0]
IC = [cos(pi*x)]
prob = PseudoSpectralProblem([U], R, [D], B, IC, n)
prob_func(prob,ctx) = remake(prob; p=Dict(D=>d[ctx.sim_id]))
output_func(sol,ctx) = (sol[U],false)
ensembleprob = EnsembleProblem(prob; prob_func=prob_func)
sol = solve(ensembleprob, ETDRK4(); tspan=(0.0,2.0), dt=dt, trajectories=length(d))
sol1=sol.u[1]