
using PseudoSpectral
using Symbolics: @variables
using OrdinaryDiffEqExponentialRK: ETDRK4
using SciMLBase: successful_retcode
using Test
##
@variables U,g0,g1,d,a,b
R = [0]
D = [1/(pi)^2] # Divide by pi^2 for a domain of size pi.
n=128
dt=0.001
X = range(0,pi,n)
B = [0,0]
IC = [cos(pi*x)]
##
prob = PseudoSpectralProblem([U], R, D, B, IC, n; p=Dict())
sol = solve(prob, ETDRK4(); tspan=(0.0,2.0), dt=dt)
@show successful_retcode(sol.sol)
@show isapprox(u, exp(-t)*cos.(X); rtol=1e-2)