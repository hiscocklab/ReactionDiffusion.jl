module PseudoSpectralTest

using PseudoSpectral
using Symbolics: @variables
using OrdinaryDiffEqExponentialRK: ETDRK4
using SciMLBase: solve, successful_retcode
using Test

@testset "heat equation" begin
    @variables U,g0,g1,d,a,b
    R = [0]
    D = [1/(pi)^2] # Divide by pi^2 for a domain of size pi.
    n=128
    dt=0.001
    X = range(0,pi,n)
    @testset "zero flux" begin
        B = [0,0]
        IC = [cos(pi*x)]
        make_prob,transform = pseudospectral_problem([U], R, D, B, IC, n)
        prob = make_prob(Dict())
        sol = solve(prob, ETDRK4(); tspan=(0.0,2.0), dt=dt)
        @test successful_retcode(sol)
        u,t = transform(sol)
        @test u ≈ exp(-t)*cos.(X) rtol=1e-2;
    end
    @testset "non-zero flux" begin
        B = [pi,pi]
        IC = [pi*x]
        make_prob,transform = pseudospectral_problem([U], R, D, B, IC, n)
        prob = make_prob(Dict())
        sol = solve(prob, ETDRK4(); tspan=(0.0,2.0), dt=dt)
        @test successful_retcode(sol)
        u,t = transform(sol)
        @test u ≈ X rtol=1e-2;
    end
    @testset "non-negative" begin
        B = [-pi,0] # Inward flux
        IC = [1.0]
        make_prob,transform = pseudospectral_problem([U], R, D, B, IC, n)
        prob = make_prob(Dict())
        sol = solve(prob, ETDRK4(); tspan=(0.0,2.0), dt=dt)
        @test successful_retcode(sol)
        u,t = transform(sol)
        @test all(>(0), u)
    end
end

end;