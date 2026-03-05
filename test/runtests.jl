module ReactionDiffusionTest

using ReactionDiffusion
using Test
include("../examples/Schnakenberg.jl")


function num_periods(u)
    deviation = sign.(u.- 0.5*(maximum(u) .+ minimum(u)))
    length(findall(!iszero,diff(deviation))) / 2
end

dynamic_range(u) = maximum(u)/minimum(u) 

##
@testset "turing_wavelength" begin
    model = Schnakenberg.model
    params = product(a = range(0.0,0.6,4), b = range(0.0,3.0,4), γ = [1.0], Dᵤ = [1.0], Dᵥ = [50.0], L = [100.0])

    λ = turing_wavelength.(model, params);
    turing_params = params[λ.>0]

    a = get.(turing_params, :a, missing)
    b = get.(turing_params, :b, missing)
    # Test whether the computed Turing parameters match the ground-truth Turing instability region
    @test a == [0.0, 0.2, 0.0, 0.2]
    @test b == [1.0, 1.0, 2.0, 2.0]
end

@testset "simulate" begin
    model = Schnakenberg.model
    params = dict(a = 0.2, b = 2.0, γ = 1.0, Dᵤ = 1.0, Dᵥ = 50.0, L=100.0, U0=1.0,V0=1.0)
    expected_periods = 1/turing_wavelength(model,params)
    u,t = simulate(model, params)

    u = u[:,1]

    # Test whether simulated PDE is 'sensible'; we evaluate the max/min value of the final pattern, and also the number of sign changes about the half maximum (both for U)
    #       note:   we give a range for both test values as we are using random initial conditions, and thus variations are to be expected
    #               (even when setting seeds, it's not clear that Pkg updates to random will conserve values).
    @test 1.5 < dynamic_range(u) < 4
    @test num_periods(u) ≈ expected_periods rtol=0.1
end

@testset "MOL" begin
    model = Schnakenberg.model
    params = dict(a = 0.2, b = 2.0, γ = 1.0, Dᵤ = 1.0, Dᵥ = 50.0, L=100.0)
    expected_periods = 1/turing_wavelength(model,params)
    u,t = simulate(model, params; seed = 10)
    u_mol,t_mol = simulate(model, params; discretisation=:mol, seed=10, reltol=1e-4)

    u = u[:,1]
    u_mol = u_mol[:,1]
    @test maximum(abs.(u_mol-u)) < 0.15
end

@testset "initial conditions" begin
    initial = @initial_conditions begin
        U0 + exp(x), U
        V0, V
    end
    model = Model(Schnakenberg.reaction, Schnakenberg.diffusion, initial)
    params = dict(a = 0.2, b = 2.0, γ = 1.0, Dᵤ = 1.0, Dᵥ = 50.0, L=100.0, U0=1.0,V0=1.0, r=0.1)
    u,t = simulate(model, params; full_solution=true)
    x = range(0,1,size(u,1))
    @test u[:,1,1] ≈ 1.0.+exp.(x) rtol=1e-4
end

@testset "boundary conditions" begin
    b0 = @reaction_network begin g0, ∅ --> U end
    b1 = @reaction_network begin g1, ∅ --> U end 
    initial = @initial_conditions begin
        x^2 * (g1-g0)/2 + x * g0, U
        0, V
    end
    model = Model(Schnakenberg.reaction, Schnakenberg.diffusion, (b0,b1), initial)
    L = 100.0
    n= 4096
    params = dict(a = 0.2, b = 2.0, γ = 1.0, Dᵤ = 1.0, Dᵥ = 50.0, L=L, r=0.1, g0=0.1, g1=0.2)
    u,t = simulate(model, params; tspan=5.0, num_verts=n)
    h = L/n
    @test (u[2,1] - u[1,1])/h  ≈ -0.1 rtol=0.1
    @test (u[end,1] - u[end-1,1])/h ≈ -0.2  rtol=0.1
end

end

module PseudoSpectralTest
include("../src/Util.jl")
include("../src/PseudoSpectral.jl")
using .PseudoSpectral
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