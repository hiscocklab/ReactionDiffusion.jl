module ReactionDiffusionTest

using ReactionDiffusion
using Test
include("../examples/Schnakenberg.jl")
using .Schnakenberg: U


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
    sol = simulate(model, params)

    u = sol[U][end]

    # Test whether simulated PDE is 'sensible'; we evaluate the max/min value of the final pattern, and also the number of sign changes about the half maximum (both for U)
    #       note:   we give a range for both test values as we are using random initial conditions, and thus variations are to be expected
    #               (even when setting seeds, it's not clear that Pkg updates to random will conserve values).
    @test 1.5 < dynamic_range(u) < 4
    @test num_periods(u) ≈ expected_periods rtol=0.1
end

# @testset "MOL" begin
#     model = Schnakenberg.model
#     params = dict(a = 0.2, b = 2.0, γ = 1.0, Dᵤ = 1.0, Dᵥ = 50.0, L=100.0)
#     expected_periods = 1/turing_wavelength(model,params)
#     sol = simulate(model, params; seed = 10)
#     u_mol,t_mol = simulate(model, params; discretisation=:mol, seed=10, reltol=1e-4)

#     u = sol[U][end]
#     u_mol = u_mol[:,1]
#     @test maximum(abs.(u_mol-u)) < 0.15
# end

@testset "initial conditions" begin
    initial = @initial_conditions begin
        U0 + exp(x), U
        V0, V
    end
    model = Model(Schnakenberg.reaction, Schnakenberg.diffusion, initial)
    params = dict(a = 0.2, b = 2.0, γ = 1.0, Dᵤ = 1.0, Dᵥ = 50.0, L=100.0, U0=1.0,V0=1.0, r=0.1)
    sol = simulate(model, params)
    u=sol[U]
    @show length(u[1])
    x = range(0,1,length(u[1]))
    @test u[1] ≈ 1.0.+exp.(x) rtol=1e-4
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
    n=512
    params = dict(a = 0.2, b = 2.0, γ = 1.0, Dᵤ = 1.0, Dᵥ = 50.0, L=L, r=0.1, g0=0.1, g1=0.2)
    sol = simulate(model, params; tspan=5.0, num_verts=n)
    u=sol[U][end]
    h = L/n
    @test (u[2] - u[1])/h  ≈ -0.1 rtol=0.25
    @test (u[end] - u[end-1])/h ≈ -0.2  rtol=0.25
end

end;