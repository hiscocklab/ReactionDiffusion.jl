# Example demonstrating the Schnakenburg model (a well-known reaction-diffusion system with analytical results for its Turing stability region)

using ReactionDiffusion

reaction = @reaction_network begin
    γ*a + γ*U^2*V,  ∅ --> U
    γ,              U --> ∅
    γ*b,            ∅ --> V
    γ*U^2,          V --> ∅
end 

diffusion = @diffusion_system L begin
    Dᵤ, U
    Dᵥ, V
end

model = Model(reaction, diffusion)
params = product(a = range(0.0,3.0,10), b = range(0.0,3.0,10), γ = 1.0, Dᵤ = 1.0, Dᵥ = 50.0, L=100.0)

##
sol=simulate(model,params; tspan=1000.0,  dt=0.01, num_verts=64, max_attempts=1)



using WGLMakie
using Makie
WGLMakie.activate!(resize_to=:body) # Make plots fill the available space.

U= reaction.U
f = Figure()
ax = Axis(f[1, 1])
for sol = sol.u
    lines!(ax,sol.x, sol[U][end])
end
timeseries_plot(model,sol.u)



f=turing_wavelength(model)


##

a = range(0,0.5;step=0.001)
b= range(0,3;step=0.001)
using Base.Threads

Λ = Matrix{Float64}(undef,length(a),length(b))



@threads for i in eachindex(a)
    for j in eachindex(b)
        p = dict(a = a[i], b = b[j], γ = 1.0, Dᵤ = 1.0, Dᵥ = 50.0, L=100.0)
        Λ[i,j] = f(p)
    end
end
heatmap(a,b,Λ)

