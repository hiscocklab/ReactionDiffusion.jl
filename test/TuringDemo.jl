
using ReactionDiffusion
using WGLMakie, Makie

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
D = range(0.0,100.0,100)
params = [dict(a = 0.2, b = 2.0, γ = 1.0, Dᵤ = 1.0, Dᵥ = d, L=100.0) for d in D]

λ = map(turing_wavelength(model), params)

f = Figure()
ax = Axis(f[1, 1],
    title = "Schnakenberg System - Dominant Wavelength",
    xlabel = "Dᵥ",
    ylabel = "λ",
)
hidedecorations!(ax; label=false, ticks=false, ticklabels=false)
lines!(ax, D, λ)