using Revise
using ReactionDiffusion
using ReactionDiffusion.Util: isnonzero, unzip_dict, lookup
using Makie, WGLMakie
using Logging: with_logger
const n = 16

##
function plot_turing(model, params)

    λ = turing_wavelength(model, params)
    a = [p[:a] for p in params]
    b = [p[:b] for p in params]

    fig = Figure()

    ax = Axis(fig[1,1], title="Schnakenberg Turing Instabilities", xlabel="a", ylabel="b")
    heatmap!(ax, a, b, λ, colormap=:grays)
    fig
end
##



##

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

cost(u,t) = sum(u[n÷2:end,1]) / sum(u[:,1])

function sample(params)
    local n = 20
    σ = 15
    params = copy.(fill(params,n))
    for p in params
        p[:L] += randn() * sigma
    end
    params
end

##
params = product(a = range(0.0,0.5,10), b = range(0.0,5.0,10), γ = [1.0], Dᵤ = [1.0], Dᵥ = [50.0], L=[50.0])
turing_params = filter_turing(model,params)
vars = lookup.([:a,:b])
a,b,c = plot_cost(model, cost, vars, turing_params; maxrepeats=0, maxiters=1e5)
paths = @tmap(turing_params[1:5:end], Vector{Tuple{Float64}}) do params
    path = optimise(model, cost, vars, params; in_domain = is_turing(model), maxsteps=100, maxrepeats=0, maxiters=1e4, savepath=true)
    [(a,b) for (a,b) in path]
end
paths = [[(a,b) for (a,b) in path] for path in paths]

fig = Figure()
ax = Axis(fig[1,1], title="Schnakenberg Gradient Descent", xlabel="a", ylabel="b")
heatmap!(ax, a, b, c, colormap=:redblue, colorrange=(0.0,1.0))
for path in paths
    lines!(ax, path)
end
fig

