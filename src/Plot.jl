module Plot
export timeseries_plot, interactive_plot
using ..Simulate
using ..Models
using LinearAlgebra: norm

using Printf: @sprintf
using Makie
using Observables

"""
    timeseries_plot(model, params; normalise=true, hide_y=true, autolimits=true, kwargs...)

Simulate and plot the results. The remaining `kwargs` are passed to `simulate`.
"""
function timeseries_plot(model, params; normalise=true, hide_y=true, autolimits=true, species=nothing, kwargs...)
    u,t=simulate(model,params; full_solution=true, kwargs...)
    timeseries_plot(model, u,t; normalise=normalise, hide_y=hide_y, autolimits=autolimits, species=species)
end


"""
    function timeseries_plot(model, u, t; normalise=true, hide_y=true, autolimits=true, kwargs...)

Display a solution in an interactive plot with a scrubber to move through time.

If `normalise` is true, values for different species will be normalised to a common scale.
"""
function timeseries_plot(model, u, t; normalise=true, hide_y=true, autolimits=true, species=nothing, kwargs...)
    model_species = Models.species(model)
    if isnothing(species)
        labels = [string(s.f) for s in model_species]
    else
        ix = [i for (i,s) in enumerate(model_species) if nameof(s.f) ∈ species]
        u = u[:,ix,:]
        labels = [string(s.f) for s in model_species[ix]]
    end

    x_steps = size(u, 1)
    x = range(0.0,1.0,length=x_steps)
	r = normalise ? norm.(eachslice(u, dims=(2,3))) : ones(size(u)[2:3])
	fig=Figure()
	ax = Axis(fig[1,1])
	hide_y && hideydecorations!(ax)
    sg = SliderGrid(fig[2,1], (label="t",range=eachindex(t), format=i->@sprintf("%.2f",t[i])))
    sl=sg.sliders[1]
	T = lift(i -> t[i], sl.value)
	U = [lift(i -> u[:,i]/r[i], sl.value) for (u,r) in zip(eachslice(u, dims=2), eachrow(r))]
	for (U,label) in zip(U,labels)
		lines!(ax,x,U, label=label)
	end
	autolimits && on(sl.value) do _
	    autolimits!(ax)
	end
	axislegend(ax)
    display(fig)
    fig
end

"""
    interactive_plot(model, param_ranges; hide_y=true, num_verts=32, kwargs...)

Generate an interactive plot of the steady state solution with sliders to adjust each of the parameters within `param_ranges`.
`param_ranges` should be a dictionary mapping parameter names to either `Range` objects or collections of possible values.
"""
function interactive_plot(model, param_ranges; normalise=true, hide_y=true, num_verts=32, kwargs...)
    simulate_ = simulate(model; num_verts=num_verts, kwargs...)
    function f(vals...)
        params = Dict(k => x isa Int ? v[x] : x for ((k,v), x) in zip(param_ranges,vals))
        u,t = parameter_set(model,params) |> simulate_
        r = norm.(eachcol(u))
        normalise ? u ./ r' : u
    end
    
	fig=Figure()
	ax = Axis(fig[1,1])
	hide_y && hideydecorations!(ax)

    param_ranges = sort(param_ranges)
    slider_specs = [eltype(v) <: AbstractFloat ? (label=string(k), range = v, format = x->@sprintf("%.2f",x)) : (label=string(k), range = 1:length(v)) for (k,v) in param_ranges]

    sg = SliderGrid(fig[1,2], slider_specs...)

    U = lift(f, (sl.value for sl in sg.sliders)...)
    U = throttle(1/120, U) # Limit update rate to 120Hz
    x = range(0,1,num_verts)
    labels = [string(s.f) for s in species(model)]
    for i in eachindex(eachcol(U[]))
        lines!(ax, x, lift(u -> u[:,i], U); label=labels[i])
    end
    on(U) do _
	    autolimits!(ax)
	end
    axislegend(ax)
    display(fig)
    fig
end
end