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
    u, t = simulate(model, params; full_solution = true, kwargs...)
    timeseries_plot(model, u, t; normalise = normalise, hide_y = hide_y, autolimits = autolimits, species = species)
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
    x = range(0.0, 1.0, length = x_steps)
	r = normalise ? maximum.(eachslice(u, dims = 2)) : ones(size(u)[2:3]) # find max value for each species across all time

	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "x / L", ylabel = "Concentration") # TODO have axis labels passed in as kwargs
	hide_y && hideydecorations!(ax)
    sg = SliderGrid(fig[2,1], (label = "t", range = eachindex(t), format = i -> @sprintf("%.2f", t[i])))

    sl = sg.sliders[1]
	T = lift(i -> t[i], sl.value)
	U = [lift(i -> u[:, i] / r, sl.value) for (u, r) in zip(eachslice(u, dims = 2), r)]
	for (U, label) in zip(U, labels)
		lines!(ax, x, U, label = label)
	end

	!normalise && autolimits && on(sl.value) do _
	    autolimits!(ax)
	end

    dy = 0.05 # padding of the plot in y
    normalise ? ylims!(ax, -dy, 1 + dy) : nothing

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
	param_ranges = sort(param_ranges)
    fig = Figure()
    layout = make_layout(fig)
    ax = Axis(layout.ax, title=name(model), xlabel = "x / L", ylabel = "Concentration")
    hide_y && hideydecorations!(ax)
    param_sliders = make_param_sliders(layout.param_sliders, param_ranges)
    save_button = Button(layout.save_button, label="🖫", font="Segoe UI Symbol")
    annotate_button = Button(layout.annotate_button, label="👁", font="Segoe UI Symbol")
    TMAX = Observable(0.0)
    t_slider_grid = SliderGrid(layout.t_slider, (label = "t", range = @lift(0.0:0.01:$TMAX), format = "{:.2f}"))
    t_slider = t_slider_grid.sliders |> only
    play_button = Button(layout.play_button; label="⏯")
    reset_button = Button(layout.reset_button; label="⏮")
    skip_button = Button(layout.skip_button; label="⏭")
    record_button = Button(layout.record_button; label="⏺", font="Segoe UI Symbol")
    capture_button = Button(layout.capture_button; label="📷", font="Segoe UI Symbol")



    get_sol! = make_integrator(model; num_verts=num_verts, kwargs...)

    function f(T, P...)
        params = parameter_set(model, Dict(k => x isa Int ? v[x] : x for ((k, v), x) in zip(param_ranges, P)))
        u = get_sol!(params, T[])
        r = maximum.(eachcol(u))
        normalise ? u ./ r' : u
    end

    
    P = (sl.value for sl in param_sliders.sliders)
    RealT = Observable(0.0)
    T = throttle(1/120, RealT)

    U = lift(f, T, P...)
    U = throttle(1/120, U) # Limit update rate to 120Hz

    x = range(0, 1, num_verts)
    labels = [string(s.f) for s in species(model)]
    for i in eachindex(eachcol(U[]))
        lines!(ax, x, lift(u -> u[:, i], U); label=labels[i])
    end

    !normalise && on(U) do _
	    autolimits!(ax)
	end
    dy = 0.05 # padding of the plot in y
    normalise ? ylims!(ax, -dy, 1 + dy) : nothing

    legend = Legend(layout.legend, ax; orientation= :horizontal)
    

    play = Observable(false)

    on(play_button.clicks) do _
        play[] = !play[]
    end

    on(events(fig).tick) do tick
        if play[]
            t = RealT[] + tick.delta_time
            RealT[] = t
            TMAX[] = max(TMAX[],T[])
            set_close_to!(t_slider, RealT[])
        end
    end

    on(t_slider.value) do t
        isapprox(t, RealT[]; atol=0.02) && return
        play[] = false
        RealT[] = t
    end

    on(reset_button.clicks) do _
        play[] = false
        RealT[] = 0
        set_close_to!(t_slider, RealT[])
    end
 

    display(fig)
    fig
end

function make_layout(fig)
    body = fig[1,1]
        plot_pane = body[1,1]
            ax = plot_pane[1,1]
            legend_bar = plot_pane[2,1]
                legend = legend_bar[1,1]
                param_buttons = legend_bar[1,2]
                    save_button = param_buttons[1,1]
                    annotate_button = param_buttons[1,2]
        param_sliders = body[1,2]
    control_bar = fig[2,1]
        t_slider = control_bar[1,1]
        control_buttons = control_bar[1,2]
            play_button = control_buttons[1,1]
            reset_button = control_buttons[1,2]
            skip_button = control_buttons[1,3]
            record_button = control_buttons[1,4]
            capture_button = control_buttons[1,5]
    (;ax, legend, save_button, annotate_button, param_sliders, t_slider,
        play_button, reset_button, skip_button, record_button, capture_button)
end

function make_param_sliders(f, param_ranges)
    slider_specs = [eltype(v) <: AbstractFloat ? (label=string(k), range = v, format = x -> @sprintf("%.2f",x)) : (label=string(k), range = 1:length(v)) for (k,v) in param_ranges]
    SliderGrid(f, slider_specs...)
end

end
