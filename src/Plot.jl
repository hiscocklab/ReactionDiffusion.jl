module Plot
export steady_state_plot, timeseries_plot, interactive_plot
using PseudoSpectralReactionDiffusion
using ..Simulate
using ..Models
using LinearAlgebra: norm

using Printf: @sprintf
using Makie
using Makie: COLOR_ACCENT
using CairoMakie
using Observables
using Random: Xoshiro, seed!
using NativeFileDialog: save_file


"""
    function steady_state_plot(model, params; normalise=true, hide_y=normalise, kwargs...)

Run a simulation and plot the final values for each species.
If `normalise` is true, values for different species will be normalised to a common scale.
"""
function steady_state_plot(model, params; normalise=false, hide_y=normalise, kwargs...)
    sol = simulate(model,params; kwargs...)
    steady_state_plot(model,sol; normalise, hide_y)
end

"""
    function steady_state_plot(model, sol; normalise=true, hide_y=normalise, kwargs...)

Plot the final value in the solution object for each species.
If `normalise` is true, values for different species will be normalised to a common scale.
"""
function steady_state_plot(model, sol::PseudoSpectralSolution; normalise=false, hide_y=normalise, species=nothing)
    u = eachcol(sol.u[end])
    model_species = Models.species(model)
    if isnothing(species)
        labels = [string(s.f) for s in model_species]
    else
        ix = [i for (i,s) in enumerate(model_species) if nameof(s.f) ∈ species]
        u = u[ix]
        labels = [string(s.f) for s in model_species[ix]]
    end

    num_verts = length(u[1])
    x = range(0.0, 1.0, num_verts)

    if normalise
        @. u = u / maximum(u)
    end

	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "x / L", ylabel = "Concentration") # TODO have axis labels passed in as kwargs
	hide_y && hideydecorations!(ax; label=false)
	for (u, label) in zip(u, labels)
		lines!(ax, x, u; label)
	end

	autolimits!(ax)

	axislegend(ax)
    display(fig)
    fig
end
"""
    timeseries_plot(model, params; normalise=true, hide_y=normalise, autolimits=true, kwargs...)

Simulate and display the results with an interactive slider to move through time. The remaining `kwargs` are passed to `simulate`.
If `normalise` is true, values for different species will be normalised to a common scale.
"""
function timeseries_plot(model, params; normalise=false, hide_y=normalise, autolimits=true, species=nothing, kwargs...)
    sol = simulate(model, params; full_solution = true, kwargs...)
    timeseries_plot(model, sol; normalise, hide_y, autolimits, species)
end

"""
    timeseries_plot(model, sol; normalise=true, hide_y=normalise, autolimits=true, kwargs...)

Display the solution with an interactive slider to move through time.
If `normalise` is true, values for different species will be normalised to a common scale.
"""
function timeseries_plot(model, sol::PseudoSpectralSolution; normalise=false, hide_y=normalise, autolimits=true, species=nothing, kwargs...)
    u = stack(sol.u)
    timeseries_plot(model, u, sol.t; normalise, hide_y, autolimits, species)
end


"""
    timeseries_plot(model, u, t; normalise=true, hide_y=normalise, autolimits=true, kwargs...)

Display the solution with an interactive slider to move through time.
If `normalise` is true, values for different species will be normalised to a common scale.
"""
function timeseries_plot(model, u, t; normalise=false, hide_y=normalise, autolimits=true, species=nothing, kwargs...)
    model_species = Models.species(model)
    if isnothing(species)
        labels = [string(s.f) for s in model_species]
    else
        ix = [i for (i,s) in enumerate(model_species) if nameof(s.f) ∈ species]
        u = u[:,ix,:]
        labels = [string(s.f) for s in model_species[ix]]
    end

    num_verts = size(u, 1)
    x = range(0.0, 1.0, num_verts)
	r = normalise ? maximum.(eachslice(u, dims = 2)) : ones(size(u)[2:3]) # find max value for each species across all time

	fig = Figure()
	ax = Axis(fig[1,1], xlabel = "x / L", ylabel = "Concentration") # TODO have axis labels passed in as kwargs
	hide_y && hideydecorations!(ax;label=false)
    sg = SliderGrid(fig[2,1], (label = "t", range = eachindex(t), format = i -> @sprintf("%.2f", t[i])))

    sl = sg.sliders[1]
	T = lift(i -> t[i], sl.value)
	U = [lift(i -> u[:, i] / r, sl.value) for (u, r) in zip(eachslice(u, dims = 2), r)]
	for (U, label) in zip(U, labels)
		lines!(ax, x, U, label = label)
	end

    if autolimits
        on(sl.value) do _
            autolimits!(ax)
        end
	end

	axislegend(ax)
    display(fig)
    fig
end

"""
    interactive_plot(model, param_ranges; normalise=false, hide_y=normalise, autolimits=true, num_verts=32, kwargs...)

Generate an interactive plot of the steady state solution with sliders to adjust each of the parameters within `param_ranges`.
`param_ranges` should be a dictionary mapping parameter names to either `Range` objects or collections of possible values.
"""
function interactive_plot(model, params; param_ranges=nothing, normalise=false, hide_y=normalise, autolimits=true, tspan=Inf64, tol=1e-5, num_verts=32, dt=0.01, seed=123, kwargs...)
    params = sort(params)
    param_ranges = something(param_ranges,
        Dict(k => range(0.0,2*v,100) for (k,v) in params)) # Parameter sliders go from 0 to 2*params by default.
    fig = Figure()
    layout = make_layout(fig)
    ax = Axis(layout.ax, title=name(model), xlabel = "x / L", ylabel = "Concentration", yticklabelspace = 50.0)
    hide_y && hideydecorations!(ax;label=false)
    param_sliders = make_param_sliders(layout.param_sliders, param_ranges, params; width=300)
    TMAX = Observable(0.0)
    state = Observable(:stop)
    recording = Observable(false)

    t_slider_grid = SliderGrid(layout.t_slider, (label = "t", range = @lift(0.0:0.01:$TMAX), format = "{:.2f}"))
    t_slider = t_slider_grid.sliders |> only
    play_button = Button(layout.play_button; label=@lift(if $state==:stop "▶" else "⏸" end), font="Segoe UI Symbol")
    reset_button = Button(layout.reset_button; label="⏮")
    skip_button = Button(layout.skip_button; label="⏭", buttoncolor=@lift(if $state==:ff COLOR_ACCENT[] else RGBf(0.94, 0.94, 0.94) end))
    record_button = Button(layout.record_button; label="⏺", buttoncolor=@lift(if $recording COLOR_ACCENT[] else RGBf(0.94, 0.94, 0.94) end), font="Segoe UI Symbol")
    capture_button = Button(layout.capture_button; label="📷", font="Segoe UI Symbol")

    

    

    function f(t)
        # step_to!(int, t)
        # sol = get_sol(int)
        # i = findfirst(>=(t), sol.t)
        u = get_u(int[],t)
        r = maximum.(eachcol(u))
        normalise ? u ./ r' : u
    end

    
    P = (throttle(1/120,sl.value) for sl in param_sliders.sliders)

    params = parameter_set(model, Dict(k => x isa Int ? v[x] : x for ((k, v), x) in zip(param_ranges, [p[] for p in P])))
    rng=Xoshiro(seed)
    prob = PseudoSpectralProblem(model, num_verts; p=params, dt, rng, tspan, kwargs...)
    callback = isinf(tspan) ? steady_state_callback(tol) : nothing # Only stop at steady state if tspan isn't specified.

    int = Ref(init(prob; callback))

    onany(P...) do p...
        params = parameter_set(model, Dict(k => x isa Int ? v[x] : x for ((k, v), x) in zip(param_ranges, p)))
        rng = seed!(rng,seed)
        int[] = remake(int[]; p=params, rng)
        U[]=f(T[])
    end

    RealT = Observable(0.0)
    TT = throttle(1/120, RealT)
    T = @lift min($TT, int[].ss, tspan)
    U = lift(f, T)

    x = range(0, 1, num_verts)
    labels = [string(s.f) for s in species(model)]
    for i in eachindex(eachcol(U[]))
        lines!(ax, x, lift(u -> u[:, i], U); label=labels[i])
    end

    if autolimits
        on(U) do _
            autolimits!(ax)
        end
    end
    legend = Legend(layout.legend, ax; orientation= :horizontal)
    



    on(play_button.clicks) do _
        if state[] == :stop
            state[] = :play
        else
            state[] = :stop
        end
    end


    on(events(fig).tick) do tick
        if state[] == :play
            RealT[] += tick.delta_time
        elseif state[] == :ff
            RealT[] += (RealT[] + 1.0)*tick.delta_time
        end
        if isassigned(video)
            recordframe!(video[])
        end
    end

    on(T) do t
        t_end = isinf(tspan) ? int[].ss : tspan
        if t >= t_end
            t = t_end
            state[]= :stop
        end
        TMAX[] = max(TMAX[],t)
        set_close_to!(t_slider, t)
    end

    on(t_slider.value) do t
        # Hack to distinguish user interaction (large movements) from ticks (small movements).
        # TODO tune this so it works more reliably.
        isapprox(t, RealT[]; atol=0.05) && return 
        state[] = :stop
        RealT[] = t
    end

    on(reset_button.clicks) do _
        state[] = :stop
        RealT[] = 0.0
        set_close_to!(t_slider, RealT[])
    end

    on(skip_button.clicks) do _
        ss= int[].ss
        if ss < Inf
            state[] = :stop
            RealT[] = ss
            set_close_to!(t_slider, RealT[])
        else
            state[] = :ff
        end 
    end

    on(capture_button.clicks) do _
        filename = save_file(pwd(); filterlist = "png;svg;pdf")
        isempty(filename) && return
        export_fig = Figure()
        export_ax = Axis(export_fig[1, 1]; title=ax.title[], xlabel=ax.xlabel[], ylabel=ax.ylabel[])
        for i in eachindex(eachcol(U[]))
            lines!(export_ax, x, U[][:, i]; label=labels[i])
        end
        save(filename, export_fig; backend=CairoMakie)
        nothing
    end

    video = Ref{VideoStream}()
 
    on(record_button.clicks) do _
        if recording[]
            filename = save_file(pwd(); filterlist = "mkv;mp4;webm;gif")
            if !isempty(filename)
                save(filename, video[])
            end
            video = Ref{VideoStream}()
            recording[]=false
            return
        end
        export_fig = Figure()
        export_ax = Axis(export_fig[1, 1]; title=ax.title[], xlabel=ax.xlabel[], ylabel=ax.ylabel[])
        for i in eachindex(eachcol(U[]))
            lines!(export_ax, x, lift(u -> u[:, i], U); label=labels[i])
        end

        video[] = VideoStream(export_fig; backend=CairoMakie)
        recording[]=true
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
        param_sliders = body[1,2]
    control_bar = fig[2,1]
        t_slider = control_bar[1,1]
        control_buttons = control_bar[1,2]
            play_button = control_buttons[1,1]
            reset_button = control_buttons[1,2]
            skip_button = control_buttons[1,3]
            record_button = control_buttons[1,4]
            capture_button = control_buttons[1,5]
    (;ax, legend, param_sliders, t_slider,
        play_button, reset_button, skip_button, record_button, capture_button)
end

function make_param_sliders(f, param_ranges, params; width=nothing)
    slider_specs = [eltype(v) <: AbstractFloat ? (label=string(k), range = v, startvalue=params[k], format = x -> @sprintf("%.2f",x)) : (label=string(k), range = 1:length(v)) for (k,v) in param_ranges]
    SliderGrid(f, slider_specs...; width=width)
end

end
