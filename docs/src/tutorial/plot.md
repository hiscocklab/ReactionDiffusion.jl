# Plot results
To visualise the results, you can use a variety of plotting packages (e.g., Makie.jl, Plots.jl). We provide a collection of interactive plotting tools built using `Makie.jl`.

To display the time series for a given parameter set,
```julia
using WGLMakie
timeseries_plot(model,params)
```
Use the slider to view the solution at different points in time.

For a fully interactive simulation, use the `interactive_plot` tool. This will allow you to adjust parameters and view the resulting solution in real time.

```julia
interactive_plot(model, params)
```

You can optionally pass ranges of possible values for some or all of the paramters.
```julia
interactive_plot(model, params; param_ranges=dict(μ₁ = range(0.1,10.0,100), μ₂ = range(0.1,10.0,5)))
```