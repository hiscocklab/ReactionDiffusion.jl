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

Each parameter in `param_ranges` can be adjusted within the given range using the sliders. Here we provide a range of functions for μ₂ which produce varying spatial gradients.
