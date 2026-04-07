
# Features

- Allow for transport processes other than classical diffusion.
- Allow models to be input as equations.

# Performance

- Try GPU acceleration
- Incorporate benchmark code

# I/O

- Support default parameters as defined with `@parameters a=2.0` etc.
- Use Catalyst syntax for ICs
- Access solutions by species rather than integer index.
- Don't rely on consistant order of species.
- Accept parameter values given as integers.
- Add output expressions combining species.
- Add per species noise level to @initial_conditions
- Move noise level to @initial_conditions.
- use @debug for output such as repeated attempts.
- Use @ivs for x.

## Plotting

- Use LaTeX to label plots.
- ~~Normalise plots by range instead of norm~~
- Combine time series and interactive plots.
- Use observables.
- Show results for different parameters on the same plot.
- Make parameter set vectors and parameter ranges for interactive plotting the same objects (somehow).
- stop plots resizing as limits change.

# Documentation

- Add x variable to tutorial.

# Misc

- Use `getname` instead of `nameof`.
