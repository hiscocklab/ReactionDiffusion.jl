# Determine pattern-forming parameter sets

Having specified the reaction-diffusion model and its parameters, we can use a single function to screen through all parameter combinations and return those that undergo a Turing instability:

```julia
turing_params = returnTuringParams(model, params);
```

*Note: as with all Julia functions, there is a significant compilation time associated with the first time you run this function. This overhead time will not be present if you re-run the simulation e.g., for different parameters*

This script returns all pattern-forming parameters found, saving the results in the `turing_params` variable, which has the following attributes:
- `steady_state_values`: denotes the (homogeneous) steady state values of each of the variables, in the absence of diffusion
- `pattern_phase`: predicts which components are in/out of phase.
- `wavelength`: denotes the approximate predicted wavelength of the pattern, assuming it is the one that is maximally unstable
- `max_real_eigval`: records the maximally unstable real eigenvalue
- `non_oscillatory`: if `true`, this is a stationary Turing pattern; if `false`, an oscillatory pattern.

`turing_params` also records the input values used to generate each pattern-forming parameter set, namely:
- `reaction_params`: records the reaction parameters 
- `diffusion_constants`: records the diffusion constants
- `initial_conditions`: records the initial condition used to compute the steady states

Individual parameters can also be accessed via helper functions:

```julia
δ₁ = get_param(model, turing_params,"δ₁","reaction")
D_COMPLEX = get_param(model, turing_params,"COMPLEX","diffusion")
```



