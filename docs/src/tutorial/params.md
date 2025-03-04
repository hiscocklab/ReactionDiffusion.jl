# Specify the parameter values

We must now define parameters, starting with the reaction terms.

We may specify constant parameter values:

```julia
params = model_parameters()

# constant parameters
params.reaction["δ₃"] = [1.0]
params.reaction["μ₁"] = [1.0]
params.reaction["μ₃"] = [1.0]
params.reaction["n₁"] = [8]
params.reaction["n₂"] = [2]
```

or ranges of parameter values:

```julia
num_params = 5
params.reaction["δ₁"] = screen_values(min = 0.1,max = 10, number=num_params)
params.reaction["δ₂"] = screen_values(min = 0.1,max = 10,number=num_params)
params.reaction["μ₂"] = screen_values(min = 0.1,max = 10, number=num_params)
params.reaction["k₊"] = screen_values(min = 10.0,max = 100.0, number=num_params)
params.reaction["k₋"] = screen_values(min = 10.0,max = 100.0, number=num_params)
params.reaction["K₁"] = screen_values(min = 0.01,max = 1,number=num_params)
params.reaction["K₂"] = screen_values(min = 0.01,max = 1, number=num_params)
params.reaction
```

Here, the function `screen_values` returns a series of `num_params` parameters that are linearly spaced between the `min` and `max` limits. The argument `mode="log"` may be used to sample in log-space instead. 

Arbitrary collections of parameters may also be specified, e.g.,

```julia
params.reaction["n₁"] = [2; 4; 8]
```

Next, we must specify the diffusion coefficients:

```julia
params.diffusion = Dict(
    "NOG"       => [1.0],
    "GDF5"      => screen_values(min = 0.1,max = 30, number=10),
    "COMPLEX"   => screen_values(min = 0.1,max = 30, number=10)
)
```

Note, any variables that are not assigned a diffusion constant are assumed to be non-diffusing (i.e., `pSMAD` in this example).

Optionally, the initial conditions used to define the steady state(s) of the system may be explicitly set:

```julia
params.initial_condition["NOG"] = [0.1;1.0,10.0]
```

This may be particularly helpful when there are more than one steady states.
