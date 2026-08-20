# Determine pattern-forming parameter sets

Having specified the reaction-diffusion model, we can now determine which combinations of parameters demonstate Turing instability and compute the dominant wavelength of the resulting pattern.

The function `turing_wavelength` will calculate the dominant wavelength if Turing instability is predicted, or return 0.0 otherwise.

```julia
λ = turing_wavelength(model, params);
```

For convenience we also provide a multithreaded filter function which will return all the elements of `params` which exhibit Turing instability.

```julia
turing_params = filter_turing(model,params)
```