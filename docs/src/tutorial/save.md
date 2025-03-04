# Save results for future analysis

Sometimes you may wish to save some results for later analysis (for example, if you perform a particularly large parameter screen). You can use `JLD2` to save these files, e.g.,

```julia
##Save solution for a later date
@save "test.jld2" model params turing_params sol
```

Then, to load the files:

```julia
@load "test.jld2"
```

