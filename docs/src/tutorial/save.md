# Save results for future analysis

Sometimes you may wish to save some results for later analysis (for example, if you perform a particularly large parameter screen). You can use `JLD2` to save these files, e.g.,

```julia
using JLD2
##Save solution for a later date
jldsave("test.jld2"; model, params, u, t) 
```

Then, to load the files:

```julia
data = load("test.jld2")
model = data["model"]
```

