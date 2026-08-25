# Installation and setup

ReactionDiffusion can be installed using the Julia package manager:

```julia
using Pkg
Pkg.add("ReactionDiffusion")
```

and then loaded ready for use:

```julia
using ReactionDiffusion
```

As with all Julia packages, upon first loading there will significant precompilation time, but this step is only required once per installation.

You may also wish to install a plotting frontend to visualise your results; we provide plotting functions which use [Makie.jl](https://makie.org).

```julia
Pkg.add("WGLMakie")
using WGLMakie
WGLMakie.activate!(resize_to=:body) # Make plots fill the available space.
```

If you don't already use Julia, this must first be installed, along with an environment to edit and run your scripts. We recommend following the guides from [SciML](https://docs.sciml.ai/Overview/stable/getting_started/installation/) to set this up. 

We also recommend using environments to manage your package dependencies; see [here](https://docs.sciml.ai/Catalyst/stable/introduction_to_catalyst/catalyst_for_new_julia_users/#catalyst_for_new_julia_users_packages) for a more detailed discussion. 

To take advantage of automatic multithreading in your simulations, ensure that the number of Julia threads is set according to the CPU on your machine (we recommend to set `Threads.nthreads()` to be equal to around double the number of cores available). In VSCode, this can be achieved by the following:
- go to Settings `Ctrl-,`
- search for "threads"
- click `edit in settings.json`
- set `"julia.NumThreads": 24` (for a machine with 12 cores)





