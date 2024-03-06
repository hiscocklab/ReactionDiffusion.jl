**NOTE: change the logo.png and favicon.ico images.** Currently these are from Catalyst.jl. Design a new logo - idea: 3 julia colors as nodes of a reaction-diffusion network?

**NOTE: one BIG think we need to do before starting the documentation**
We need to work out how to handle plotting.
We have two options
1. write a separate package `ReactionDiffusionPlots` using `Makie.jl`. **PREFERRED**
2. somehow copy and paste the helper plotting scripts that you use, along with `Plots.jl`, and just say on the documentation "write your own plotting scripts"
    - I think it's helpful to fully document the simulation outputs so people know what they are and how they would plot it

The advantage of number 2 is that it is more customizable. The disadvantage is that it is more complicated to explain! It would be nice if we could just `add ReactionDiffusionPlots` and it worked. However, `Plots.jl` is notoriously fiddly with dependencies - I have not managed to get it working as a package dependency. `Makie.jl` is the future! So we would need to convert our helper functions into `Makie.jl` and make all the plots using these helper functions. I think this wouldn't be too difficult, but we can discuss.

# ReactionDiffusion.jl for modelling biological systems

- insert a short, punchy description of the package
    - easy to use
    - very fast
    - intuitive and biological
- provide a trimmed down illustration of the package functionality, using a biologically-inspired case study as an example.
    - how about a gdf5_mRNA, GDF5, nog_mRNA, NOG, pSMAD model, motivated from Grall et al. 2024?
    - we should get across:
        - intuitive @reaction_network model specification
        - single line of code screens millions of parameters
        - single line of code simulates the PDE
        - show images and a movie!
    - in this introductory section, we do not provide complete code - just snippets to emphasize how straightforward it is to use. 



## Quick start guide

Here we should provide as simple as possible example of using the scripts. This should be the same example as above, but now complete end-to-end code. Most users may just look at the quick start guide - it should be enough to reproduce a simulation, but not too much that it sounds too complicated!

## Citation and future developments

- provide citation to our preprint
    - also to Catalyst.jl and DiffEq.jl
- say that more functionality is planned in the future, but we are currently a small team of developers




