# Simulate pattern formation in 1D

Having now screened through parameter sets, we can now simulate the corresponding PDEs. The simulation is performed on a 1D domain with reflective boundary conditions.  

To simulate the system for a single parameters set, simply run:

```julia
sol = simulate(model, turing_params[1])
```
This will return a solution object which can be indexed by species (`sol[NOG]`) to return a vector of concentrations for NOG at each time step, or by time index (`sol[end]`) to return a vector of num_verts x num_species matrices.

For instance if `species(model) == [GDF5, NOG, COMPLEX, pSMAD]`, then `sol[end][1,3]` is the concetration of `COMPLEX` at the left-most end of the domain at the conclusion of the simulation, and `sol.t[end]` will be the time at which this occurs. Alternatively the same value may be accessed to by species, as `sol[COMPLEX][end][end]`.


A collection of parameter sets can be simulated as an ensemble
```julia
sols = simulate(model, turing_params)
```
This will run in parallel and avoids repeating expensive setup computation for each parameter set. 

