module Simulate
export simulate

using PseudoSpectralReactionDiffusion
using ..Models
using ..Util: issingle
import SciMLBase
using SciMLBase.ReturnCode: Default
using OrdinaryDiffEqExponentialRK: ETDRK4
# using OrdinaryDiffEqTsit5: Tsit5 # temp
using OrdinaryDiffEqSDIRK: KenCarp3

using ProgressMeter: Progress, BarGlyphs, next!

using Symbolics:Num #temp

using Logging: with_logger, ConsoleLogger, stderr, Error
using Pipe: @pipe
using Random: seed!
"""
    simulate(model, params; output_func=nothing, full_solution=false, alg=ETDRK4(), num_verts=64, dt=0.1, max_attempts = 4, tol=1e-4, kwargs...)


Simulate `model` for the parameters and initial conditions given in `params`, stopping when a steady state is reached. Returns `(u,t)` with the solution values and time.


# Arguments
- `model`: `Model` object containing the system to be simulated.
- `params`: Either a single parameter set or a vector of parameter sets to be solved as an ensemble. Parameter sets can be created manually with parameter_set or supplied as a dict or collection of pairs in which case defaults will be used for any missed values and low-level noise added to initial conditions. Parameters values may be either single numbers which are replicated homogenously over the domain, or functions mapping the interval [0.0,1.0] to values for the corresponding point in space. 
- `output_func(u, t)`: Function to transform output values. 
- `full_solution`: Return a vector of values at each time point if true, instead of just the steady-state solution.
- `max_attempts`: Number of times to retry with reduced dt before giving up if the solution fails to converge. 
- `num_verts`: Number of points in spatial discretisation.
For other keyword arguments see https://docs.sciml.ai/DiffEqDocs/stable/basics/common_solver_opts/.
"""
simulate(model, params; kwargs...) = simulate(model; kwargs...)(params)

"""
    function simulate(model; output_func=nothing, full_solution=false, alg=ETDRK4(), num_verts=64, dt=0.1, max_attempts = 4, tol=1e-5, noise=1e-4, kwargs...)

Partially applied version of `simulate` to avoid repeating expensive setup when simulating the same model reapeatedly.
"""
function simulate(model; discretisation=:pseudospectral, seed=nothing, kwargs...)
    seed!(seed)
    if discretisation==:pseudospectral
        simulate_pseudospectral(model; kwargs...)
    elseif discretisation==:mol
        simulate_mol(model; kwargs...)
    else
        error("Unsupported discretisation: $(discretisation).")
    end
end

function simulate_pseudospectral(model; output_func=nothing, alg=ETDRK4(), num_verts=64, dt=0.1, max_attempts = 4, tol=1e-5, noise=1e-4, kwargs...)
    prob = PseudoSpectralProblem(model, num_verts; noise=noise)

    f(params) = f([params]).u |> only # Accept a single parameter set instead of a vector.
    f(params::AbstractVector) = f(parameter_set.(model, params))
    function f(params::Vector{ParameterSet})
        isempty(params) && return PseudoSpectralSolution(species(model),[], [], Default) # Handle an empty collection of parameter sets.
        
        progress = Progress(length(params); desc="Simulating parameter sets: ", dt=0.1, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow)

        function _output_func(sol,ctx)
            if successful_retcode(sol)
                out = isnothing(output_func) ? sol : output_func(sol)
                repeat = false
                next!(progress) # Advance progress bar.
            else
                out = missing
                repeat = ctx.repeat < max_attempts
            end
            (out, repeat)
        end
            
        function _prob_func(prob, ctx)
            p = params[ctx.sim_id]
            _dt = dt/2^(ctx.repeat-1) # halve dt if solve was unsuccessful.
            @warn "retry with dt <- dt/2"
            prob = remake(prob; p, dt=_dt, rng=ctx.rng)
        end

        ensemble_prob = EnsembleProblem(prob; prob_func=_prob_func, output_func=_output_func)
        
        # with_logger(ConsoleLogger(stderr, Error)) do
        solve(ensemble_prob, alg; trajectories=length(params), callback=steady_state_callback(tol), kwargs...)
        # end
    end
end


function simulate_mol(model; output_func=tuple, full_solution=false, alg=KenCarp3(), num_verts=64, tol=1e-5, noise=1e-4, kwargs...)
    prob = mol_problem(model, num_verts)
    prob = remake(prob; u0=prob.u0+noise*abs.(randn(length(prob.u0))))
    sps = species(model)
    f(params) = f([params]) |> only # Accept a single parameter set instead of a vector.
    f(params::AbstractVector) = f(parameter_set.(model, params))
    function f(params::Vector{ParameterSet})
        isempty(params) && return SciML.EnsembleSolution([], 0.0, false) # Handle an empty collection of parameter sets.

        progress = Progress(length(params); desc="Simulating parameter sets: ", dt=0.1, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow)
        function _output_func(sol,i)
            u = @pipe sol.u |> values |> stack |> permutedims(_,[2,3,1])
            t = sol.t
            if !full_solution
                u = u[:,:,end]
                t = t[end]
            end
            if successful_retcode(sol)
                out = output_func(u,t)
                next!(progress) # Advance progress bar.
            else
                out = missing
            end
            (out, false)
        end
        
        prob_func(prob, i, attempt) = SciML.remake(prob; p=params[i])
        ensemble_prob = SciML.EnsembleProblem(prob; output_func=_output_func, prob_func=prob_func)
        
        with_logger(ConsoleLogger(stderr, Error)) do
            SciML.solve(ensemble_prob, alg; trajectories=length(params), callback=steady_state_callback(tol), verbose=false, maxiters=1e6, kwargs...)
        end
    end
end

end


