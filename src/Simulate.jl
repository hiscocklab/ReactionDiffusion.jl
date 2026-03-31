module Simulate
export simulate

using ..Models
using ..PseudoSpectral
using ..Util: issingle
using SciMLBase: solve, successful_retcode, EnsembleProblem, EnsembleSolution, DiscreteCallback, terminate!, get_du
using OrdinaryDiffEqExponentialRK: ETDRK4
using ProgressMeter: Progress, BarGlyphs, next!

using Symbolics:Num #temp

using Logging: with_logger, ConsoleLogger, stderr, Error

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
function simulate(model; output_func=tuple, full_solution=false, alg=ETDRK4(), num_verts=64, dt=0.1, max_attempts = 4, tol=1e-5, noise=1e-4, kwargs...)

    make_prob, transform = pseudospectral_problem(model, num_verts; noise=noise)

    f(params) = f([params]) |> only # Accept a single parameter set instead of a vector.
    f(params::AbstractVector) = f(parameter_set.(model, params))
    function f(params::Vector{ParameterSet})
        isempty(params) && return EnsembleSolution([], 0.0, false) # Handle an empty collection of parameter sets.

        progress = Progress(length(params); desc="Simulating parameter sets: ", dt=0.1, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow)

        function _output_func(sol,i)
            attempt = sol.prob.p.attempt
            u,t = transform(sol; full_solution=full_solution)

            if successful_retcode(sol)
                out = output_func(u, t)
                repeat = false
                next!(progress) # Advance progress bar.
            else
                out = missing
                repeat = attempt < max_attempts
            end
            (out, repeat)
        end
            
        function prob_func(prob, i, attempt)
            p = params[i]
            dt′ = dt/2^(attempt-1) # halve dt if solve was unsuccessful.
            prob = make_prob(p, attempt; dt=dt′)
        end

        ensemble_prob = EnsembleProblem(make_prob(params[1]); output_func=_output_func, prob_func=prob_func)
        
        with_logger(ConsoleLogger(stderr, Error)) do
            solve(ensemble_prob, alg; trajectories=length(params), callback=steady_state_callback(tol), verbose=false, kwargs...)
        end
    end
end


function steady_state_callback(tol=1e-4)
    condition(u,t,integrator) = isapprox(get_du(integrator), zero(u); atol=tol)
    DiscreteCallback(condition, terminate!)
end

end


