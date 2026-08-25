module Optimise
export optimise, plot_cost

using ..Util: zip_dict, tmap
using ..Simulate
using FiniteDiff
using StatsBase: sample
using Pipe: @pipe
using Optim: optimize, SAMIN, Options


function optimise(model, cost, vars, params_min, params_max, params0; in_domain=x->true, sample=nothing, max_steps=10000, verbosity=1, kwargs...)
    _simulate = simulate(model; kwargs...)
    function _cost(sol)
        (isempty(sol) || any(ismissing, sol.u)) && return 1.0
        isnothing(sample) ? cost(only(sol.u)...) : cost(sol.u)
    end

    _sample = something(sample, x->[x])

    __cost(p) = @pipe p |> zip_dict(vars,_) |> merge(params0, _) |> _sample |> filter(in_domain,_) |> _simulate |> _cost

    p_min = [params_min[v] for v in vars]
    p_max = [params_max[v] for v in vars]
    p0 = [params0[v] for v in vars]
    callback(state) = state.value <= 0.1
    optimize(__cost, p_min, p_max, p0, SAMIN(verbosity=verbosity), Options(iterations=max_steps, callback=s -> s.f_x <= 0.0))
end


function plot_cost(model, cost, vars, params; in_domain=x->true, sample=nothing, kwargs...)
    _simulate = simulate(model; kwargs...)
    
    function _cost(sol)
        isempty(sol.u) || any(ismissing, sol.u) && return Inf
        isnothing(sample) ? cost(only(sol.u)...) : cost(sol.u)
    end

    _sample = something(sample, x->[x])

    __cost(p) = @pipe p |> zip_dict(vars,_) |> merge(params[1], _) |> _sample |> filter(in_domain,_) |> _simulate |> _cost

    p = [[params[v] for v in vars] for params in params]
    c = tmap(__cost, Float64, p)
    a = [first(x) for x in p]
    b = [last(x) for x in p] # TODO Proper unzip
    (a,b,c)
end

end