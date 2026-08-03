module PseudoSpectral
export PseudoSpectralProblem, PseudoSpectralSolution, PseudoSpectralIntegrator, steady_state_callback, x, step!, step_to!, get_sol, get_u

import Base: getindex, eachindex, lastindex
export getindex, eachindex, lastindex

import SciMLBase: EnsembleProblem, solve, remake, successful_retcode, DEIntegrator
export EnsembleProblem, solve, remake, successful_retcode

import SciMLBase
using SciMLBase: SplitODEProblem, ODEProblem, ODESolution, DiagonalOperator, ODEFunction, update_coefficients!, ReturnCode, DiscreteCallback, terminate!, get_du, init
using SciMLBase.ReturnCode: Terminated
using OrdinaryDiffEqExponentialRK: ETDRK4
using FFTW: plan_r2r!, REDFT00, MEASURE, ScaledPlan
using Symbolics: variable, @variables, Num, sparsejacobian, build_function, substitute, get_variables
using Random: default_rng, AbstractRNG

const x = variable(:x) |> Num

 # BC is Nothing for homogeneous for BC or Function for Heterogeneous BC.
 mutable struct PseudoSpectralProblem
    ode_problem::ODEProblem
    dims::Tuple{Int,Int}
    species::Vector{Num}
    reaction_params::Vector{Num}
    diffusion_params::Vector{Num}
    boundary_params::Vector{Num}
    initial_params::Vector{Num}
    plan::ScaledPlan
    initial_function::Function
    lifting_function::Union{Nothing, Function}
    rng::AbstractRNG
end

struct PseudoSpectralSolution
    species::Vector{Num}
    u::Vector{Matrix{Float64}}
    x:: Vector{Float64}
    t::Vector{Float64}
    retcode::ReturnCode.T
end

struct Parameters
    u :: Matrix{Float64} # Working array for dct.
    r :: Vector{Float64} # Reaction parameters.
    d :: Vector{Float64} # Diffusion parameters.
    ϕ :: Matrix{Float64} # Boundary lifting function
    Δϕ :: Matrix{Float64}
    attempt :: Int64 # Track number of attempts at solution. TODO: remove this.
end


"""
Construct a SplitODEProblem to solve a reaction diffusion system with reflective boundaries.

Returns the SplitODEProblem with solutions in the frequency (DCT-1) domain and a FFTW plan to transform solutions back to the spatial domain.
"""
function PseudoSpectralProblem(species, reaction_rates, diffusion_rates, boundary_conditions, initial_conditions, num_verts; p=nothing, noise=1e-4, rng=default_rng(), kwargs...)
    n = num_verts
    m = length(species)
    
    # Collect parameter symbols. 
    rs, ds, bs, is = (setdiff(collect_variables(exprs), x, species) for exprs in (reaction_rates, diffusion_rates, vec(boundary_conditions), initial_conditions))
    p = something(p, Dict(q => 0 for q in union(rs,ds,bs,is)))

    u = Matrix{Float64}(undef, n, m)
    plan = 1/sqrt(2*(n-1)) * plan_r2r!(u, REDFT00, 1; flags=MEASURE)

    fu0 = make_initial_function(initial_conditions, is, noise, n)
    lf = make_lifting_function(boundary_conditions, diffusion_rates, bs,ds, n)
    
    R = reaction_operator(species, reaction_rates, rs, plan, Val(!isnothing(lf)))
    D = diffusion_operator(diffusion_rates, ds, n)
    odeprob = SplitODEProblem(D, R, vec(u), Inf, nothing; kwargs...)
    prob = PseudoSpectralProblem(odeprob, (n,m), species, rs, ds, bs, is, plan, fu0, lf, rng)
    remake(prob; p)
end


function remake(prob::PseudoSpectralProblem; p=nothing, rng=nothing, attempt=1, kwargs...)
    if isnothing(p)
        prob.ode_problem = remake(prob.ode_problem; kwargs...)
        return prob
    end
    if !isnothing(rng)
        prob.rng=rng
    end
    r = Float64[p[k] for k in prob.reaction_params]
    d = Float64[p[k] for k in prob.diffusion_params]
    b = Float64[p[k] for k in prob.boundary_params]
    i = Float64[p[k] for k in prob.initial_params]

    w = Matrix{Float64}(undef,prob.dims...) # Allocate working memory for FFTW.
    u0 = prob.initial_function(i,prob.rng)
    lf = prob.lifting_function
    if !isnothing(lf)
        ϕ, Δϕ = lf(d,b)
        u0 .-= ϕ
    else
        ϕ = Δϕ = Matrix{Float64}(undef,0,0)
    end
    p = Parameters(w,r,d,ϕ,Δϕ,attempt)
    prob.plan * u0
    u0 = vec(u0)
    update_coefficients!(prob.ode_problem.f.f1.f, nothing, p, nothing) # Set parameter values in diffusion operator.
    prob.ode_problem = remake(prob.ode_problem; u0, p, kwargs...) # Set parameter values in SplitODEProblem.
    prob
end


function solve(prob::PseudoSpectralProblem, alg; kwargs...)
    odesol = SciMLBase.solve(prob.ode_problem, alg; kwargs...)
    PseudoSpectralSolution(prob, odesol)
end

# Separate constructor so we can use it both with solve and as an output function for EnsmbleProblem.
function PseudoSpectralSolution(prob::PseudoSpectralProblem, sol::ODESolution)
    u = [transform(prob, u) for u in sol.u]
    PseudoSpectralSolution(prob.species, u, range(0.0,1.0,prob.dims[1]), sol.t, sol.retcode)
end

function transform(prob::PseudoSpectralProblem, u)
    u = reshape(u, prob.dims)
    prob.plan * u
    if !isnothing(prob.lifting_function)
        u .+= prob.ode_problem.p.ϕ
    end
    u
end


function make_lifting_function(boundary_conditions, diffusion_rates, boundary_params,diffusion_params, n)
    iszero(boundary_conditions) && return nothing
    a,b = eachrow(boundary_conditions)
    X = range(0.0,1.0,n)
    ϕ = X.^2 * (b'-a')/2 + X * a'
    Δϕ = ((b-a).*diffusion_rates)'
    fϕ,_ = build_function(ϕ, boundary_params; expression=Val{false})
    fΔϕ,_ = build_function(Δϕ, diffusion_params, boundary_params; expression=Val{false})
    (d, b) -> (fϕ(b), fΔϕ(d,b))
end


function make_initial_function(initial_conditions, initial_params, initial_noise, n)
    m = length(initial_conditions)
    u0 = [substitute(ic, x=>X) for X in range(0,1,n), ic in initial_conditions]
    f,_= build_function(u0, initial_params; expression=Val{false})
    function (p,rng)
        noise = initial_noise * abs.(randn(rng,n,m))
        f(p) + noise
    end
end


"Build function for the reaction component, with `f(v+ϕ) + Δϕ` offset for non-zero-flux BCs."
function reaction_operator(species, reaction_rates, rs, plan!, ::Val{BC}) where BC
    n,m = size(plan!)
    @variables u[1:n, 1:m]
    # TODO: Clever things to make only spatially varying parameters expand?
    # Build an nxm matrix of derivatives, substituting reactants for u[i,j] and parameters for p[k,l].
    du = [substitute(expr, Dict([x=>X, zip(species,v)...])) for (v,X) in zip(eachrow(collect(u)), range(0,1,n)), expr in reaction_rates]
    _, f! = build_function(du, u, rs; expression=Val{false})
    
    function f̂!(du,u,p,t)
        du = reshape(du,n,m)
        copyto!(p.u, u)
        plan! * p.u
        BC && (p.u .+= p.ϕ)
        f!(du, p.u, p.r)
        BC && (du .+= p.Δϕ)
        plan! * du
        nothing
    end
    ODEFunction(f̂!)
end

"Build linear operator for the diffusion component."
function diffusion_operator(diffusion_rates, ps, n)
    k = 0:n-1 # Wavenumbers
    h = 1/(n-1)

    ## 2nd order Fourier differentiation coefficients.
    # For a continuous FT this would be σ² = -(k/2pi h)^2, but corrected for
    # the discrete transform this becomes:
    σ² = @. -(4/h^2) * sin(k*pi/(2*(n-1)))^2

    λ = vec(σ² * diffusion_rates') |> collect
    (f,f!) = build_function(λ, ps; expression=Val{false})
    λ0 = similar(λ, Float64)
    update!(λ,u,p,t) = f!(λ, p.d)
    DiagonalOperator(λ0; update_func! = update!)
end


function getindex(sol::PseudoSpectralSolution, species::Num)
    name=nameof(species)
    i = findfirst(s -> nameof(s.val.f)===name, sol.species)
    [vec(u[:,i]) for u in sol.u]
end

getindex(sol::PseudoSpectralSolution) = getindex(sol.u)
getindex(sol::PseudoSpectralSolution, i::Union{Int,CartesianIndex{1}}) = sol.u[i]

eachindex(sol::PseudoSpectralSolution) = eachindex(sol.u)
# (sol::PseudoSpectralSolution)(t) = transform(sol, sol.sol(t))

lastindex(sol::PseudoSpectralSolution) = lastindex(sol.u)

successful_retcode(sol::PseudoSpectralSolution) = SciMLBase.successful_retcode(sol.retcode)

function EnsembleProblem(prob::PseudoSpectralProblem, params; output_func=nothing, kwargs...)
    prob_func(_prob,ctx) = remake(_prob; p=params[ctx.sim_id], rng=ctx.rng)
    EnsembleProblem(prob; prob_func,output_func, trajectories=length(params), kwargs...)
end

function EnsembleProblem(prob::PseudoSpectralProblem; prob_func, output_func=nothing, kwargs...)
    _prob_func(_prob, ctx) = prob_func(prob, ctx).ode_problem
    function _output_func(sol, ctx) 
        ps_sol = PseudoSpectralSolution(prob, sol)
        isnothing(output_func) ?  (ps_sol,false) : output_func(ps_sol,ctx)
    end
    SciMLBase.EnsembleProblem(prob.ode_problem; prob_func=_prob_func, output_func=_output_func, kwargs...)
end


## Integrator interface
mutable struct PseudoSpectralIntegrator
    integrator::DEIntegrator
    prob::PseudoSpectralProblem
    ss::Float64
end

function PseudoSpectralIntegrator(prob::PseudoSpectralProblem; alg=ETDRK4(), kwargs...)
    integrator=init(prob.ode_problem, alg; kwargs...)
    PseudoSpectralIntegrator(integrator, prob, Inf)
end

get_sol(integrator::PseudoSpectralIntegrator) = PseudoSpectralSolution(integrator.prob, integrator.integrator.sol)

function get_u(integrator::PseudoSpectralIntegrator, t) 
    step_to!(integrator, t)
    u = integrator.integrator.sol(t)
    transform(integrator.prob, u)
end

function step!(integrator::PseudoSpectralIntegrator, dt=nothing, stop_at_tdt=false)
    SciMLBase.step!(integrator.integrator, dt, stop_at_tdt)
    if integrator.integrator.sol.retcode == Terminated
        integrator.ss = integrator.integrator.t
    end
end

function step_to!(integrator::PseudoSpectralIntegrator, t, stop_at_tdt=false)
    dt = max(0.0, t - integrator.integrator.t)
    step!(integrator, dt, stop_at_tdt)
end

function remake(integrator::PseudoSpectralIntegrator; kwargs...)
    prob = remake(integrator.prob; kwargs...)
    PseudoSpectralIntegrator(prob)
end

function steady_state_callback(tol=1e-4)
    condition(u,t,integrator) = isapprox(get_du(integrator), zero(u); atol=tol)
    DiscreteCallback(condition, terminate!)
end


## Symbolics utility functions
"Sort parameters by name."
sort_variables(p) = sort(p, by=_nameof)
#_nameof(v) = isspecies(v) ? nameof(v.f) : nameof(v)
_nameof(v) = try nameof(v); catch e nameof(v.f) end # TODO: Something less hacky.


"Extract variables from a (possibly nested) collection of expressions and sort them by name."
collect_variables(exprs...) = collect_variables(exprs) # Combine multiple arguments.
collect_variables(exprs::Union{Tuple,Vector}) = exprs .|> collect_variables |> splat(union) |> sort_variables
collect_variables(expr) = get_variables(expr) |> collect # Call recursively until we get down to a single expression.

end

