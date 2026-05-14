module PseudoSpectral
export pseudospectral_problem, x

using SciMLBase: SplitODEProblem, ODEProblem, DiagonalOperator, ODEFunction, update_coefficients!, remake
using FFTW: plan_r2r!, REDFT00, MEASURE, r2rFFTWPlan
using Symbolics: variable, @variables, Num, sparsejacobian, build_function, substitute, get_variables
using Random:seed!

const x = variable(:x) |> Num

 # BC is Nothing for homogeneous for BC or Function for Heterogeneous BC.
 struct PseudoSpectralProblem
    problem::SplitODEProblem
    dims::Tuple{Int,Int}
    species::Vector{Num}
    reaction_params::Vector{Num}
    diffusion_params::Vector{Num}
    boundary_params::Vector{Num}
    initial_params::Vector{Num}
    plan::r2rFFTWPlan
    lifting_function::Union{Nothing, Function}
end

struct PseudoSpectralSolution
    sol::ODESolution
    species::Vector{Num}
    plan::r2rFFTWPlan
end

function make_lifting_function(boundary_conditions, diffusion_rates, boundary_params,diffusion_params, n)
    iszero(boundary_conditions) && return nothing
    a,b = eachrow(boundary_conditions)
    X = range(0.0,1.0,n)
    ϕ = X.^2 * (b'-a')/2 + X * a'
    Δϕ = ((b-a).*diffusion_rates)' |> collect
    fϕ,_ = build_function(ϕ, bs; expression=Val{false})
    fΔϕ,_ = build_function(Δϕ, ds, bs; expression=Val{false})
    (d, b) -> (fϕ(b), fΔϕ(d,b))
end

function make_initial_function(initial_conditions, initial_params, initial_noise, n; seed=nothing)
    u0 = [substitute(ic, x=>X) for X in range(0,1,n), ic in initial_conditions]
    f,_= build_function(u0, initial_params; expression=Val{false})
    initial_noise = initial_noise * abs.(randn(n,m))
    function (p)
        seed!(seed)
        noise = initial_noise * abs.(randn(n,m))
        f(p) + noise
    end
end


"""
Construct a SplitODEProblem to solve a reaction diffusion system with reflective boundaries.

Returns the SplitODEProblem with solutions in the frequency (DCT-1) domain and a FFTW plan to transform solutions back to the spatial domain.
"""
function PseudoSpectralProblem(species, reaction_rates, diffusion_rates, boundary_conditions, initial_conditions, num_verts; noise=1e-4, seed=nothing, kwargs...)
    n = num_verts
    m = length(species)
    
    # Collect parameter symbols. 
    rs, ds, bs, is = (setdiff(collect_variables(exprs), x, species) for exprs in (reaction_rates, diffusion_rates, vec(boundary_conditions), initial_conditions))

    u = Matrix{Float64}(undef, n, m)
    plan = 1/sqrt(2*(n-1)) * plan_r2r!(u, REDFT00, 1; flags=MEASURE)

    fu0 = make_initial_function(initial_conditions, is, noise, n; seed=seed)
    lf = make_lifting_functions(boundary_conditions, diffusion_rates, bs,ds, n)
    
    R = reaction_operator(species, reaction_rates, rs, plan, Val(!isnothing(lf)))
    D = diffusion_operator(diffusion_rates, ds, n)
    prob = SplitODEProblem(D, R, vec(u), Inf, nothing; kwargs...)
    PseudoSpectralProblem(prob, (n,m), species, rs, ds, bs, is, plan, lf)
end

function remake!(prob::PseudoSpectralProblem; p=nothing, kwargs...)
    if isnothing(p)
        prob.prob = remake(prob.prob; kwargs...)
        return prob
    end
    r = Float64[params[k] for k in prob.rs]
    d = Float64[params[k] for k in prob.ds]
    b = Float64[params[k] for k in prob.bs]
    i = Float64[params[k] for k in prob.is]

    w = Matrix{Float64}(undef,n,m) # Allocate working memory for FFTW.
    u0 = prob.fu0(i)
    lf = prob.lifting_function
    if !isnothing(lf)
        ϕ, Δϕ = lf(d,b)
        u0 .-= ϕ
    else
        ϕ = Δϕ = Matrix{Float64}(undef,0,0)
    end

    p = Parameters(w,r,d,ϕ,Δϕ,attempt)

    plan * u0
    u0 = vec(u0)
    update_coefficients!(prob.prob.f.f1.f, nothing, p, nothing) # Set parameter values in diffusion operator.
    prob.prob = remake(prob.prob; p, kwargs...) # Set parameter values in SplitODEProblem.
    prob
end

function solve(prob::PseudoSpectralProblem, alg; kwargs...)
    sol = solve(prob, alg; kwargs)
    PseudoSpectralSolution(sol, prob.species, prob.plan)
end
    

function transform(sol, u)
    dims = size(sol.plan)
    u = reshape(u, dims)
    plan * u
    if !isnothing(sol.ϕ)
        u .+= sol.ϕ
    end
    u
end

getindex(sol::PseudoSpectralSolution, i::Union{Int,CartesianIndex{1}}) = transform(sol, sol.sol.u[i])
(sol::PseudoSpectralSolution)(t) = transform(sol, sol.sol(t))
get_u(sol::PseudoSpectralSolution) = stack(sol[i] for in in eachindex(sol))
eachindex(sol::PseudoSpectralSolution) = eachindex(sol.sol)



"Build function for the reaction component, with `f(v+ϕ) + Δϕ` offset for non-zero-flux BCs."
function reaction_operator(species, reaction_rates, rs, plan!, ::Val{BC})
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

    λ = vec(σ² * diffusion_rates')
    (f,f!) = build_function(λ, ps; expression=Val{false})
    λ0 = similar(λ, Float64)
    update!(λ,u,p,t) = f!(λ, p.d)
    DiagonalOperator(λ0; update_func! = update!)
end

struct Parameters{T}
    u :: Matrix{Float64} # Working array for dct.
    r :: Vector{Float64} # Reaction parameters.
    d :: Vector{Float64} # Diffusion parameters.
    ϕ :: Matrix{Float64} # Boundary lifting function
    Δϕ :: Matrix{Float64}
    attempt :: Int64 # Track number of attempts at solution.
end

## Symbolics utility functions
"Sort parameters by name."
sort_variables(p) = sort(p, by=_nameof)
_nameof(v) = isspecies(v) ? nameof(v.f) : nameof(v)


"Extract variables from a (possibly nested) collection of expressions and sort them by name."
collect_variables(exprs...) = collect_variables(exprs) # Combine multiple arguments.
collect_variables(exprs::Union{Tuple,Vector}) = exprs .|> collect_variables |> splat(union) |> sort_variables
collect_variables(expr) = get_variables(expr) |> collect # Call recursively until we get down to a single expression.

end

