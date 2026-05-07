module PseudoSpectral
export pseudospectral_problem, x

using ..Util: collect_variables, safe_stack
using SciMLBase: SplitODEProblem, ODEProblem, DiagonalOperator, ODEFunction, update_coefficients!, remake
using FFTW: plan_r2r!, REDFT00, MEASURE
using Symbolics: @variables, sparsejacobian, build_function, substitute
using Random:seed!

const x = only(@variables(x))

"""
Construct a SplitODEProblem to solve a reaction diffusion system with reflective boundaries.

Returns the SplitODEProblem with solutions in the frequency (DCT-1) domain and a FFTW plan to transform solutions back to the spatial domain.
"""
function pseudospectral_problem(species, reaction_rates, diffusion_rates, boundary_conditions, initial_conditions, num_verts; noise=1e-4, seed=nothing, kwargs...)
    n = num_verts
    m = length(species)
    
    # Collect parameter symbols. 
    rs, ds, bs, is = (setdiff(collect_variables(exprs), x, species) for exprs in (reaction_rates, diffusion_rates, vec(boundary_conditions), initial_conditions))

    u = Matrix{Float64}(undef, n, m)
    plan = 1/sqrt(2*(n-1)) * plan_r2r!(u, REDFT00, 1; flags=MEASURE)


    u0 = [substitute(ic, x=>X) for X in range(0,1,n), ic in initial_conditions]
    _fu0,_= build_function(u0, is; expression=Val{false})

    seed!(seed)
    initial_noise = noise * abs.(randn(n,m))

    fu0(i) = _fu0(i) + initial_noise

    # Dispatch on `iszero(boundary_conditions)`
    function dispatch_bcs(::Val{BC}) where BC
        ## Offsets for constant, non-zero flux boundary conditions.
        # For u′(0) = a, u′(1) = b,
        # define ϕ as a smooth function such that ϕ′(0) = a, ϕ′(1) = b, and write v = u - ϕ.
        # Then v′(0) = 0, v′(1) = 0, so we can solve for v using DCT-I.
        if BC
            a,b = eachrow(boundary_conditions)
            X = range(0.0,1.0,n)
            ϕ = X.^2 * (b'-a')/2 + X * a'
            Δϕ = ((b-a).*diffusion_rates)' |> collect
            fϕ,_ = build_function(ϕ, bs; expression=Val{false})
            fΔϕ,_ = build_function(Δϕ, ds,bs; expression=Val{false})
        end
        
        R = reaction_operator(species, reaction_rates, rs, plan, Val(BC))
        D = diffusion_operator(diffusion_rates, ds, n)
        prob = SplitODEProblem(D, R, vec(u), Inf, nothing; kwargs...)

        function make_problem(params; attempt=1, kwargs...)
            r = Float64[params[k] for k in rs]
            d = Float64[params[k] for k in ds]
            b = Float64[params[k] for k in bs]
            i = Float64[params[k] for k in is]

            local ϕ, Δϕ
            if BC
                ϕ = fϕ(b)
                Δϕ = fΔϕ(d,b)
            else
                ϕ = Δϕ = Matrix{Float64}(undef,n,m)
            end

            w = Matrix{Float64}(undef,n,m) # Allocate working memory for FFTW.

            p = Parameters(w,r,d,ϕ,Δϕ,attempt)

            local u0
            u0 = fu0(i)
            BC && (u0 .-= p.ϕ)
            plan * u0
            u0 = vec(u0)
            
            update_coefficients!(prob.f.f1.f, nothing, p, nothing) # Set parameter values in diffusion operator.
            remake(prob; p, u0, kwargs...) # Set parameter values in SplitODEProblem.
        end

        # Function to transform output back to spatial domain.
        # TODO fix code duplication.
        function transform(sol, t::Float64)
            u = sol(t)
            u = reshape(u,n,m)
            plan * u
            BC && (u .+= sol.prob.p.ϕ)
            u
        end

        function transform(sol, i::Int)
            u = sol[i]
            u = reshape(u,n,m)
            plan .* u
            BC && (u .+= sol.prob.p.ϕ)
            u
        end

        transform(sol) = transform.(sol, eachindex(sol))
        

        make_problem, transform
    end
    
    dispatch_bcs(Val(!iszero(boundary_conditions)))
end

"Build function for the reaction component, with `f(v+ϕ) + Δϕ` offset for non-zero-flux BCs."
function reaction_operator(species, reaction_rates, rs, plan!, ::Val{BC}) where BC
    n,m = size(plan!)
    @variables u[1:n, 1:m]
    # TODO: Clever things to make only spatially varying parameters expand?
    # Build an nxm matrix of derivatives, substituting reactants for u[i,j] and parameters for p[k,l].
    du = [substitute(expr, Dict([x=>X, zip(species,v)...])) for (v,X) in zip(eachrow(u), range(0,1,n)), expr in reaction_rates]
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

struct Parameters
    u :: Matrix{Float64} # Working array for dct.
    r :: Vector{Float64} # Reaction parameters.
    d :: Vector{Float64} # Diffusion parameters.
    ϕ :: Matrix{Float64} # Boundary lifting function
    Δϕ :: Matrix{Float64}
    attempt :: Int64 # Track number of attempts at solution.
end

end

