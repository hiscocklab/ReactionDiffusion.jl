module Turing
export turing_wavelength, is_turing, filter_turing, turing_wavelength_problem

using ..Models
using ..Util: issingle, tmap, tfilter

using Groebner
using Symbolics: jacobian, symbolic_solve, substitute, build_function, Num
using SciMLBase: remake, SteadyStateProblem, solve
using SteadyStateDiffEq: SSRootfind
using LinearAlgebra: diagm, eigvals
using Base.Threads: @threads


"""
    turing_wavelength(model, params; k=logrange(0.01,1000,1000))

Compute dominant wavelengths of Turing instabilities for each of `params`.
Returns 0.0 for parameter sets for which Turing instability does not occur.
"""
turing_wavelength(model, params; k=logrange(0.01,1000,1000)) = turing_wavelength(model;k=k)(params)

function turing_wavelength(model; k=logrange(0.01,1000,1000))   
    du = reaction_rates(model)
    u = species(model) .|> Num # For some unfathomable reason, symbolic_solve complains about missing Groebner if we don't convert u to Num.
    ps = parameters(model)
    jac = jacobian(du,u; simplify=true)
    d = diffusion_rates(model)/domain_size(model)^2 |> collect
    (fd,fd!) = build_function(diagm(d), ps; expression=Val{false})
    k² = k.^2
    
    ss = symbolic_solve(du, u; warns=false)
    
    if isnothing(ss)
        # Fall back to numerical solution.
        @warn "Symbolics.symbolic_solve can not solve this input currently. Falling back to finding steady state solutions numerically."
        (f,f!) = build_function(du, u, ps, (); expression=Val{false})
        (fjac,fjac!) = build_function(jac, u, ps; expression=Val{false})
        u0 = zeros(num_species(model))
        ss_prob = SteadyStateProblem(f!, u0)
        function (params)
            params = parameter_set(model, params)
            p = [params[key] for key in ps]
            prob = remake(ss_prob; p=p)
            ss = solve(prob, SSRootfind())
            J = fjac(ss.u, p)
            
            all(<(eps(J[1])), real(eigvals(J))) || return 0.0
            D = fd(p)
            real_max, i = findmax(real(eigvals(J - D * k²)[end]) for k² in k²)
            real_max > 0.0 ? 2pi/k[i] : 0.0
        end
    else
        jac_ss = substitute(jac, only(ss)) # TODO: Handle multiple steady states.
        (fjac,fjac!) = build_function(jac_ss, ps; expression=Val{false})
        function (params)
            params = parameter_set(model, params)
            p = [params[key] for key in ps]
            J = fjac(p)
            all(<(eps(J[1])), real(eigvals(J))) || return 0.0
            D = fd(p)
            real_max, i = findmax(real(eigvals(J - D * k²)[end]) for k² in k²)
            real_max > 0.0 ? 2pi/k[i] : 0.0
        end
    end
end

"""
    is_turing(model,params)

Test whether `params` exhibit Turing instability.
"""
is_turing(model, params) = is_turing(model)(params)

function is_turing(model)
    f = turing_wavelength(model)
    params -> any(>(0.0), f(params))
end
"""
    filter_turing(model,params)

Return only `params` which demonstrate Turing instability. Multithreaded.
"""
function filter_turing(model,params)
    f = turing_wavelength(model)
    tfilter(params) do p
        f(p) > 0.0
    end
end

end