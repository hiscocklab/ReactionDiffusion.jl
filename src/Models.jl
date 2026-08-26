module Models
export Model, name, species, parameters, reaction_parameters, boundary_parameters, diffusion_parameters,
    num_species, num_params, num_reaction_params, num_diffusion_params,
    domain_size, initial_conditions, noise,
    reaction_rates, diffusion_rates,
    @diffusion_system, @initial_conditions,
    parameter_set, ParameterSet,
    mol_problem

## Extended methods
import PseudoSpectralReactionDiffusion: PseudoSpectralProblem
export PseudoSpectralProblem

import ModelingToolkit: ODESystem
export ODESystem
##

using Symbolics: Num, value, get_variables, @variables, getname, substitute
import Catalyst # Catalyst.species and Catalyst.parameters would conflict with our functions.
using Catalyst: numspecies, numparams, assemble_oderhs, @species, @parameters, @reaction_network, ExprValues, get_usexpr, get_psexpr, esc_dollars!, find_parameters_in_rate!, forbidden_symbol_check, DEFAULT_IV_SYM, default_t, setmetadata, ReactionSystem, independent_variable
import ModelingToolkit # Needed for internal Catalyst functions.
using ..Util: subst, ensure_function, zip_dict
using Pipe


## MOL
using ModelingToolkit: Differential, PDESystem, @named
using SciMLBase: discretize
using MethodOfLines: MOLFiniteDifference#
##

# TODO CHECK for unnecessary Num conversions! Alternatively add needed Num conversions (and remove from Turing.jl)

# Types

# SpeciesValues
#___________________________________________________________________________________________________________________________________________________________________________________

SpeciesValues = Dict{Num,Num}

# DiffusionSystem
#___________________________________________________________________________________________________________________________________________________________________________________

struct DiffusionSystem
    domain_size::Num
    rates::SpeciesValues
end

"""
    @diffusion_system L begin D, species;... end

Define a spatial domain of length `L` and a set of diffusion rates and Neumann boundary conditions for the given species.
The boundary conditions `uₓ(0)=a` and `uₓ(L)=b` default to 0 if ommitted.

# Example
```
@diffusion_system L begin
    0.5,             U
    Dᵥ, (0.0, 0.5),  V
    Dᵣ/k, (a, a*s),  R
end
```
"""
macro diffusion_system(L, body)
    diffusion_system(L, body, __source__)
end

macro diffusion_system(body)
    diffusion_system(1, body, __source__)
end

function diffusion_system(L, body, source)
    species, parameters, pairs = parse_body(body, source)
    rexpr = dict_expr(pairs)
    L = parse_expr!(parameters, L)
    forbidden_symbol_check(parameters)
    psexpr = get_psexpr(parameters, [], Dict{Symbol, Expr}()) # @parameters
    iv = :($(DEFAULT_IV_SYM) = default_t()) # t
    sexpr = get_usexpr(species, Dict{Symbol, Expr}()) # @species
    dsexpr = :(DiffusionSystem($L, $rexpr))
    quote
        $psexpr
        $iv
        $sexpr
        $dsexpr
    end
end

parameters(ds::DiffusionSystem) = union(get_variables(ds.domain_size), parameters(ds.rates))
parameters(v::SpeciesValues) = @pipe v |> values .|> get_variables |> union(_..., []) |> Num.(_)

# Model
#___________________________________________________________________________________________________________________________________________________________________________________

"""
    Model(reaction, diffusion, boundary_conditions, initial_conditions)

An object containing a mathematical description of a reaction diffusion system to be simulated, independent of parameter values.

# Fields
- `reaction::ReactionSystem`
- `diffusion::DiffusionSystem`
- `boundary_flux::(ReactionSystem, ReactionSystem)`
- `initial_conditions::SpeciesValues`

"""
struct Model
    name
    reaction
    diffusion
    boundary_flux
    initial_conditions
end

Model(reaction, diffusion; name="") = Model(name, reaction, diffusion, (@reaction_network, @reaction_network), SpeciesValues())
Model(reaction, diffusion, initial::SpeciesValues; name="") = Model(name, reaction, diffusion, (@reaction_network, @reaction_network), initial)
Model(reaction, diffusion, boundary::Tuple{ReactionSystem, ReactionSystem}; name="") = Model(name, reaction, diffusion, boundary, SpeciesValues())
Model(reaction, diffusion, boundary::Tuple{ReactionSystem, ReactionSystem}, initial::SpeciesValues; name="") = Model(name, reaction, diffusion, boundary, initial)
# Don't try to broadcast over a model.
Base.broadcastable(model::Model) = Ref(model)

# ODESystem(model::Model)
#___________________________________________________________________________________________________________________________________________________________________________________

ODESystem(model::Model) = convert(ODESystem, model.reaction)

# Model getters
#___________________________________________________________________________________________________________________________________________________________________________________

# TODO Eliminate unused getters.
name(model::Model) = model.name
species(model::Model) = Catalyst.species(model.reaction)
parameters(model::Model) = union(reaction_parameters(model), diffusion_parameters(model), initial_condition_parameters(model), boundary_parameters(model))

reaction_parameters(model::Model) = Catalyst.parameters(model.reaction)
diffusion_parameters(model::Model) = parameters(model.diffusion)
initial_condition_parameters(model::Model) = parameters(model.initial_conditions)
boundary_parameters(model::Model) = union(Catalyst.parameters.(model.boundary_flux)...)

reaction_rates(model) = assemble_oderhs(model.reaction, species(model))
diffusion_rates(model::Model, default=0.0) = [get(model.diffusion.rates, s, default) for s in species(model)]
initial_conditions(model::Model, default=0.0) = [get(model.initial_conditions, s, default) for s in species(model)]

function boundary_flux(model::Model)
    b0, b1 = model.boundary_flux
    s = species(model)
    vcat(assemble_oderhs(b0, s)', -assemble_oderhs(b1, s)') .|> Num
end

num_species(model::Model) = numspecies(model.reaction)
num_params(model::Model) = num_reaction_params(model) + num_diffusion_params(model)
num_reaction_params(model::Model) = numparams(model.reaction)
num_diffusion_params(model::Model) = length(diffusion_parameters(model))

domain_size(model::Model) = model.diffusion.domain_size
function domain_size(model::Model, params)
    L = domain_size(model)
    L isa Num ? params[nameof(L)] : L
end

is_fixed_size(model::Model) = typeof(domain_size(model)) != Num # TODO use type system. 

reaction_parameters(model::Model, params, default=0.0) = subst(reaction_parameters(model), params, default)
#diffusion_parameters(model::Model, params, default=0.0) = get_vector(params, diffusion_parameters(model), default)

function diffusion_rates(model::Model, params::Dict{Symbol, Float64}, default=0.0) # wrong and bad
    syms = Dict(nameof(p) => p for p in parameters(model))
    params = Dict(syms[k] => v for (k, v) in params)
    [(substitute(D, params)) for D in diffusion_rates(model, default)]
end

# Initialisers
#___________________________________________________________________________________________________________________________________________________________________________________

"""
    @initial_conditions begin IC, species;... end

Define a set of initial conditions for the given species. IC may depend on arbitrary parameters and additionally the spatial variable `x`.

# Example
```
@inital_conditions begin
    U0,             U
    V0 + exp(x),    V
end
```
"""
macro initial_conditions(body)
    species, parameters, pairs = parse_body(body, __source__)
    icexpr = dict_expr(pairs)
    psexpr = get_psexpr(parameters, [], Dict{Symbol, Expr}()) # @parameters
    iv = :($(DEFAULT_IV_SYM) = default_t()) # t
    sexpr = get_usexpr(species, Dict{Symbol, Expr}()) # @species
    quote
        $psexpr
        $iv
        $sexpr
        $icexpr
    end
end

# parameters
#___________________________________________________________________________________________________________________________________________________________________________________

ParameterSet = Dict{Num, Float64}

"""
    function parameter_set(model, params)

Create a set of parameter values and initial conditions for `model`.
Defaults are used for values missing from `params`.
"""
function parameter_set(model, params=Dict())
    set = ParameterSet()
    
    for rs in reaction_parameters(model)
        set[rs] = get(params, nameof(rs), 1.0)
    end

    for ds in diffusion_parameters(model)
        set[ds] = get(params, nameof(ds), 0.0)
    end

    for ds in boundary_parameters(model)
        set[ds] = get(params, nameof(ds), 0.0)
    end
    
    for is in initial_condition_parameters(model)
        set[is] = get(params, nameof(is), 0.0)
    end

    # Domain size
    if !is_fixed_size(model)
        L = domain_size(model)
        set[L] = get(params, nameof(L), 1.0)
    end

    set
end

parameter_set(params::ParameterSet) = params


# Constructors for different spatial discretisations
#___________________________________________________________________________________________________________________________________________________________________________________

"""
Construct a SplitODEProblem to solve a reaction diffusion system with reflective boundaries.

Returns the SplitODEProblem with solutions in the frequency (DCT-1) domain and a FFTW plan to transform solutions back to the spatial domain.
"""
function PseudoSpectralProblem(model, num_verts; kwargs...)
    L = domain_size(model)
    S = species(model)
    R = reaction_rates(model)
    D = diffusion_rates(model)/L^2
    B = -L * boundary_flux(model)
    I = initial_conditions(model)
    PseudoSpectralProblem(S, R, D, B, I, num_verts; kwargs...)
end
                              
function mol_problem(model, num_verts; noise=1e-4, kwargs...)
    @parameters t x
    Dt = Differential(t)
    Dx = Differential(x)
    Dxx = Differential(x)^2
    L = domain_size(model)
    S= species(model)
    R = reaction_rates(model)
    D = diffusion_rates(model)/L^2
    B = -L * boundary_flux(model) ./ diffusion_rates(model)'
    I = initial_conditions(model)
    ps = parameters(model)

    names = model |> species .|> getname
    S_f = [only(@variables $n(..)) for n in names]
    S_tx = [u(t,x) for u in S_f]
    R_tx = [substitute(r, zip_dict(S,S_tx)) for r in R]

    eqs = [Dt(u) ~ r + d * Dxx(u) for (u,r,d) in zip(S_tx, R_tx, D)]
    ics = [u(0,x) ~ ic for (u,ic) in zip(S_f,I)]
    bcs0 = [Dx(u(t,0)) ~ bc for (u,bc) in zip(S_f,B[1,:])]
    bcs1 = [Dx(u(t,1)) ~ bc for (u,bc) in zip(S_f,B[2,:])]
    bcs = [ics; bcs0; bcs1]
   
    domains = [t ∈ (0.0, 5e5), x ∈ (0.0, 1.0)]
    p0 = zip_dict(ps, zeros(length(ps)))
    @named pdesys = PDESystem(eqs, bcs, domains, [t, x], S_tx, ps; defaults=p0)
    discretize(pdesys, MOLFiniteDifference([x => 1/(num_verts-1)], t))
end

# Macro functions
# ________________________________________________________________________________________________________

function parse_body(body, source)
    Base.remove_linenums!(body)
    parameters = ExprValues[]
    species = ExprValues[]
    pairs = Pair{ExprValues, ExprValues}[]

    for b in body.args
        r, s = b.args
        # Handle interpolation of variables
        r = parse_expr!(parameters,r)
        s = esc_dollars!(s)
        push!(pairs, s => r)
        push!(species, s)
    end

    forbidden_symbol_check(species)
    forbidden_symbol_check(parameters)
    species, parameters, pairs
end

function parse_expr!(parameters, x)
    esc_dollars!(x)
    find_parameters_in_rate!(parameters, x)
    x
end

dict_expr(pairs) = :(SpeciesValues($([:($k => $v) for (k, v) in pairs]...)))

end