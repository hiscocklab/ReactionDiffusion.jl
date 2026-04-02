module Util
using Pipe: @pipe
using Catalyst: @species, @parameters, isspecies, ExprValues
using Symbolics: get_variables, Num
using Base.Threads: @threads

## Generic Dictionaries
"Convenience function to construct a dict using (k=v, ...) syntax"
dict(;kwargs...) = Dict(kwargs)

"Return a vector of dictionaries containing the cartesian product of the given values."
product(;kwargs...) = vec([Dict(zip(keys(kwargs), vals)) for vals in Iterators.product(values(kwargs)...)])
product(dict::Dict) = product(;dict...)

"Replace each element of keys with either the corresponding value in dict or default."
function subst(keys, dict, default)
    v = Vector(undef,0)
    for k in keys
        val = get(dict, k, default)
        push!(v, val)
    end
    v
end

"Zip a keys vector and a values vector into a dictionary."
zip_dict(keys, values) = Dict(zip(keys,values))
"Unzip a dictionary into a keys vector and a values vector."
unzip_dict(dict) = (collect(keys(dict)), collect(values(dict)))



## Symbolics 
"Sort parameters by name."
sort_variables(p) = sort(p, by=_nameof)
_nameof(v) = isspecies(v) ? nameof(v.f) : nameof(v)


"Extract variables from a (possibly nested) collection of expressions and sort them by name."
collect_variables(exprs...) = collect_variables(exprs) # Combine multiple arguments.
collect_variables(exprs::Union{Tuple,Vector}) = @pipe exprs .|> collect_variables |> splat(union) |> sort_variables
collect_variables(expr) = get_variables(expr) # Call recursively until we get down to a single expression.


## Subscripted symbols
"Map integers to subscript characters."
sub(i) = join(Char(0x2080 + d) for d in reverse!(digits(i)))

"Subscript a symbol with `i...` separated by `_`."
subscript(X, i...) = Symbol(X, join(sub.(i), "_"))

"Build an array of subscripted symbols."
function subscripts(name, dim)
    dim = tuple(dim...) # ensure tuple
    ixs = Iterators.product(range.(1,dim)...)
    [subscript(name, r...) for r in ixs]
end

## Subscripted Catalyst variable arrays.
"Define a set of subscripted species and return them as a vector."
function defspecies(name, t, n)
    names = subscripts(name, n)
    [only(@species $n(t)) for n in names]
end
"Define a set of subscripted species and return them as a vector."
function defparams(name, dim)
    names = subscripts(name, dim)
    [only(@parameters $n) for n in names]
end
"Define a subscripted parameter."
function defparam(name, i...)
    name = subscript(name, i...)
    only(@parameters $name)
end


## Threads
function tmap(f, T, a)
    b = similar(a, T)
    @threads for i in eachindex(a)
        b[i] = f(a[i])
    end
    b
end

tfilter(f, a) = a[tmap(f,Bool,a)]


## Misc
issingle(x) = !(x isa AbstractVector)
isnonzero(x) = !(ismissing(x) || iszero(x))
ensure_vector(v::AbstractVector) = v
ensure_vector(x) = [x]
ensure_function(f::Function) = f
ensure_function(x) = _ -> x

decrange(start,stop) = logrange(10.0^start,10.0^stop,stop-start+1)


# Arrays
"""
Modified version of `stack` which returns stacks of zero depth when `iter` is empty instead of throwing an error.
`size` should be the size of each item in `iter`.
"""
function safe_stack(iter::Union{AbstractVector{T}, Base.Generator{<:AbstractVector{T}, S}}, size; dims=2) where T where S
    if isempty(iter)
        size = [s for s in size]
        if isempty(Base.size(size))
            size = reshape(size,1)
        end
        insert!(size, dims, 0)
        Array{T}(undef, size...)
    else
        stack(iter; dims=dims)
    end
end

end
