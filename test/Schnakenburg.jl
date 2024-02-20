# Test scripts on Schnakenburg model (a well-known reaction-diffusion system with analytical results for its Turing stability region)

model = @reaction_network begin
    γ*a + γ*U^2*V,  ∅ --> U
    γ,              U --> ∅
    γ*b,            ∅ --> V
    γ*U^2,          V --> ∅
end 


params = model_parameters()

params.reaction["a"] = screen_values(min = 0, max = 0.6, mode = "linear", number = 4)
params.reaction["b"] = screen_values(min = 0.0, max = 3.0, mode = "linear", number = 4)
params.reaction["γ"] = [1.0]

params.diffusion["U"] = [1.0] 
params.diffusion["V"] = [50.0]


turing_params = returnTuringParams(model, params);
a = get_param(model, turing_params,"a","reaction")
b = get_param(model, turing_params,"b","reaction")

# Test whether the computed Turing parameters match the ground-truth Turing instability region
@test a == [0.0; 0.2; 0.0; 0.2]
@test b == [1.0; 1.0; 2.0; 2.0]


param1 = get_params(model, turing_params[4])

sol = simulate(model,param1)

U_final = last(sol)[:,1]
dynamicRange = maximum(U_final)/minimum(U_final) 
deviation = sign.(U_final.- 0.5*(maximum(U_final) .+ minimum(U_final)))
halfMaxSignChanges = length(findall(!iszero,diff(deviation)))

# Test whether simulated PDE is 'sensible'; we evaluate the max/min value of the final pattern, and also the number of sign changes about the half maximum (both for U)
#       note:   we give a range for both test values as we are using random initial conditions, and thus variations are to be expected
#               (even when setting seeds, it's not clear that Pkg updates to random will conserve values).
@test dynamicRange > 1.5 && dynamicRange < 4
@test halfMaxSignChanges > 3 && halfMaxSignChanges < 7