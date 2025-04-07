# Schnakenberg model

The Schnakenberg model is an activator-inhibitor reaction-diffusion model capable of generating periodic stable patterns. It is relatively simple, yet chemically sensible, and therefore has been widely used and studied.

The analytical solution of this system is known. Hence, we use it here as a case study to show agreement between our numerical pipeline and the analytical solutions.

The Schnakenberg system can be expressed as a system of two PDEs:

$\frac{\partial U}{\partial t} = D_U\nabla^2U+\gamma\left(a-U+U^2V\right)$

$\frac{\partial V}{\partial t} =D_V\nabla^2 V+\gamma \left(b-U^2V\right)$
## Modelling Schnakenberg model
Implementation of the Schnakenberg model into our pipeline is straightforward: 
 
 $U$ is produced at a rate of $γa + γU^2V$ and degraded at a rate of $γ$;
 
  $V$ is produced at a rate $γb$ and degraded at a rate $γU^2$:



```@example schnakenberg
using ReactionDiffusion

model = @reaction_network begin
    γ*a + γ*U^2*V,  ∅ --> U
    γ,              U --> ∅
    γ*b,            ∅ --> V
    γ*U^2,          V --> ∅
end
```


## Specifying model parameters

Next, we define the parameter ranges. We define  
 $0\leq a\leq 0.5$,  
$0\leq b\leq 3$,  
$γ=1$,    
$D_u=1$,  
$D_v=50.$

```@example schnakenberg; output=false
params = model_parameters()
N_screen = 100
params.reaction["a"] = screen_values(min = 0, max = 0.5, mode = "linear", number = N_screen)
params.reaction["b"] = screen_values(min = 0, max = 3, mode = "linear", number = N_screen)
params.reaction["γ"] = [1]

params.diffusion["U"] = [1.0] 
params.diffusion["V"] = [50.0]
```
## Finding pattern-forming parameter sets
To find all the pattern-forming parameter sets we use built-in function `returnTuringParams`:
```@example schnakenberg
@time turing_params = returnTuringParams(model, params);
```
## Analytical solution
We want to compare our numerical solution to the general analytical solution.

Considering the Schnakenberg system near the steady state

$u_0=a+b,~~~~~v_0=\frac{b}{u_0^2}=\frac{u_0-a}{u_0^2},$

we can perform a standard analysis of this linearized system, testing for Turing instability. For details, see for example [C. Beentjes, 2014](https://cbeentjes.github.io/files/Ramblings/PatternFormationSchnakenberg.pdf).

It has been shown that the lower bound satisfying the Turing instability conditions may be expressed as:

$a_{\ell}=-b+\frac{1}{3^{1/3}\left(-9b+\sqrt{3+81b^2}\right)^{1/3}}-\frac{\left(-9b+\sqrt{3+81b^2}\right)^{1/3}}{3^{2/3}},$

and the upper bound as 

$a_{\mathcal{U}}=\frac{1}{3}\left[-3b-2\sqrt{d}+\frac{d}{\left(27bd+d^{3/2}+\sqrt{729b^2d^2+54bd^{5/2}}\right)^{1/3}}+\left(27bd+d^{3/2}+\sqrt{729b^2d^2+54bd^{5/2}}\right)^{1/3}\right],$

where $d:=D_U/D_V$. Notice that in our case $d=50$.

## Plotting analytical and numerical solution
Now we plot our numerical and analytical solutions for comparison.
To plot the solution, we use the 'Plots' package. To access the steady state values, we need to use `RecursiveArrayTools`
```@example schnakenberg
using Plots
using RecursiveArrayTools
```

We call the steady states and the parameters. For generating all the pattern-forming parameter sets, we use build-in function `get_param`:
```@example schnakenberg; output=false
U₀ = VectorOfArray(turing_params.steady_state_values)[1,:]
V₀ = VectorOfArray(turing_params.steady_state_values)[2,:]
γ =get_param(model, turing_params,"γ","reaction")
a = get_param(model, turing_params,"a","reaction")
b = get_param(model, turing_params,"b","reaction")
```

Let's start by plotting the analytical solution. We plot the lower bound $a_{\ell}$ and upper bound $a_{\mathcal{U}}$ which depend on $d$ and $b$, as discussed above.
```@example schnakenberg
b̃ = range(0,stop=3,length=100)
# Define the diffusion coefficient d
d = params.diffusion["V"]/params.diffusion["U"]

b̃H = (1 ./ (1.44 .* (-9*b̃ .+ (sqrt.(3 .+ 81*b̃.^2))).^(1/3))) .- b̃ .- ((-9*b̃ .+ (sqrt.(3 .+ 81*b̃.^2))).^(1/3)) ./ (2.08)

h = (27*b̃.*d .+ d.^(3/2) .+ sqrt.(729*b̃.^2 .* d.^2 .+ 54*b̃.*d.^(5/2))).^(1/3)
b̃T = (1/3).*(-3*b̃ .- 2*sqrt.(d) .+ d./h .+ h )

Plots.plot(b̃,[b̃H,b̃T], labels=["lower bound (analytical)" "upper bound (analytical)"], ylimits =(0,0.45), xlimits =(0, 3), xlabel="b", ylabel="a", linewidth = 5)
```

As a last step we plot all the parameter sets which are Turing-pattern-forming:
```@example schnakenberg
Plots.scatter!(b,a,strokewidth = 0,label = "numerical solution", color = :black, markersize=1)
title!("Analytical and numerical sol. - Schnakenberg model")
```
The numerical solutions nicely fill the whole region defined by the analytical boundaries, showing that the `ReactionDiffusion` package correctly identifies pattern-forming parameter sets.



