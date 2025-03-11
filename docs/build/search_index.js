var documenterSearchIndex = {"docs":
[{"location":"tutorial/model.html#Specify-the-reaction-diffusion-model","page":"Specify the reaction-diffusion model","title":"Specify the reaction-diffusion model","text":"","category":"section"},{"location":"tutorial/model.html","page":"Specify the reaction-diffusion model","title":"Specify the reaction-diffusion model","text":"In the following sections, we go step-by-step through the code used in the Quick Start Guide.","category":"page"},{"location":"tutorial/model.html","page":"Specify the reaction-diffusion model","title":"Specify the reaction-diffusion model","text":"We first specify the reaction part of the reaction-diffusion system, using the intuitive syntax developed in Catalyst.jl. This allows models to be written in a very natural way which reflects the (bio)-chemical interactions involved in the system: ","category":"page"},{"location":"tutorial/model.html","page":"Specify the reaction-diffusion model","title":"Specify the reaction-diffusion model","text":"model = @reaction_network begin\n    # complex formation\n    (k₊, k₋),               GDF5 + NOG <--> COMPLEX \n    # degradation\n    δ₁,                     GDF5 --> ∅\n    δ₂,                     NOG --> ∅\n    δ₃,                     pSMAD --> ∅\n    # transcriptional feedbacks (here: repressive hill functions)\n    hillr(pSMAD,μ₁,K₁,n₁),  ∅ --> GDF5\n    hillr(pSMAD,μ₂,K₂,n₂),  ∅ --> NOG\n    # signalling\n    μ₃*GDF5,                ∅ --> pSMAD\nend  ","category":"page"},{"location":"tutorial/model.html","page":"Specify the reaction-diffusion model","title":"Specify the reaction-diffusion model","text":"Here, reaction rates are assumed to follow mass action kinetics according to the stoichiometries of the reactants. So, for example, the (k₊, k₋), GDF5 + NOG <--> COMPLEX term represents the binding/unbinding of GDF5 and NOG to form a COMPLEX. The rate of the forward reaction (i.e., the rate of complex formation) is k_+ mathrmGDF5 mathrmNOG, and the rate of the reverse reaction (i.e., the rate of complex dissociation) is k_- mathrmCOMPLEX.","category":"page"},{"location":"tutorial/model.html","page":"Specify the reaction-diffusion model","title":"Specify the reaction-diffusion model","text":"Reaction rates can also be overriden by user-specified functions, for example to denote regulatory feedbacks. In this example, the effect of pSMAD on GDF5 expression is captured by the repressive hill function hillr: fracmu_1 1 + (mathrmpSMADK_1)^n_1. ","category":"page"},{"location":"tutorial/model.html","page":"Specify the reaction-diffusion model","title":"Specify the reaction-diffusion model","text":"Arbitrary user-defined functions may also be used directly, for example the following code would reproduce the repressive hill function interaction:","category":"page"},{"location":"tutorial/model.html","page":"Specify the reaction-diffusion model","title":"Specify the reaction-diffusion model","text":"function myOwnHillrFunction(input,μ,K,n)\n    return μ/((input/K)^n + 1)\nend\n\nmodel = @reaction_network begin\n# ...\n    myOwnHillrFunction(pSMAD,μ₁,K₁,n₁),  ∅ --> GDF5\n# ...\nend","category":"page"},{"location":"tutorial/screen.html#Determine-pattern-forming-parameter-sets","page":"Determine pattern-forming parameter sets","title":"Determine pattern-forming parameter sets","text":"","category":"section"},{"location":"tutorial/screen.html","page":"Determine pattern-forming parameter sets","title":"Determine pattern-forming parameter sets","text":"Having specified the reaction-diffusion model and its parameters, we can use a single function to screen through all parameter combinations and return those that undergo a Turing instability:","category":"page"},{"location":"tutorial/screen.html","page":"Determine pattern-forming parameter sets","title":"Determine pattern-forming parameter sets","text":"turing_params = returnTuringParams(model, params);","category":"page"},{"location":"tutorial/screen.html","page":"Determine pattern-forming parameter sets","title":"Determine pattern-forming parameter sets","text":"Note: as with all Julia functions, there is a significant compilation time associated with the first time you run this function. This overhead time will not be present if you re-run the simulation e.g., for different parameters","category":"page"},{"location":"tutorial/screen.html","page":"Determine pattern-forming parameter sets","title":"Determine pattern-forming parameter sets","text":"This script returns all pattern-forming parameters found, saving the results in the turing_params variable, which has the following attributes:","category":"page"},{"location":"tutorial/screen.html","page":"Determine pattern-forming parameter sets","title":"Determine pattern-forming parameter sets","text":"steady_state_values: denotes the (homogeneous) steady state values of each of the variables, in the absence of diffusion\npattern_phase: predicts which components are in/out of phase.\nwavelength: denotes the approximate predicted wavelength of the pattern, assuming it is the one that is maximally unstable\nmax_real_eigval: records the maximally unstable real eigenvalue\nnon_oscillatory: if true, this is a stationary Turing pattern; if false, an oscillatory pattern.","category":"page"},{"location":"tutorial/screen.html","page":"Determine pattern-forming parameter sets","title":"Determine pattern-forming parameter sets","text":"turing_params also records the input values used to generate each pattern-forming parameter set, namely:","category":"page"},{"location":"tutorial/screen.html","page":"Determine pattern-forming parameter sets","title":"Determine pattern-forming parameter sets","text":"reaction_params: records the reaction parameters \ndiffusion_constants: records the diffusion constants\ninitial_conditions: records the initial condition used to compute the steady states","category":"page"},{"location":"tutorial/screen.html","page":"Determine pattern-forming parameter sets","title":"Determine pattern-forming parameter sets","text":"Individual parameters can also be accessed via helper functions:","category":"page"},{"location":"tutorial/screen.html","page":"Determine pattern-forming parameter sets","title":"Determine pattern-forming parameter sets","text":"δ₁ = get_param(model, turing_params,\"δ₁\",\"reaction\")\nD_COMPLEX = get_param(model, turing_params,\"COMPLEX\",\"diffusion\")","category":"page"},{"location":"tutorial/simulate.html#Simulate-pattern-formation-in-1D","page":"Simulate pattern formation in 1D","title":"Simulate pattern formation in 1D","text":"","category":"section"},{"location":"tutorial/simulate.html","page":"Simulate pattern formation in 1D","title":"Simulate pattern formation in 1D","text":"Having now screened through parameter sets, we can now pick a single one and simulate the corresponding PDE. The simulation is performed on a 1D domain with reflective boundary conditions.  ","category":"page"},{"location":"tutorial/simulate.html","page":"Simulate pattern formation in 1D","title":"Simulate pattern formation in 1D","text":"We first select the parameter values from turing_params; here we choose the 1000th parameter set found. ","category":"page"},{"location":"tutorial/simulate.html","page":"Simulate pattern formation in 1D","title":"Simulate pattern formation in 1D","text":"param1 = get_params(model, turing_params[1000])","category":"page"},{"location":"tutorial/simulate.html","page":"Simulate pattern formation in 1D","title":"Simulate pattern formation in 1D","text":"Optionally, we can manually change the parameters too:","category":"page"},{"location":"tutorial/simulate.html","page":"Simulate pattern formation in 1D","title":"Simulate pattern formation in 1D","text":"param1.reaction[\"n1\"] = 10","category":"page"},{"location":"tutorial/simulate.html","page":"Simulate pattern formation in 1D","title":"Simulate pattern formation in 1D","text":"The size of the 1D domain over which the simulation is performed is chosen according to the turing_params.wavelength value by default, although you can change this manually too:","category":"page"},{"location":"tutorial/simulate.html","page":"Simulate pattern formation in 1D","title":"Simulate pattern formation in 1D","text":"param1.domain_size = 100","category":"page"},{"location":"tutorial/simulate.html","page":"Simulate pattern formation in 1D","title":"Simulate pattern formation in 1D","text":"To simulate the script, you then simply run:","category":"page"},{"location":"tutorial/simulate.html","page":"Simulate pattern formation in 1D","title":"Simulate pattern formation in 1D","text":"sol = simulate(model,param1)","category":"page"},{"location":"tutorial/simulate.html","page":"Simulate pattern formation in 1D","title":"Simulate pattern formation in 1D","text":"Note: as with all Julia functions, there is a significant compilation time associated with the first time you run this function. This overhead time will not be present if you re-run the simulation e.g., for different parameters","category":"page"},{"location":"tutorial/simulate.html","page":"Simulate pattern formation in 1D","title":"Simulate pattern formation in 1D","text":"The resulting PDE dynamics are represented by the solution object sol from DifferentialEquations.jl. To visualize the results, you can use a variety of plotting packages (e.g., Makie.jl, Plots.jl). We provide simple helper functions for the Plots.jl package.","category":"page"},{"location":"tutorial/simulate.html","page":"Simulate pattern formation in 1D","title":"Simulate pattern formation in 1D","text":"To visualize the final pattern (once steady state has been reached), use: endpoint() ","category":"page"},{"location":"tutorial/simulate.html","page":"Simulate pattern formation in 1D","title":"Simulate pattern formation in 1D","text":"using Plots\nplot(endpoint(),model,sol)","category":"page"},{"location":"tutorial/simulate.html","page":"Simulate pattern formation in 1D","title":"Simulate pattern formation in 1D","text":"To visualize intermediate timepoints (e.g., for making a movie of the dynamics), use timepoint(), e.g.,:","category":"page"},{"location":"tutorial/simulate.html","page":"Simulate pattern formation in 1D","title":"Simulate pattern formation in 1D","text":"@gif for t in 0.0:0.01:1\n    plot(timepoint(),model,sol,t)\nend fps=20","category":"page"},{"location":"tutorial/simulate.html","page":"Simulate pattern formation in 1D","title":"Simulate pattern formation in 1D","text":"Here, for example, plot(timepoint(),model,sol,0.1) plots the solution at a time that is 10% through the simulation. ","category":"page"},{"location":"API/api.html#API-for-ReactionDiffusion.jl","page":"API","title":"API for ReactionDiffusion.jl","text":"","category":"section"},{"location":"API/api.html","page":"API","title":"API","text":"ReactionDiffusion aims to be an easy-to-use and computationally-efficient pipeline to simulate biologically-inspired reaction-diffusion models. It is our hope that models can be built with just a few lines of code and solved without the user having any knowledge of PDE solver methods. ","category":"page"},{"location":"API/api.html","page":"API","title":"API","text":"This is achieved by drawing from a range of numerical routines from the SciML packages, including Catalyst.jl, Symbolics.jl, ModelingToolkit.jl, and DifferentialEquations.jl. ","category":"page"},{"location":"API/api.html","page":"API","title":"API","text":"For users more familiar with PDE solvers, it is possible to specify optional arguments to e.g., control the specific solver algorithm used (see the API below). However, ReactionDiffusion does not aim to be a fully customizable PDE solving package that covers all bases; our focus is to make something that is relatively easy-to-use and performant but only for a restricted number of use cases.","category":"page"},{"location":"API/api.html","page":"API","title":"API","text":"If you require more customization or a more flexible PDE-solving framework, we highly recommend ModelingToolkit.jl and/or Catalyst.jl.","category":"page"},{"location":"API/api.html#Data-structures","page":"API","title":"Data structures","text":"","category":"section"},{"location":"API/api.html","page":"API","title":"API","text":"We introduce two structures","category":"page"},{"location":"API/api.html","page":"API","title":"API","text":"save_turing is an object that records parameter sets that undergo a Turing instability. It has the following fields:\nsteady_state_values: The computed steady state values of each variable\nreaction_params: The reaction parameters\ndiffusion_constants: The diffusion constants\ninitial_conditions: The initial conditions of each variable used to compute the steady state values\npattern_phase: The predicted phase of each of the variables in the final pattern.[1 1]would be in-phase,[1 -1]` would be out-of-phase\nwavelength: The wavelength that is maximally unstable\nmax_real_eigval: The maximum real eigenvalue associated with the Turing instability\nnon_oscillatory: If true, this parameter set represents a stationary Turing pattern. If false, the unstable mode has complex eigenvalues and thus may be oscillatory.\nmodel_parameters is an object that records parameter values for model simulations. It has the following fields:\n-reaction: reaction parameters   -diffusion: diffusion constants   -initial_condition: values of each of the variables used for homogeneous initial conditions   -initial_noise: level of (normally distributed) noise added to the homogeneous initial conditions   -domain_size: size of 1D domain   -random_seed: seed associated with random number generation; only set this when you need to reproduce exact simulations each run","category":"page"},{"location":"API/api.html#Functions","page":"API","title":"Functions","text":"","category":"section"},{"location":"API/api.html","page":"API","title":"API","text":"screen_values\nget_params\nget_param\nreturnTuringParams\nsimulate","category":"page"},{"location":"API/api.html#ReactionDiffusion.screen_values","page":"API","title":"ReactionDiffusion.screen_values","text":"screen_values(;min=0,max=1,number=10, mode=\"linear\")\n\nReturns a series of num_params parameter values that are linearly spaced between the min and max limits. The argument mode=\"log\" may be used to sample in log-space instead. \n\n\n\n\n\n","category":"function"},{"location":"API/api.html#ReactionDiffusion.get_params","page":"API","title":"ReactionDiffusion.get_params","text":"get_params(model, turing_param)\n\nFor a given model and a single pattern-forming parameter set, turing_param, this function creates a  corresponding model_parameters variable. This sets, by default, the model_parameters fields:\n\ndomain_size: chosen to be 3x the computed pattern wavelength, i.e., 3*turing_param.wavelength\ninitial_condition: chosen to be the computed steady state values, i.e.,  turing_param.steady_state_values\ninitial_noise = 0.01: the magnitude of noise (normally distributed random numbers) added to the steady state values to define the initial conditions.\n\n\n\n\n\n","category":"function"},{"location":"API/api.html#ReactionDiffusion.get_param","page":"API","title":"ReactionDiffusion.get_param","text":"get_param(model, turing_params, name, type)\n\nFor a given model and a (potentially large) number of pattern-forming parameter sets, turing_params, this function extracts the parameter values prescribed by the input name. For reaction parameters, used type=\"reaction\", for diffusion constants, use type=\"diffusion\".\n\nExample:\n\nδ₁ = get_param(model, turing_params,\"δ₁\",\"reaction\")\nD_COMPLEX = get_param(model, turing_params,\"COMPLEX\",\"diffusion\")\n\n\n\n\n\n","category":"function"},{"location":"API/api.html#ReactionDiffusion.returnTuringParams","page":"API","title":"ReactionDiffusion.returnTuringParams","text":"returnTuringParams(model, params; maxiters = 1e3,alg=Rodas5(),abstol=1e-8, reltol=1e-6, tspan=1e4,ensemblealg=EnsembleThreads(),batch_size=1e4)\n\nReturn a save_turing object of parameters that are predicted to be pattern forming.\n\nRequired inputs:\n\nmodel: specified via the @reaction_network macro\nparams: all reaction and diffusion parameters, in a model_parameters object\n\nOptional inputs:\n\nbatch_size: the number of parameter sets to consider at once. Increasing/decreasing from the default value may improve speed.  \n\nInputs carried over from DifferentialEquations.jl; see here for further details:\n\nmaxiters: maximum number of iterations to reach steady state (otherwise simulation terminates)\nalg: ODE solver algorithm\nabstol and reltol: tolerance levels of solvers\ntspan: maximum time allowed to reach steady state (otherwise simulation terminates)\nensemblealg: ensemble simulation method\n\n\n\n\n\n","category":"function"},{"location":"API/api.html#ReactionDiffusion.simulate","page":"API","title":"ReactionDiffusion.simulate","text":"simulate(model,param;alg=KenCarp4(),reltol=1e-6,abstol=1e-8, dt = 0.1, maxiters = 1e3, save_everystep = true)\n\nSimulate model for a single parameter set param.\n\nRequired inputs:\n\nmodel: specified via the @reaction_network macro\nparam: all reaction and diffusion parameters, in a model_parameters object. This must be a single parameter set only \n\nInputs carried over from DifferentialEquations.jl; see here for further details:\n\nmaxiters: maximum number of iterations to reach steady state (otherwise simulation terminates)\nalg: solver algorithm\nabstol and reltol: tolerance levels of solvers\ndt: initial value for timestep\nsave_everystep: controls whether all timepoints are saved, defaults to true\n\n\n\n\n\n","category":"function"},{"location":"examples/cima.html#CIMA-model","page":"CIMA model","title":"CIMA model","text":"","category":"section"},{"location":"examples/cima.html","page":"CIMA model","title":"CIMA model","text":"Show code and show agreement with theory","category":"page"},{"location":"tutorial/params.html#Specify-the-parameter-values","page":"Specify the parameter values","title":"Specify the parameter values","text":"","category":"section"},{"location":"tutorial/params.html","page":"Specify the parameter values","title":"Specify the parameter values","text":"We must now define parameters, starting with the reaction terms.","category":"page"},{"location":"tutorial/params.html","page":"Specify the parameter values","title":"Specify the parameter values","text":"We may specify constant parameter values:","category":"page"},{"location":"tutorial/params.html","page":"Specify the parameter values","title":"Specify the parameter values","text":"params = model_parameters()\n\n# constant parameters\nparams.reaction[\"δ₃\"] = [1.0]\nparams.reaction[\"μ₁\"] = [1.0]\nparams.reaction[\"μ₃\"] = [1.0]\nparams.reaction[\"n₁\"] = [8]\nparams.reaction[\"n₂\"] = [2]","category":"page"},{"location":"tutorial/params.html","page":"Specify the parameter values","title":"Specify the parameter values","text":"or ranges of parameter values:","category":"page"},{"location":"tutorial/params.html","page":"Specify the parameter values","title":"Specify the parameter values","text":"num_params = 5\nparams.reaction[\"δ₁\"] = screen_values(min = 0.1,max = 10, number=num_params)\nparams.reaction[\"δ₂\"] = screen_values(min = 0.1,max = 10,number=num_params)\nparams.reaction[\"μ₂\"] = screen_values(min = 0.1,max = 10, number=num_params)\nparams.reaction[\"k₊\"] = screen_values(min = 10.0,max = 100.0, number=num_params)\nparams.reaction[\"k₋\"] = screen_values(min = 10.0,max = 100.0, number=num_params)\nparams.reaction[\"K₁\"] = screen_values(min = 0.01,max = 1,number=num_params)\nparams.reaction[\"K₂\"] = screen_values(min = 0.01,max = 1, number=num_params)\nparams.reaction","category":"page"},{"location":"tutorial/params.html","page":"Specify the parameter values","title":"Specify the parameter values","text":"Here, the function screen_values returns a series of num_params parameters that are linearly spaced between the min and max limits. The argument mode=\"log\" may be used to sample in log-space instead. ","category":"page"},{"location":"tutorial/params.html","page":"Specify the parameter values","title":"Specify the parameter values","text":"Arbitrary collections of parameters may also be specified, e.g.,","category":"page"},{"location":"tutorial/params.html","page":"Specify the parameter values","title":"Specify the parameter values","text":"params.reaction[\"n₁\"] = [2; 4; 8]","category":"page"},{"location":"tutorial/params.html","page":"Specify the parameter values","title":"Specify the parameter values","text":"Next, we must specify the diffusion coefficients:","category":"page"},{"location":"tutorial/params.html","page":"Specify the parameter values","title":"Specify the parameter values","text":"params.diffusion = Dict(\n    \"NOG\"       => [1.0],\n    \"GDF5\"      => screen_values(min = 0.1,max = 30, number=10),\n    \"COMPLEX\"   => screen_values(min = 0.1,max = 30, number=10)\n)","category":"page"},{"location":"tutorial/params.html","page":"Specify the parameter values","title":"Specify the parameter values","text":"Note, any variables that are not assigned a diffusion constant are assumed to be non-diffusing (i.e., pSMAD in this example).","category":"page"},{"location":"tutorial/params.html","page":"Specify the parameter values","title":"Specify the parameter values","text":"Optionally, the initial conditions used to define the steady state(s) of the system may be explicitly set:","category":"page"},{"location":"tutorial/params.html","page":"Specify the parameter values","title":"Specify the parameter values","text":"params.initial_condition[\"NOG\"] = [0.1;1.0,10.0]","category":"page"},{"location":"tutorial/params.html","page":"Specify the parameter values","title":"Specify the parameter values","text":"This may be particularly helpful when there are more than one steady states.","category":"page"},{"location":"tutorial/save.html#Save-results-for-future-analysis","page":"Save results for future analysis","title":"Save results for future analysis","text":"","category":"section"},{"location":"tutorial/save.html","page":"Save results for future analysis","title":"Save results for future analysis","text":"Sometimes you may wish to save some results for later analysis (for example, if you perform a particularly large parameter screen). You can use JLD2 to save these files, e.g.,","category":"page"},{"location":"tutorial/save.html","page":"Save results for future analysis","title":"Save results for future analysis","text":"##Save solution for a later date\n@save \"test.jld2\" model params turing_params sol","category":"page"},{"location":"tutorial/save.html","page":"Save results for future analysis","title":"Save results for future analysis","text":"Then, to load the files:","category":"page"},{"location":"tutorial/save.html","page":"Save results for future analysis","title":"Save results for future analysis","text":"@load \"test.jld2\"","category":"page"},{"location":"examples/schnakenburg.html#Schnakenburg-model","page":"Schnakenburg model","title":"Schnakenburg model","text":"","category":"section"},{"location":"examples/schnakenburg.html","page":"Schnakenburg model","title":"Schnakenburg model","text":"Show code and show agreement with theory","category":"page"},{"location":"tutorial/installation.html#Installation-and-setup","page":"Installation and setup","title":"Installation and setup","text":"","category":"section"},{"location":"tutorial/installation.html","page":"Installation and setup","title":"Installation and setup","text":"ReactionDiffusion can be installed using the Julia package manager:","category":"page"},{"location":"tutorial/installation.html","page":"Installation and setup","title":"Installation and setup","text":"using Pkg\nPkg.add(\"ReactionDiffusion\")","category":"page"},{"location":"tutorial/installation.html","page":"Installation and setup","title":"Installation and setup","text":"and then loaded ready for use:","category":"page"},{"location":"tutorial/installation.html","page":"Installation and setup","title":"Installation and setup","text":"using ReactionDiffusion","category":"page"},{"location":"tutorial/installation.html","page":"Installation and setup","title":"Installation and setup","text":"As with all Julia packages, upon first loading there will significant precompilation time, but this step is only required once per installation.","category":"page"},{"location":"tutorial/installation.html","page":"Installation and setup","title":"Installation and setup","text":"You may also wish to install a plotting package to visualize your results; we provide plotting scripts which use Plots.jl.","category":"page"},{"location":"tutorial/installation.html","page":"Installation and setup","title":"Installation and setup","text":"Pkg.add(\"Plots\")\nusing Plots","category":"page"},{"location":"tutorial/installation.html","page":"Installation and setup","title":"Installation and setup","text":"If you don't already use Julia, this must first be installed, along with an environment to edit and run your scripts. We recommend following the guides from SciML to set this up. ","category":"page"},{"location":"tutorial/installation.html","page":"Installation and setup","title":"Installation and setup","text":"We also recommend using environments to manage your package dependencies; see here for a more detailed discussion. ","category":"page"},{"location":"tutorial/installation.html","page":"Installation and setup","title":"Installation and setup","text":"To take advantage of automatic multithreading in your simulations, ensure that the number of Julia threads is set according to the CPU on your machine (we recommend to set Threads.nthreads() to be equal to the number of cores available). In VSCode, this can be achieved by the following:","category":"page"},{"location":"tutorial/installation.html","page":"Installation and setup","title":"Installation and setup","text":"go to Settings Ctrl-,\nsearch for \"threads\"\nclick edit in settings.json\nset \"julia.NumThreads\": 8 (choose the number of cores on your machine)","category":"page"},{"location":"index.html#ReactionDiffusion.jl-for-modelling-pattern-formation-in-biological-systems","page":"Home","title":"ReactionDiffusion.jl for modelling pattern formation in biological systems","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"Reaction-diffusion dynamics are present across many areas of the physical and natural world, and allow complex spatiotemporal patterns to self-organize de novo. ReactionDiffusion.jl aims to be an easy-to-use and performant pipeline to simulate reaction-diffusion PDEs of arbitrary complexity, with a focus on pattern formation in biological systems. Using this package, complex, biologically-inspired reaction-diffusion models can be:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"specified using an intuitive syntax\nscreened across millions of parameter combinations to identify pattern-forming networks (i.e., those that undergo a Turing instability)\nrapidly simulated to predict spatiotemporal patterns","category":"page"},{"location":"index.html#Quick-start-guide","page":"Home","title":"Quick start guide","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"Here we show how ReactionDiffusion.jl can be used to simulate a biologically-inspired reaction-diffusion system, responsible for generating evenly spaced joints along the length of your fingers and toes (from Grall et el, 2024).","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"We begin by specifying the reaction-diffusion dynamics via the intuitive syntax developed in Catalyst.jl, which naturally mirrors biochemical feedbacks and interactions.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"using ReactionDiffusion\n\nmodel = @reaction_network begin\n    # complex formation\n    (k₊, k₋),               GDF5 + NOG <--> COMPLEX \n    # degradation\n    δ₁,                     GDF5 --> ∅\n    δ₂,                     NOG --> ∅\n    δ₃,                     pSMAD --> ∅\n    # transcriptional feedbacks (here: repressive hill functions)\n    hillr(pSMAD,μ₁,K₁,n₁),  ∅ --> GDF5\n    hillr(pSMAD,μ₂,K₂,n₂),  ∅ --> NOG\n    # signalling\n    μ₃*GDF5,                ∅ --> pSMAD\nend  ","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"We can then specify values for each parameter:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"params = model_parameters()\n\n# constant parameters\nparams.reaction[\"δ₃\"] = [1.0]\nparams.reaction[\"μ₁\"] = [1.0]\nparams.reaction[\"μ₃\"] = [1.0]\nparams.reaction[\"n₁\"] = [8]\nparams.reaction[\"n₂\"] = [2]\n\n# varying parameters\nnum_params = 5\nparams.reaction[\"δ₁\"] = screen_values(min = 0.1,max = 10, number=num_params)\nparams.reaction[\"δ₂\"] = screen_values(min = 0.1,max = 10,number=num_params)\nparams.reaction[\"μ₂\"] = screen_values(min = 0.1,max = 10, number=num_params)\nparams.reaction[\"k₊\"] = screen_values(min = 10.0,max = 100.0, number=num_params)\nparams.reaction[\"k₋\"] = screen_values(min = 10.0,max = 100.0, number=num_params)\nparams.reaction[\"K₁\"] = screen_values(min = 0.01,max = 1,number=num_params)\nparams.reaction[\"K₂\"] = screen_values(min = 0.01,max = 1, number=num_params)\nparams.reaction\n","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"We must then also specify the diffusion coefficients for each variable:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"params.diffusion = Dict(\n    \"NOG\"       => [1.0],\n    \"GDF5\"      => screen_values(min = 0.1,max = 30, number=10),\n    \"COMPLEX\"   => screen_values(min = 0.1,max = 30, number=10)\n)\nparams.diffusion","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Then, with a single line of code, we can perform a Turing instability analysis across all combinations of parameters:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"turing_params = returnTuringParams(model, params);","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"This returns all parameter combinations that can break symmetry from a homogeneous initial condition. We take advantage of the highly performant numerical solvers in DifferentialEquations.jl to be able to simulate millions of parameter sets per minute on a standard laptop. ","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"We may then take a single parameter set and simulate its spatiotemporal dynamics directly, using Plots.jl to visualize the resulting pattern:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"param1 = get_params(model, turing_params[1000])\nsol = simulate(model,param1)\n\nusing Plots\nplot(endpoint(),model,sol)","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"We may also view the full spatiotemporal dynamics:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"@gif for t in 0.0:0.01:1\n    plot(timepoint(),model,sol,t)\nend fps=20","category":"page"},{"location":"index.html#Support,-citation-and-future-developments","page":"Home","title":"Support, citation and future developments","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"If you find ReactionDiffusion.jl helpful in your research, teaching, or other activities, please star the repository and consider citing this paper:","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"@article{TBD,\n doi = {TBD},\n author = {Muzatko, Daniel AND Daga, Bijoy AND Hiscock, Tom W.},\n journal = {biorXiv},\n publisher = {TBD},\n title = {TBD},\n year = {TBD},\n month = {TBD},\n volume = {TBD},\n url = {TBD},\n pages = {TBD},\n number = {TBD},\n}","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"We are a small team of academic researchers from the Hiscock Lab, who build mathematical models of developing embryos and tissues. We have found these scripts helpful in our own research, and make them available in case you find them helpful in your research too. We hope to extend the functionality of ReactionDiffusion.jl as our future projects, funding and time allows.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"This work is supported by ERC grant SELFORG-101161207, and UK Research and Innovation (Biotechnology and Biological Sciences Research Council, grant number BB/W003619/1) ","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"Funded by the European Union. Views and opinions expressed are however those of the author(s) only and do not necessarily reflect those of the European Union or the European Research Council Executive Agency. Neither the European Union nor the granting authority can be held responsible for them","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"(Image: ERC_logo)","category":"page"},{"location":"examples/gm.html#Gierer-Meinhardt-model","page":"Gierer-Meinhardt model","title":"Gierer-Meinhardt model","text":"","category":"section"},{"location":"examples/gm.html","page":"Gierer-Meinhardt model","title":"Gierer-Meinhardt model","text":"Show code and show agreement with theory","category":"page"}]
}
