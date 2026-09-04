using ReactionDiffusion, WGLMakie
WGLMakie.activate!(resize_to=:body) 

reaction = @reaction_network begin
    # complex formation
    (k₊, k₋),               GDF5 + NOG <--> COMPLEX 
    # degradation
    δ₁,                     GDF5 --> ∅
    δ₂,                     NOG --> ∅
    δ₃,                     pSMAD --> ∅
    # transcriptional feedbacks (here: repressive hill functions)
    hillr(pSMAD,μ₁,K₁,n₁),  ∅ --> GDF5
    hillr(pSMAD,μ₂,K₂,n₂),  ∅ --> NOG
    # signalling
    μ₃*GDF5,                ∅ --> pSMAD
end

diffusion = @diffusion_system L begin
    D₁,     GDF5
    D₂,      NOG
    D₃,  COMPLEX
end

model = Model(reaction, diffusion; name="BMP Signalling Dynamics")

params = dict(
    μ₁ = 1.0,
    μ₂ = 10.0,
    k₊ = 40.0,
    k₋ = 40.0,
    μ₃ = 1.0,
    δ₁ = 0.1,
    δ₂ = 6.7,
    δ₃ = 1.0,
    K₁ = 0.01,
    K₂ = 0.01,
    n₁ = 8.0,
    n₂ = 2.0,
    D₁ = 1.0,
    D₂ = 1.0,
    D₃ = 30.0,
    L  = 25.0
)


interactive_plot(model, params; dt=0.02, num_verts=64, tspan=700, normalise=true)
