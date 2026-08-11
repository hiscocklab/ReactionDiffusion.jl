#push!(LOAD_PATH,"../src/")

using Documenter, ReactionDiffusion, PseudoSpectralReactionDiffusion, Catalyst

makedocs(
    sitename="ReactionDiffusion.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets=["assets/favicon.ico"]
    ),
    pages = [
       "Home" => "index.md",
       "Tutorial" => ["tutorial/installation.md", "tutorial/model.md", "tutorial/params.md", "tutorial/turing.md", "tutorial/simulate.md", "tutorial/simulate.md", "tutorial/save.md"],
       # "Examples" => Any["examples/cima.md","examples/gm.md","examples/schnakenburg.md"],
       "API" => "API/api.md",
      #  "Cheatsheet" => "cheatsheet.md"
    ],
    doctest = false,
    clean=true,
    modules=[ReactionDiffusion, PseudoSpectralReactionDiffusion, Catalyst],
    checkdocs=:none,
    checkdocs_ignored_modules = [Catalyst]
)

# deploydocs(repo = "github.com/hiscocklab/ReactionDiffusion.jl.git";
#     push_preview = true)


#run in Julia repl using: include("make.jl"), in the docs environment
#if adding docstrings, these are only updated once you restart Julia (maybe Revise.jl has taken care of this?)