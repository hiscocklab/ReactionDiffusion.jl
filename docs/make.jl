push!(LOAD_PATH,"../src/")

using Documenter, ReactionDiffusion

pages = Any["Home" => "index.md",
            "Tutorial" => Any["tutorial/installation.md",
            "tutorial/model.md","tutorial/params.md","tutorial/screen.md","tutorial/simulate.md","tutorial/save.md",],
            "Examples" => Any["examples/schnakenberg.md", "examples/nodal.md"],
            "API" => "API/api.md"
            ]

makedocs(sitename="ReactionDiffusion.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        assets=["assets/favicon.ico"]
    ),
    pages=pages,
    doctest =false,
    clean=true,
    modules=[ReactionDiffusion],

)

deploydocs(repo = "github.com/hiscocklab/ReactionDiffusion.jl.git";
    push_preview = true)


#run in Julia repl using: include("make.jl"), in the docs environment
#if adding docstrings, these are only updated once you restart Julia (maybe Revise.jl has taken care of this?)