push!(LOAD_PATH,"../src/")

using Documenter, ReactionDiffusion

pages = Any["Home" => "index.md",
            "Tutorial" => Any["tutorial/installation.md",
            "tutorial/model.md","tutorial/params.md","tutorial/screen.md","tutorial/simulate.md","tutorial/save.md",],
            "Examples" => Any["examples/cima.md","examples/gm.md","examples/schnakenburg.md"],
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


#run in Julia repl using: include("make.jl")
#if adding docstrings, these are only updated once you restart Julia