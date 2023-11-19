using BioprocessingModelLibrary
using Documenter

DocMeta.setdocmeta!(BioprocessingModelLibrary, :DocTestSetup, :(using BioprocessingModelLibrary); recursive=true)

makedocs(;
    modules=[BioprocessingModelLibrary],
    authors="Fabian Müller",
    repo="https://github.com/dfabianus/BioprocessingModelLibrary.jl/blob/{commit}{path}#{line}",
    sitename="BioprocessingModelLibrary.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://dfabianus.github.io/BioprocessingModelLibrary.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/dfabianus/BioprocessingModelLibrary.jl",
    devbranch="master",
)
