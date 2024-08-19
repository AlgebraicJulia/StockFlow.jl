using Documenter
using Literate

const literate_dir = joinpath(@__DIR__, "literate")
const generated_dir = joinpath(@__DIR__, "src", "generated")

@info "Loading StockFlow.jl"
using StockFlow

const no_literate = "--no-literate" in ARGS
if !no_literate
  @info "Building Literate.jl docs"

  # XXX: Work around old LaTeX distribution in GitHub CI.
  if haskey(ENV, "GITHUB_ACTIONS")
    import TikzPictures
    TikzPictures.standaloneWorkaround(true)
  end

  # Set Literate.jl config if not being compiled on recognized service.
  config = Dict{String,String}()
  if !(haskey(ENV, "GITHUB_ACTIONS") || haskey(ENV, "GITLAB_CI"))
    config["nbviewer_root_url"] = "https://nbviewer.jupyter.org/github/AlgebraicJulia/StockFlow.jl/blob/gh-pages/dev"
    config["repo_root_url"] = "https://github.com/AlgebraicJulia/StockFlow.jl/blob/main/docs"
  end

  for (root, dirs, files) in walkdir(literate_dir)
    out_dir = joinpath(generated_dir, relpath(root, literate_dir))
    for file in files
      if last(splitext(file)) == ".jl"
        Literate.markdown(joinpath(root, file), out_dir;
                          config=config, documenter=true, credit=false)
        Literate.notebook(joinpath(root, file), out_dir;
                          execute=true, documenter=true, credit=false)
      end
    end
  end
end

@info "Building Documenter.jl docs"
makedocs(
  modules     = [StockFlow],
  format      = Documenter.HTML(
                                assets = ["assets/analytics.js"],
                               ),
  sitename    = "StockFlow.jl",
  doctest     = false,
  checkdocs   = :none,
  pages       = Any[
    "StockFlow.jl" => "index.md",
    "Domain Specific Languages" => "DSLs.md",
    "Serialization" => "json.md",
    "Example" => Any[
      "generated/Covid19_composition_model_in_paper.md",
      "generated/full_fledged_schema_examples_new/CasualLoopDiagrams/convert_from_SEIR_stockFlowDiagram.md",
      "generated/full_fledged_schema_examples_new/CasualLoopDiagrams/Causal_Loop_Operations.md",
      "generated/full_fledged_schema_examples_new/composition/composed_open_population_SIRV_model.md",
      "generated/full_fledged_schema_examples_new/composition/COVID_full_model.md",
      "generated/full_fledged_schema_examples_new/composition/curable_sexually_transmitted_diseases_model.md",
      "generated/full_fledged_schema_examples_new/composition/diabetes_model.md",
      "generated/full_fledged_schema_examples_new/composition/SEIR_full_model_measles_chickenpox.md",
      "generated/full_fledged_schema_examples_new/stratification/diabetes_diagnose.md",
      "generated/full_fledged_schema_examples_new/stratification/sir_linear_stratification.md",
      "generated/full_fledged_schema_examples_new/stratification/sir_standard_stratification.md",
    ],
    "Exercises" => Any[
      "practices/SEIRVD/SEIRVD_model_hard.md",
      "practices/SIRV/SIRV_composition_model_simple.md"
    ]
  ]
)

@info "Deploying do cs"
deploydocs(
  target = "build",
  repo   = "github.com/AlgebraicJulia/StockFlow.jl.git",
  branch = "gh-pages"
)
