using Documenter, NumericalAnalysis

makedocs(modules = [NumericalAnalysis], sitename = "NumericalAnalysis.jl")

deploydocs(
    repo = "github.com/ZhouZhuofei/NumericalAnalysis.jl",
)
