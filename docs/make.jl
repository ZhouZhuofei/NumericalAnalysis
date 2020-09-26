using Documenter, NumericalAnalysis

makedocs(
    modules = [NumericalAnalysis],
    clean = false,
    sitename = "NumericalAnalysis.jl",
    pages = Any[
        "Home" => "index.md",
        "Functions" => Any[
            "Basic function" => "m/basic.md",
            "SEq in one variable" => "m/SEq1.md",
        ],
    ],
)
