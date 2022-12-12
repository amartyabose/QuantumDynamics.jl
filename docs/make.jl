using Documenter
using QuantumDynamics

makedocs(
    modules = [QuantumDynamics],
    sitename = "QuantumDynamics",
    pages = [
        "Introduction" => "index.md",
        "Tutorial" => [
            "QuAPI" => "./tutorial/QuAPI.md"
        ],
        "Documentation" => [
            "Spectral Densities" => "./documentation/SpectralDensities.md",
        ]
    ]
)
deploydocs(
    repo = "github.com/amartyabose/QuantumDynamics.git"
)