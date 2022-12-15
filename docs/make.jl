using Documenter
using QuantumDynamics

makedocs(
    modules=[QuantumDynamics],
    sitename="QuantumDynamics",
    pages=[
        "Introduction" => "index.md",
        "Tutorial" => [
            "QuAPI" => "./tutorial/QuAPI.md"
        ],
        "Documentation" => [
            "Spectral Densities" => "./documentation/SpectralDensities.md",
            "Eta Coefficients" => "./documentation/EtaCoefficients.md",
            "QuAPI" => "./documentation/QuAPI.md",
            "Blip Decomposition" => "./documentation/Blip.md",
            "Bare System Propagation" => "./documentation/Blip.md",
        ]
    ]
)
deploydocs(
    repo="github.com/amartyabose/QuantumDynamics.git"
)