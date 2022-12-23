using Documenter
using QuantumDynamics

makedocs(
    modules=[QuantumDynamics],
    sitename="QuantumDynamics",
    pages=[
        "Introduction" => "index.md",
        "Tutorial" => [
            "Path Integrals" => "./tutorial/BasicPI.md"
            "Dynamics under External Fields" => "./tutorial/ExternalFieldDynamics.md"
            "Bloch-Redfield Master Equation" => "./tutorial/Bloch-Redfield.md"
        ],
        "Documentation" => [
            "Bare System Propagation" => "./documentation/Bare.md",
            "Spectral Densities" => "./documentation/SpectralDensities.md",
            "Eta Coefficients" => "./documentation/EtaCoefficients.md",
            "QuAPI" => "./documentation/QuAPI.md",
            "Blip Decomposition" => "./documentation/Blip.md",
        ]
    ]
)
deploydocs(
    repo="github.com/amartyabose/QuantumDynamics.git"
)