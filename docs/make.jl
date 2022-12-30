using Documenter, DocumenterCitations
using QuantumDynamics

bib = CitationBibliography("library.bib", sorting=:nyt)

makedocs(
    bib,
    modules=[QuantumDynamics],
    sitename="QuantumDynamics",
    pages=[
        "Introduction" => "index.md",
        "Tutorial" => [
            "Path Integrals" => "./tutorial/BasicPI.md"
            "Bloch-Redfield Master Equation" => "./tutorial/Bloch-Redfield.md"
            "Dynamics under External Fields" => "./tutorial/ExternalFieldDynamics.md"
        ],
        "Documentation" => [
            "Bare System Propagation" => "./documentation/Bare.md",
            "Spectral Densities" => "./documentation/SpectralDensities.md",
            "Bloch-Redfield Master Equation" => "./documentation/BlochRedfield.md",
            "Eta Coefficients" => "./documentation/EtaCoefficients.md",
            "QuAPI" => "./documentation/QuAPI.md",
            "Blip Decomposition" => "./documentation/Blip.md",
            "Utilities" => "./documentation/Utilities.md",
        ],
        "References" => "references.md"
    ]
)
deploydocs(
    repo="github.com/amartyabose/QuantumDynamics.git"
)