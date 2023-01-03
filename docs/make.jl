using Documenter
using QuantumDynamics

# bib = CitationBibliography("library.bib", sorting=:nyt)

makedocs(
    # bib,
    modules=[QuantumDynamics],
    sitename="QuantumDynamics",
    pages=[
        "Introduction" => "index.md",
        "Tutorial" => [
            "Empirical Approaches" => "./tutorial/EmpiricalApproaches.md"
            "Path Integrals" => "./tutorial/BasicPI.md"
            "Quantum-Classical Path Integral" => "./tutorial/QCPI.md"
            "Bloch-Redfield Master Equation" => "./tutorial/Bloch-Redfield.md"
            "Dynamics under External Fields" => "./tutorial/ExternalFieldDynamics.md"
            "Hierarchy Equations of Motion" => "./tutorial/HEOM.md"
        ],
        "Documentation" => [
            "Bare System Propagation" => "./documentation/Bare.md",
            "Spectral Densities" => "./documentation/SpectralDensities.md",
            "Bloch-Redfield Master Equation" => "./documentation/BlochRedfield.md",
            "Eta Coefficients" => "./documentation/EtaCoefficients.md",
            "Quasi-Adiabatic Propagator Path Integral" => "./documentation/QuAPI.md",
            "Blip Decomposition" => "./documentation/Blip.md",
            "Quantum-Classical Path Integral" => "./documentation/QCPI.md",
            "Tensor Network Path Integral" => "./documentation/TNPI.md",
            "Hierarchy Equations of Motion" => "./documentation/HEOM.md",
            "Utilities" => "./documentation/Utilities.md",
        ],
        # "References" => "references.md"
    ]
)
deploydocs(
    repo="github.com/amartyabose/QuantumDynamics.git"
)