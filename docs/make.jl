using Documenter
using QuantumDynamics

makedocs(
    modules=[QuantumDynamics],
    sitename="QuantumDynamics.jl",
    pages=[
        "Introduction" => "index.md",
        "Tutorial" => [
            # "Empirical Approaches" => "./tutorial/EmpiricalApproaches.md"
            # "Path Integrals" => "./tutorial/BasicPI.md"
            "Bloch-Redfield Master Equation" => "./tutorial/Bloch-Redfield.md"
            # "Quantum-Classical Path Integral" => "./tutorial/QCPI.md"
            # "Dynamics under External Fields" => "./tutorial/ExternalFieldDynamics.md"
            # "Hierarchical Equations of Motion" => "./tutorial/HEOM.md"
        ],
        "Documentation" => [
            "Bare System Propagation" => "./documentation/Bare.md",
            "Spectral Densities" => "./documentation/SpectralDensities.md",
            "Bloch-Redfield Master Equation" => "./documentation/BlochRedfield.md",
            "Incoherent Forster Theory" => "./documentation/Forster.md",
            "Eta Coefficients" => "./documentation/EtaCoefficients.md",
            "Quasi-Adiabatic Propagator Path Integral" => "./documentation/QuAPI.md",
            "Blip Decomposition" => "./documentation/Blip.md",
            "Quantum-Classical Path Integral" => "./documentation/QCPI.md",
            "Time-Evolving Matrix Product Operators" => "./documentation/TEMPO.md",
            "Pairwise-Connected Tensor Network Path Integral" => "./documentation/PCTNPI.md",
            "Transfer Tensor Method" => "./documentation/TTM.md",
            "Hierarchical Equations of Motion" => "./documentation/HEOM.md",
            "Utilities" => "./documentation/Utilities.md",
        ],
    ]
)
deploydocs(
    repo="github.com/amartyabose/QuantumDynamics.jl.git"
)