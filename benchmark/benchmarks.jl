using Pkg

tempdir = mktempdir()
Pkg.activate(tempdir)
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
Pkg.add(["BenchmarkTools", "PkgBenchmark"])

using QuantumDynamics
using BenchmarkTools

const SUITE = BenchmarkGroup()
SUITE["quapi"] = BenchmarkGroup(["path-integral", "quapi"])
SUITE["tempo"] = BenchmarkGroup(["path-integral", "tempo"])
SUITE["heom"] = BenchmarkGroup(["heom"])
SUITE["brme"] = BenchmarkGroup(["approximate", "brme"])
SUITE["lindblad"] = BenchmarkGroup(["approximate", "lindblad"])

function quapi_direct1()
    H0 = Utilities.create_tls_hamiltonian(; ϵ=0.0, Δ=2.0)
    Jw = SpectralDensities.ExponentialCutoff(; ξ=0.1, ωc=7.5)
    β = 5.0
    dt = 0.25
    ntimes = 100
    fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt, ntimes)
    ρ0 = Matrix{ComplexF64}([
        1.0 0.0
        0.0 0.0
    ])
    t, ρs = QuAPI.propagate(; fbU, Jw=[Jw], β, ρ0, dt, ntimes, kmax=10)
end

function quapi_TTM1()
    H0 = Utilities.create_tls_hamiltonian(; ϵ=0.0, Δ=2.0)
    Jw = SpectralDensities.ExponentialCutoff(; ξ=0.1, ωc=7.5)
    β = 5.0
    dt = 0.25
    ntimes = 100
    fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt, ntimes)
    ρ0 = Matrix{ComplexF64}([
        1.0 0.0
        0.0 0.0
    ])
    t, ρs = TTM.propagate(; fbU, Jw=[Jw], β, ρ0, dt, ntimes, rmax=10, path_integral_routine=QuAPI.build_augmented_propagator, extraargs=QuAPI.QuAPIArgs())
end

function tempo_direct1()
    H0 = Utilities.create_tls_hamiltonian(; ϵ=0.0, Δ=2.0)
    Jw = SpectralDensities.ExponentialCutoff(; ξ=0.1, ωc=7.5)
    β = 5.0
    dt = 0.25
    ntimes = 100
    fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt, ntimes)
    ρ0 = Matrix{ComplexF64}([
        1.0 0.0
        0.0 0.0
    ])
    t, ρs = TEMPO.propagate(; fbU, Jw=[Jw], β, ρ0, dt, ntimes, kmax=10, extraargs=TEMPO.TEMPOArgs(; cutoff=1e-15, maxdim=1000))
end

function tempo_TTM1()
    H0 = Utilities.create_tls_hamiltonian(; ϵ=0.0, Δ=2.0)
    Jw = SpectralDensities.ExponentialCutoff(; ξ=0.1, ωc=7.5)
    β = 5.0
    dt = 0.25
    ntimes = 100
    fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt, ntimes)
    ρ0 = Matrix{ComplexF64}([
        1.0 0.0
        0.0 0.0
    ])
    t, ρs = TTM.propagate(; fbU, Jw=[Jw], β, ρ0, dt, ntimes, rmax=10, path_integral_routine=TEMPO.build_augmented_propagator, extraargs=TEMPO.TEMPOArgs(; cutoff=1e-15, maxdim=1000))
end

SUITE["quapi"]["direct1"] = @benchmarkable quapi_direct1
SUITE["quapi"]["TTM1"] = @benchmarkable quapi_TTM1
SUITE["tempo"]["direct1"] = @benchmarkable tempo_direct1
SUITE["tempo"]["TTM1"] = @benchmarkable tempo_TTM1

tune!(SUITE)
run(SUITE, verbose=true)
