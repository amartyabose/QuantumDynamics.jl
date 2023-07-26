using QuantumDynamics
using Test
using Plots


@testset "QuAPI and TEMPO" begin
    H = Utilities.create_tls_hamiltonian(; ϵ=0.0, Δ=1.0)
    ρ0 = [1.0+0.0im 0.0; 0.0 0.0]

    ξ = 0.09
    ωc = 2.5
    
    dt = 0.25
    ntimes = 100
    
    β = 0.1
    
    Jw = SpectralDensities.ExponentialCutoff(; ξ, ωc, n=1)
    
    ω, c = SpectralDensities.discretize(Jw, 100)
    hb = Solvents.HarmonicBath(β, ω, c, [1.0, -1.0], 1000)


    barefbU = Propagators.calculate_bare_propagators(; Hamiltonian=H, dt=dt, ntimes=ntimes)
    
    t, ρs = QuAPI.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, kmax=3)

    t2, ρs2 = TEMPO.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, kmax=3)

    @test maximum(abs.(real.(ρs[:,1,1] - ρs2[:,1,1]))) < 1e-3
    
end

