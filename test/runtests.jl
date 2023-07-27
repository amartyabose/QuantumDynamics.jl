using QuantumDynamics
using Test
using Plots


@testset "QuAPI and TEMPO" begin
    
    # Model from Gao, et al. 
    
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

    @test maximum(abs.(real.(ρs[:,1,1] - ρs2[:,1,1]))) < 0.1
    
end

@testset "FMO HEOM" begin
    
    # FMO complex from Bose JCP paper

    num_modes = 2
    Lmax = 3

    thz2au = 0.0001519828500716
    invcm2au = 4.55633e-6
    au2fs = 0.02418884254
    
    β = 1 / (77 * 3.16683e-6)

    H = Matrix{ComplexF64}([
        12410 -87.7 5.5 -5.9 6.7 -13.7 -9.9
        -87.7 12530 30.8 8.2 0.7 11.8 4.3
        5.5 30.8 12210 -53.5 -2.2 -9.6 6.0
        -5.9 8.2 -53.5 12320 -70.7 -17.0 -63.3
        6.7 0.7 -2.2 -70.7 12480 81.1 -1.3
        -13.7 11.8 -9.6 -17.0 81.1 12630 39.7
        -9.9 4.3 6.0 -63.3 -1.3 39.7 12440
    ]) * invcm2au

    nsteps = 500
    dt = 1000 / au2fs / nsteps
    ρ0 = Matrix{ComplexF64}(zeros(7, 7))
    ρ0[1, 1] = 1

    λs = repeat([35.0], 7) * invcm2au
    γs = 1 ./ (repeat([50.0], 7) ./ au2fs)
    Jw = Vector{SpectralDensities.DrudeLorentz}()
    sys_ops = Vector{Matrix{ComplexF64}}()
    for (j, (λ, γ)) in enumerate(zip(λs, γs))
        push!(Jw, SpectralDensities.DrudeLorentz(; λ, γ, Δs=1.0))
        op = zeros(7, 7)
        op[j, j] = 1.0
        push!(sys_ops, op)
    end
    
    threshold=1e-5

    @time t, ρ = HEOM.propagate(; Hamiltonian=H, ρ0=ρ0, Jw, β, ntimes=nsteps, dt, sys_ops, num_modes, Lmax, extraargs=Utilities.DiffEqArgs(; reltol=1e-6, abstol=1e-6))

end

