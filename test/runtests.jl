using QuantumDynamics
using Test
using DelimitedFiles, TOML

@testset "QuAPI and TEMPO" begin
    
    # Spin Boson Model A from Gao, et al. 
    
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

@testset "QCPI" begin
    # Spin Boson Model B from Gao et al.
    H = Utilities.create_tls_hamiltonian(; ϵ=0.0, Δ=1.0)
    ρ0 = [1.0+0.0im 0.0; 0.0 0.0]
    
    ξ = 0.09
    ωc = 2.5
    
    Jw = SpectralDensities.ExponentialCutoff(; ξ, ωc, n=1)
    
    dt = 0.25
    ntimes = 100
    
    β = 5.0

    ω, c = SpectralDensities.discretize(Jw, 100)
    hb = Solvents.HarmonicBath(β, ω, c, [1.0, -1.0], 1000)
    
    barefbU = Propagators.calculate_bare_propagators(; Hamiltonian=H, dt=dt, ntimes=ntimes)

    @time t3, ρs3 = QCPI.propagate(; Hamiltonian=H, Jw, solvent=hb, ρ0, classical_dt = dt/100, dt, ntimes, kmax=3, extraargs=QuAPI.QuAPIArgs(), path_integral_routine=QuAPI.propagate)

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

@testset "TTM FMO" begin
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
    for (j, (λ, γ)) in enumerate(zip(λs, γs))
        push!(Jw, SpectralDensities.DrudeLorentz(; λ, γ, Δs=1.0))
    end

    threshold=1e-5 
    
    barefbU = Propagators.calculate_bare_propagators(; Hamiltonian=H, dt=dt, ntimes=nsteps)

    svec = [1.0 -1.0 ; 1.0 -1.0; 1.0 -1.0; 1.0 -1.0 ; 1.0 -1.0 ; 1.0 -1.0; 1.0 -1.0]
    
    @time times, ρ = TTM.propagate(; fbU=barefbU, ρ0=ρ0, Jw=Jw, β, ntimes=nsteps, dt, svec, rmax=1, extraargs=QuAPI.QuAPIArgs(), path_integral_routine=QuAPI.build_augmented_propagator)
end

@testset verbose = true "Path Integrals" begin
    Hamiltonian = Utilities.create_tls_hamiltonian(; ϵ=0.0, Δ=2.0)
    Jw = SpectralDensities.ExponentialCutoff(; ξ=0.1, ωc=7.5)
    β = 5.0
    dt = 0.1
    kmax = 7

    @testset "η-coefficients" begin
        correct_etas = TOML.parsefile("correct_etas.out")
        ηs = EtaCoefficients.calculate_η(Jw; β, dt, kmax)
        @test ηs.η00 .≈ correct_etas["eta_00"][1] + 1im * correct_etas["eta_00"][2] atol = 1e-5
        @test ηs.ηmm .≈ correct_etas["eta_mm"][1] + 1im * correct_etas["eta_mm"][2] atol = 1e-5
        for k = 1:kmax
            @test ηs.ηmn[k] .≈ correct_etas["eta_mn"][k][1] + 1im * correct_etas["eta_mn"][k][2] atol = 1e-5
            @test ηs.η0e[k] .≈ correct_etas["eta_0e"][k][1] + 1im * correct_etas["eta_0e"][k][2] atol = 1e-5
            @test ηs.η0m[k] .≈ correct_etas["eta_0m"][k][1] + 1im * correct_etas["eta_0m"][k][2] atol = 1e-5
        end
    end

    ntimes = 7

    fbU = Propagators.calculate_bare_propagators(; Hamiltonian, dt, ntimes)

    dat = readdlm("correct_output.out"; skipstart=1)
    t_dat = dat[:, 1]
    re_dat = dat[:, 2:2:8]
    im_dat = dat[:, 3:2:9]
    complex_dat = reshape(re_dat .+ 1im * im_dat, (length(t_dat), 2, 2))

    ρ0 = [1.0+0.0im 0.0; 0.0 0.0]
    svec = [1.0 -1.0]

    @testset "QuAPI" begin
        t, ρs = QuAPI.propagate(; fbU, Jw=[Jw], β, ρ0, dt, ntimes, kmax, svec)

        @test all(ρs[:, 1, 2] .≈ conj(ρs[:, 2, 1]))
        @test all(ρs[:, 1, 1] .+ ρs[:, 2, 2] .≈ 1.0 + 0.0im)
        @test all(norm.(ρs[:, 1, 1] .- complex_dat[1:ntimes+1, 1, 1]) .< 1e-2)
        @test all(norm.(ρs[:, 1, 2] .- complex_dat[1:ntimes+1, 2, 1]) .< 1e-2)
    end

    @testset "TEMPO" begin
        t, ρs = TEMPO.propagate(; fbU, Jw=[Jw], β, ρ0, dt, ntimes, kmax, svec, extraargs=TEMPO.TEMPOArgs(; cutoff=1e-15, maxdim=1000))

        @test all(ρs[:, 1, 2] .≈ conj(ρs[:, 2, 1]))
        @test all(ρs[:, 1, 1] .+ ρs[:, 2, 2] .≈ 1.0 + 0.0im)
        @test all(norm.(ρs[:, 1, 1] .- complex_dat[1:ntimes+1, 1, 1]) .< 1e-2)
        @test all(norm.(ρs[:, 1, 2] .- complex_dat[1:ntimes+1, 2, 1]) .< 1e-2)
    end

    @testset "PCTNPI" begin
        U = PCTNPI.build_augmented_propagator(; fbU, Jw=[Jw], β, dt, ntimes, svec)
        t, ρs = Utilities.apply_propagator(; propagators=U, ρ0, ntimes, dt)

        @test all(ρs[:, 1, 2] .≈ conj(ρs[:, 2, 1]))
        @test all(ρs[:, 1, 1] .+ ρs[:, 2, 2] .≈ 1.0 + 0.0im)
        @test all(norm.(ρs[:, 1, 1] .- complex_dat[1:ntimes+1, 1, 1]) .< 1e-2)
        @test all(norm.(ρs[:, 1, 2] .- complex_dat[1:ntimes+1, 2, 1]) .< 1e-2)
    end

    @testset "Blip" begin
        U = Blip.build_augmented_propagator(; fbU, Jw=[Jw], β, dt, ntimes, svec)
        t, ρs = Utilities.apply_propagator(; propagators=U, ρ0, ntimes, dt)

        @test all(ρs[:, 1, 2] .≈ conj(ρs[:, 2, 1]))
        @test all(ρs[:, 1, 1] .+ ρs[:, 2, 2] .≈ 1.0 + 0.0im)
        @test all(norm.(ρs[:, 1, 1] .- complex_dat[1:ntimes+1, 1, 1]) .< 1e-2)
        @test all(norm.(ρs[:, 1, 2] .- complex_dat[1:ntimes+1, 2, 1]) .< 1e-2)
    end

    ntimes = 20
    @testset "QuAPI-TTM" begin
        t, ρs = TTM.propagate(; fbU, Jw=[Jw], β, ρ0, dt, ntimes, rmax=kmax, svec, path_integral_routine=QuAPI.build_augmented_propagator, extraargs=QuAPI.QuAPIArgs())

        @test all(ρs[:, 1, 2] .≈ conj(ρs[:, 2, 1]))
        @test all(ρs[:, 1, 1] .+ ρs[:, 2, 2] .≈ 1.0 + 0.0im)
        @test all(norm.(ρs[:, 1, 1] .- complex_dat[1:ntimes+1, 1, 1]) .< 1e-2)
        @test all(norm.(ρs[:, 1, 2] .- complex_dat[1:ntimes+1, 2, 1]) .< 1e-2)
    end

    ntimes = 20
    @testset "Blip-TTM" begin
        t, ρs = TTM.propagate(; fbU, Jw=[Jw], β, ρ0, dt, ntimes, rmax=kmax, svec, path_integral_routine=Blip.build_augmented_propagator, extraargs=Blip.BlipArgs())

        @test all(ρs[:, 1, 2] .≈ conj(ρs[:, 2, 1]))
        @test all(ρs[:, 1, 1] .+ ρs[:, 2, 2] .≈ 1.0 + 0.0im)
        @test all(norm.(ρs[:, 1, 1] .- complex_dat[1:ntimes+1, 1, 1]) .< 1e-2)
        @test all(norm.(ρs[:, 1, 2] .- complex_dat[1:ntimes+1, 2, 1]) .< 1e-2)
    end

    ntimes = 20
    @testset "TEMPO-TTM" begin
        t, ρs = TTM.propagate(; fbU, Jw=[Jw], β, ρ0, dt, ntimes, rmax=kmax, svec, path_integral_routine=TEMPO.build_augmented_propagator, extraargs=TEMPO.TEMPOArgs(; cutoff=1e-15, maxdim=1000))

        @test all(ρs[:, 1, 2] .≈ conj(ρs[:, 2, 1]))
        @test all(ρs[:, 1, 1] .+ ρs[:, 2, 2] .≈ 1.0 + 0.0im)
        @test all(norm.(ρs[:, 1, 1] .- complex_dat[1:ntimes+1, 1, 1]) .< 1e-2)
        @test all(norm.(ρs[:, 1, 2] .- complex_dat[1:ntimes+1, 2, 1]) .< 1e-2)
    end

    ntimes = 20
    @testset "PCTNPI-TTM" begin
        t, ρs = TTM.propagate(; fbU, Jw=[Jw], β, ρ0, dt, ntimes, rmax=kmax, svec, path_integral_routine=PCTNPI.build_augmented_propagator, extraargs=PCTNPI.PCTNPIArgs())

        @test all(ρs[:, 1, 2] .≈ conj(ρs[:, 2, 1]))
        @test all(ρs[:, 1, 1] .+ ρs[:, 2, 2] .≈ 1.0 + 0.0im)
        @test all(norm.(ρs[:, 1, 1] .- complex_dat[1:ntimes+1, 1, 1]) .< 1e-2)
        @test all(norm.(ρs[:, 1, 2] .- complex_dat[1:ntimes+1, 2, 1]) .< 1e-2)
    end
end
