@testmodule pisetup begin
    using QuantumDynamics, DelimitedFiles

    Hamiltonian = Utilities.create_tls_hamiltonian(; ϵ=0.0, Δ=2.0)
    Jw = SpectralDensities.ExponentialCutoff(; ξ=0.1, ωc=7.5)
    β = 5.0
    dt = 0.1
    kmax = 7
    
    ntimes = 7
    ntimes_long = 20
    ρ0 = [1.0+0.0im 0.0; 0.0 0.0]
    svec = [1.0 -1.0]

    fbU = Propagators.calculate_bare_propagators(; Hamiltonian, dt, ntimes=ntimes+1)
    U = Propagators.calculate_bare_propagators(; Hamiltonian, dt, ntimes=ntimes+1, forward_backward=false)

    dat = readdlm("correct_output.out"; skipstart=1)
    t_dat = dat[:, 1]
    re_dat = dat[:, 2:2:8]
    im_dat = dat[:, 3:2:9]
    complex_dat = reshape(re_dat .+ 1im * im_dat, (length(t_dat), 2, 2))
end

@testitem "η-coefficients" setup=[pisetup] begin
    using TOML
    correct_etas = TOML.parsefile("correct_etas.out")
    ηs = EtaCoefficients.calculate_η(pisetup.Jw; pisetup.β, pisetup.dt, pisetup.kmax)
    @test ηs.η00 .≈ correct_etas["eta_00"][1] + 1im * correct_etas["eta_00"][2] atol = 1e-5
    @test ηs.ηmm .≈ correct_etas["eta_mm"][1] + 1im * correct_etas["eta_mm"][2] atol = 1e-5
    for k = 1:pisetup.kmax
        @test ηs.ηmn[k] .≈ correct_etas["eta_mn"][k][1] + 1im * correct_etas["eta_mn"][k][2] atol = 1e-5
        @test ηs.η0e[k] .≈ correct_etas["eta_0e"][k][1] + 1im * correct_etas["eta_0e"][k][2] atol = 1e-5
        @test ηs.η0m[k] .≈ correct_etas["eta_0m"][k][1] + 1im * correct_etas["eta_0m"][k][2] atol = 1e-5
    end
end

@testitem "QuAPI" setup=[pisetup] begin
    Us = QuAPI.build_augmented_propagator(; pisetup.fbU, Jw=[pisetup.Jw], pisetup.β, pisetup.dt, pisetup.ntimes, pisetup.kmax, pisetup.svec)
    t_prop, ρs_prop = Utilities.apply_propagator(; propagators=Us, pisetup.ρ0, pisetup.dt, pisetup.ntimes)
    t, ρs = QuAPI.propagate(; pisetup.fbU, Jw=[pisetup.Jw], pisetup.β, pisetup.ρ0, pisetup.dt, pisetup.ntimes, pisetup.kmax, pisetup.svec)

    @test all(ρs .≈ ρs_prop)
    @test all(ρs[:, 1, 2] .≈ conj(ρs[:, 2, 1]))
    @test all(ρs[:, 1, 1] .+ ρs[:, 2, 2] .≈ 1.0 + 0.0im)
    @test all(norm.(ρs[:, 1, 1] .- pisetup.pisetup.complex_dat[1:pisetup.ntimes+1, 1, 1]) .< 1e-2)
    @test all(norm.(ρs[:, 1, 2] .- pisetup.pisetup.complex_dat[1:pisetup.ntimes+1, 2, 1]) .< 1e-2)
end

@testitem "adaptive-kink-QuAPI" setup=[pisetup] begin
    t, ρs = QuAPI.propagate_kink(; fbU=pisetup.U, Jw=[pisetup.Jw], pisetup.β, pisetup.ρ0, pisetup.dt, pisetup.ntimes, pisetup.svec)

    @test all(ρs[:, 1, 2] .≈ conj(ρs[:, 2, 1]))
    @test all(ρs[:, 1, 1] .+ ρs[:, 2, 2] .≈ 1.0 + 0.0im)
    @test all(norm.(ρs[:, 1, 1] .- pisetup.complex_dat[1:pisetup.ntimes+1, 1, 1]) .< 1e-2)
    @test all(norm.(ρs[:, 1, 2] .- pisetup.complex_dat[1:pisetup.ntimes+1, 2, 1]) .< 1e-2)
end

@testitem "TEMPO" setup=[pisetup] begin
    Us = TEMPO.build_augmented_propagator(; pisetup.fbU, Jw=[pisetup.Jw], pisetup.β, pisetup.dt, pisetup.ntimes, pisetup.kmax, pisetup.svec, extraargs=TEMPO.TEMPOArgs(; cutoff=1e-15, maxdim=1000))
    t_prop, ρs_prop = Utilities.apply_propagator(; propagators=Us, pisetup.ρ0, pisetup.dt, pisetup.ntimes)
    t, ρs = TEMPO.propagate(; pisetup.fbU, Jw=[pisetup.Jw], pisetup.β, pisetup.ρ0, pisetup.dt, pisetup.ntimes, pisetup.kmax, pisetup.svec, extraargs=TEMPO.TEMPOArgs(; cutoff=1e-15, maxdim=1000))

    @test all(ρs .≈ ρs_prop)
    @test all(ρs[:, 1, 2] .≈ conj(ρs[:, 2, 1]))
    @test all(ρs[:, 1, 1] .+ ρs[:, 2, 2] .≈ 1.0 + 0.0im)
    @test all(norm.(ρs[:, 1, 1] .- pisetup.complex_dat[1:pisetup.ntimes+1, 1, 1]) .< 1e-2)
    @test all(norm.(ρs[:, 1, 2] .- pisetup.complex_dat[1:pisetup.ntimes+1, 2, 1]) .< 1e-2)
end

@testitem "PCTNPI" setup=[pisetup] begin
    U = PCTNPI.build_augmented_propagator(; pisetup.fbU, Jw=[pisetup.Jw], pisetup.β, pisetup.dt, pisetup.ntimes, pisetup.svec)
    t, ρs = Utilities.apply_propagator(; propagators=U, pisetup.ρ0, pisetup.ntimes, pisetup.dt)

    @test all(ρs[:, 1, 2] .≈ conj(ρs[:, 2, 1]))
    @test all(ρs[:, 1, 1] .+ ρs[:, 2, 2] .≈ 1.0 + 0.0im)
    @test all(norm.(ρs[:, 1, 1] .- pisetup.complex_dat[1:pisetup.ntimes+1, 1, 1]) .< 1e-2)
    @test all(norm.(ρs[:, 1, 2] .- pisetup.complex_dat[1:pisetup.ntimes+1, 2, 1]) .< 1e-2)
end

@testitem "Blip" setup=[pisetup] begin
    U = Blip.build_augmented_propagator(; pisetup.fbU, Jw=[pisetup.Jw], pisetup.β, pisetup.dt, pisetup.ntimes, pisetup.svec)
    t, ρs = Utilities.apply_propagator(; propagators=U, pisetup.ρ0, pisetup.ntimes, pisetup.dt)

    @test all(ρs[:, 1, 2] .≈ conj(ρs[:, 2, 1]))
    @test all(ρs[:, 1, 1] .+ ρs[:, 2, 2] .≈ 1.0 + 0.0im)
    @test all(norm.(ρs[:, 1, 1] .- pisetup.complex_dat[1:pisetup.ntimes+1, 1, 1]) .< 1e-2)
    @test all(norm.(ρs[:, 1, 2] .- pisetup.complex_dat[1:pisetup.ntimes+1, 2, 1]) .< 1e-2)
end

@testitem "QuAPI-TTM" setup=[pisetup] begin
    t, ρs = TTM.propagate(; pisetup.fbU, Jw=[pisetup.Jw], pisetup.β, pisetup.ρ0, pisetup.dt, ntimes=pisetup.ntimes_long, rmax=pisetup.kmax + 1, pisetup.svec, path_integral_routine=QuAPI.build_augmented_propagator, extraargs=QuAPI.QuAPIArgs())

    @test all(ρs[:, 1, 2] .≈ conj(ρs[:, 2, 1]))
    @test all(ρs[:, 1, 1] .+ ρs[:, 2, 2] .≈ 1.0 + 0.0im)
    @test all(norm.(ρs[:, 1, 1] .- pisetup.complex_dat[1:pisetup.ntimes_long+1, 1, 1]) .< 1e-2)
    @test all(norm.(ρs[:, 1, 2] .- pisetup.complex_dat[1:pisetup.ntimes_long+1, 2, 1]) .< 1e-2)
end

@testitem "adaptive_kink-QuAPI-TTM" setup=[pisetup] begin
    t, ρs = TTM.propagate(; fbU=pisetup.U, Jw=[pisetup.Jw], pisetup.β, pisetup.ρ0, pisetup.dt, ntimes=pisetup.ntimes_long, rmax=pisetup.kmax + 1, pisetup.svec, path_integral_routine=QuAPI.build_augmented_propagator_kink, extraargs=QuAPI.QuAPIArgs(), forward_backward=false)

    @test all(ρs[:, 1, 2] .≈ conj(ρs[:, 2, 1]))
    @test all(ρs[:, 1, 1] .+ ρs[:, 2, 2] .≈ 1.0 + 0.0im)
    @test all(norm.(ρs[:, 1, 1] .- pisetup.complex_dat[1:pisetup.ntimes_long+1, 1, 1]) .< 1e-2)
    @test all(norm.(ρs[:, 1, 2] .- pisetup.complex_dat[1:pisetup.ntimes_long+1, 2, 1]) .< 1e-2)
end

@testitem "Blip-TTM" setup=[pisetup] begin
    t, ρs = TTM.propagate(; pisetup.fbU, Jw=[pisetup.Jw], pisetup.β, pisetup.ρ0, pisetup.dt, ntimes=pisetup.ntimes_long, rmax=pisetup.kmax + 1, pisetup.svec, path_integral_routine=Blip.build_augmented_propagator, extraargs=Blip.BlipArgs())

    @test all(ρs[:, 1, 2] .≈ conj(ρs[:, 2, 1]))
    @test all(ρs[:, 1, 1] .+ ρs[:, 2, 2] .≈ 1.0 + 0.0im)
    @test all(norm.(ρs[:, 1, 1] .- pisetup.complex_dat[1:pisetup.ntimes_long+1, 1, 1]) .< 1e-2)
    @test all(norm.(ρs[:, 1, 2] .- pisetup.complex_dat[1:pisetup.ntimes_long+1, 2, 1]) .< 1e-2)
end

@testitem "TEMPO-TTM" setup=[pisetup] begin
    t, ρs = TTM.propagate(; pisetup.fbU, Jw=[pisetup.Jw], pisetup.β, pisetup.ρ0, pisetup.dt, ntimes=pisetup.ntimes_long, rmax=pisetup.kmax + 1, pisetup.svec, path_integral_routine=TEMPO.build_augmented_propagator, extraargs=TEMPO.TEMPOArgs(; cutoff=1e-15, maxdim=1000))

    @test all(ρs[:, 1, 2] .≈ conj(ρs[:, 2, 1]))
    @test all(ρs[:, 1, 1] .+ ρs[:, 2, 2] .≈ 1.0 + 0.0im)
    @test all(norm.(ρs[:, 1, 1] .- pisetup.complex_dat[1:pisetup.ntimes_long+1, 1, 1]) .< 1e-2)
    @test all(norm.(ρs[:, 1, 2] .- pisetup.complex_dat[1:pisetup.ntimes_long+1, 2, 1]) .< 1e-2)
end

@testitem "PCTNPI-TTM" setup=[pisetup] begin
    t, ρs = TTM.propagate(; pisetup.fbU, Jw=[pisetup.Jw], pisetup.β, pisetup.ρ0, pisetup.dt, ntimes=pisetup.ntimes_long, rmax=pisetup.kmax + 1, pisetup.svec, path_integral_routine=PCTNPI.build_augmented_propagator, extraargs=PCTNPI.PCTNPIArgs())

    @test all(ρs[:, 1, 2] .≈ conj(ρs[:, 2, 1]))
    @test all(ρs[:, 1, 1] .+ ρs[:, 2, 2] .≈ 1.0 + 0.0im)
    @test all(norm.(ρs[:, 1, 1] .- pisetup.complex_dat[1:pisetup.ntimes_long+1, 1, 1]) .< 1e-2)
    @test all(norm.(ρs[:, 1, 2] .- pisetup.complex_dat[1:pisetup.ntimes_long+1, 2, 1]) .< 1e-2)
end