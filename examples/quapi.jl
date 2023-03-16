using QuantumDynamics
using LaTeXStrings
import PyPlot;
const plt = PyPlot;

const thz2au = 0.0001519828500716
const invcm2au = 4.55633e-6
const au2fs = 0.02418884254

function new_figure(plot_type="full")
    if plot_type == "full"
        fig = plt.figure(; figsize=(3.175, 2.25))
    elseif plot_type == "half"
        fig = plt.figure(; figsize=(1.56, 1.4))
    end
    ax = fig.add_subplot(111)

    for axis in ["top", "bottom", "left", "right"]
        ax.spines[axis].set_linewidth(0.1)
        ax.spines[axis].set_color("gray")
    end
    ax.tick_params(width=0.1, direction="in", color="gray")

    fig, ax
end

function spin_boson(ξ, ωc, β, dt, kmax, ntimes, method, name, svec=[1.0 -1.0])
    H = Utilities.create_tls_hamiltonian(; ϵ=0.0, Δ=2.0)
    ρ0 = [1.0+0.0im 0.0; 0.0 0.0]
    Jw = SpectralDensities.ExponentialCutoff(; ξ, ωc)

    barefbU = Propagators.calculate_bare_propagators(; Hamiltonian=H, dt, ntimes)
    t, ρs = method(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, kmax, svec)

    fig, ax = new_figure("full")
    plt.plot(t, real.(ρs[:, 1, 1]), lw=0.75, label=L"\Re\rho_{1,1}(t)")
    plt.plot(t, real.(ρs[:, 2, 2]), lw=0.75, label=L"\Re\rho_{2,2}(t)")
    plt.plot(t, real.(ρs[:, 1, 2]), "--", lw=0.75, label=L"\Re\rho_{1,2}(t)")
    plt.plot(t, imag.(ρs[:, 1, 2]), "--", lw=0.75, label=L"\Im\rho_{1,2}(t)")
    plt.legend()
    plt.xlabel(L"\Omega t")
    plt.ylabel(L"P(t)")
    plt.savefig("$(name)_$(β)_$(ξ)_$(ωc)_$(dt)_$(kmax).pdf"; bbox_inches="tight")
end

function molecular_dimer(bo::Real, β::Real, dt=0.125, kmax=20, ntimes=100)
    # define Hamiltonian, initial RDM and simulation parameters
    H = Matrix{ComplexF64}([
        20.0 0.0 0.0 0.0
        0.0 10.0 -1.0 0.0
        0.0 -1.0 10.0 0.0
        0.0 0.0 0.0 0.0
    ])
    ρ0 = Matrix{ComplexF64}([
        0.0 0.0 0.0 0.0
        0.0 1.0 0.0 0.0
        0.0 0.0 0.0 0.0
        0.0 0.0 0.0 0.0
    ])

    # define the jump operators corresponding to the decohering effects of the
    # individual Born-Oppenheimer surfaces.
    svec = Vector{Matrix{ComplexF64}}()
    jw1 = SpectralDensities.DrudeLorentz(; λ=bo, γ=5.0, Δs=1.0)
    svec = [1.0 1.0 -1.0 -1.0; 1.0 -1.0 1.0 -1.0]

    barefbU = Propagators.calculate_bare_propagators(; Hamiltonian=H, dt, ntimes)
    @time times, ρs = TTM.propagate(; fbU=barefbU, ρ0=ρ0, Jw=[jw1, jw1], β, ntimes, dt,
        svec, rmax=kmax, extraargs=QuAPI.QuAPIArgs(), path_integral_routine=QuAPI.build_augmented_propagator)

    fig, ax = new_figure("full")
    plt.plot(times, real.(ρs[:, 1, 1]), lw=0.75, label=L"\ket{ee}")
    plt.plot(times, real.(ρs[:, 2, 2]), lw=0.75, label=L"\ket{ge}")
    plt.plot(times, real.(ρs[:, 3, 3]), lw=0.75, label=L"\ket{eg}")
    plt.plot(times, real.(ρs[:, 4, 4]), lw=0.75, label=L"\ket{gg}")
    plt.plot(times, real.(ρs[:, 1, 1] + ρs[:, 2, 2] + ρs[:, 3, 3] + ρs[:, 4, 4]), lw=0.75, label="Total")
    plt.legend()
    plt.xlabel(L"t")
    plt.ylabel(L"P(t)")
    plt.savefig("molecular_dimer_$(β)_$(dt)_$(kmax).pdf"; bbox_inches="tight")
    plt.close()
end

function ishizaki_fleming_dimer(ntimes=200, dt=200.0, rmax=20, λ=20, γ=53.08)
    H = Matrix{ComplexF64}([
        50 100.0
        100.0 -50
    ] .* invcm2au)
    fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H, ntimes, dt)

    J1 = SpectralDensities.DrudeLorentz(; λ=λ * invcm2au, γ=γ * invcm2au, Δs=1.0)
    J2 = SpectralDensities.DrudeLorentz(; λ=λ * invcm2au, γ=γ * invcm2au, Δs=1.0)
    svec = [1.0 0.0; 0.0 1.0]
    β = 1052.0

    ρ0 = [1.0+0.0im 0.0; 0.0 0.0]
    @time times, ρs = TTM.propagate(; fbU=fbU, ρ0=ρ0, Jw=[J1, J2], β, ntimes, dt,
        svec, rmax=rmax, extraargs=TEMPO.TEMPOArgs(; cutoff=1e-13, maxdim=10000), path_integral_routine=TEMPO.build_augmented_propagator, verbose=true)

    fig, ax = new_figure("half")
    plt.plot(times .* au2fs, real.(ρs[:, 1, 1]), lw=0.75)
    plt.legend().set_visible(false)
    plt.xlabel(L"t (\unit{\fs})")
    plt.ylabel(L"P(t)")
    plt.savefig("ishizaki_fleming_dimer_$(λ)_$(γ)_$(β)_$(dt)_$(rmax).pdf"; bbox_inches="tight")
    plt.close()
end

println("Ishizaki Fleming dimer runs")
ishizaki_fleming_dimer(200, 200.0, 75, 20, 53.08)
ishizaki_fleming_dimer(200, 200.0, 75, 100, 53.08)

println("Spin-boson examples")
spin_boson(0.1, 7.5, 5, 0.25, 6, 100, QuAPI.propagate, "quapi")
spin_boson(2.0, 1.0, 1, 0.125, 150, 200, TEMPO.propagate, "tempo", [0.0 -2.0])