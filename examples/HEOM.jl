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

function FMO(num_modes, Lmax, β, threshold=1e-10, scaled=true)
    # C. tepidum from Adolphs J, Renger T (2006) How proteins trigger excitation
    # energy transfer in the FMO complex of green sulfur bacteria. Biophys J
    # 91:2778–2797.
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
    Jw = Vector{SpectralDensities.SpectralDensity}()
    sys_ops = Vector{Matrix{ComplexF64}}()
    for (j, (λ, γ)) in enumerate(zip(λs, γs))
        push!(Jw, SpectralDensities.DrudeLorentz(; λ, γ, Δs=1.0))
        op = zeros(7, 7)
        op[j, j] = 1.0
        push!(sys_ops, op)
    end

    @time t, ρ = HEOM.propagate(; Hamiltonian=H, ρ0=ρ0, Jw, β, ntimes=nsteps, dt, sys_ops, num_modes, Lmax, scaled, threshold, extraargs=Utilities.DiffEqArgs(; reltol=1e-6, abstol=1e-6))
    t .*= au2fs
    fig, ax = new_figure("full")
    for j = 1:7
        plt.plot(t, real.(ρ[:, j, j]), lw=0.75, label="Site $(j)")
    end
    plt.legend()
    plt.xlabel(L"t (\unit{\fs})")
    plt.ylabel(L"P(t)")
    plt.savefig("heom_fmo_$(round(β;digits=2))_scaled$(scaled)_$(num_modes)_$(Lmax)_$(threshold).pdf"; bbox_inches="tight")
    plt.close()
end

function HEOM_spontaneous_emission(bo::Real, se::Real, num_modes::Int, Lmax::Int, β::Real)
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
    dt = 0.125
    ntimes = 100

    # define the jump operators corresponding to the decohering effects of the
    # individual Born-Oppenheimer surfaces.
    # bo is the coupling strength
    Jw = Vector{SpectralDensities.DrudeLorentz}()
    svec = Vector{Matrix{ComplexF64}}()
    σz = Matrix{ComplexF64}([
        1.0 0.0
        0.0 -1.0
    ])
    id = Matrix{ComplexF64}([
        1.0 0.0
        0.0 1.0
    ])
    jw1 = SpectralDensities.DrudeLorentz(; λ=bo, γ=5.0, Δs=1.0)
    bo1 = kron(σz, id)
    bo2 = kron(id, σz)

    # define the jump operators corresponding to spontaneous emission
    # se is the coupling strength
    σx = Matrix{ComplexF64}([
        0.0 1.0
        1.0 0.0
    ])
    jw3 = SpectralDensities.DrudeLorentz(; λ=se, γ=5.0, Δs=1.0)
    se3 = kron(σx, id) + kron(id, σx)

    display(dt)
    @time times, ρs = HEOM.propagate(; Hamiltonian=H, ρ0=ρ0, Jw=[jw1, jw1, jw3], β, ntimes, dt,
        sys_ops=[bo1, bo2, se3], num_modes, Lmax)

    fig, ax = new_figure("full")
    plt.plot(times, real.(ρs[:, 1, 1]), lw=0.75, label=L"\ket{ee}")
    plt.plot(times, real.(ρs[:, 2, 2]), lw=0.75, label=L"\ket{ge}")
    plt.plot(times, real.(ρs[:, 3, 3]), lw=0.75, label=L"\ket{eg}")
    plt.plot(times, real.(ρs[:, 4, 4]), lw=0.75, label=L"\ket{gg}")
    plt.plot(times, real.(ρs[:, 1, 1] + ρs[:, 2, 2] + ρs[:, 3, 3] + ρs[:, 4, 4]), lw=0.75, label="Total")
    plt.legend()
    plt.xlabel(L"t")
    plt.ylabel(L"P(t)")
    plt.savefig("heom_$(bo)_$(se)_$(β)_$(num_modes)_$(Lmax).pdf"; bbox_inches="tight")
    plt.close()
end

# println("FMO 77K")
# FMO(2, 3, 1 / (77 * 3.16683e-6), 1e-10)
# FMO(2, 3, 1 / (77 * 3.16683e-6), 1e-5)
# FMO(2, 3, 1 / (77 * 3.16683e-6), 1e-3)
# FMO(2, 3, 1 / (77 * 3.16683e-6), 1e-1)
println("FMO 300")
FMO(2, 3, 1 / (300 * 3.16683e-6), 1e-10)
FMO(2, 3, 1 / (300 * 3.16683e-6), 1e-5)
FMO(2, 3, 1 / (300 * 3.16683e-6), 1e-3)
FMO(2, 3, 1 / (300 * 3.16683e-6), 1e-1)
#println("Dimer non-commuting baths")
#HEOM_spontaneous_emission(0.25, 0.1, 3, 2, 1.0)
