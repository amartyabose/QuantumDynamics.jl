using QuantumDynamics
using LaTeXStrings
import PyPlot;
const plt = PyPlot;

thz2au = 0.0001519828500716
invcm2au = 4.55633e-6
au2fs = 0.02418884254

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

function brme_example()
    H = Matrix{ComplexF64}([
        0.0 -1.0
        -1.0 0.0
    ])
    β = 5.0
    dt = 0.25
    ntimes = 100
    ρ0 = Matrix{ComplexF64}([
        1.0 0.0
        0.0 0.0
    ])
    Jw = SpectralDensities.ExponentialCutoff(; ξ=0.1, ωc=7.5)
    time, ρs = BlochRedfield.propagate(; Hamiltonian=H, Jw=[Jw], β, ρ0, dt, ntimes, sys_ops=[[1.0+0.0im 0.0; 0.0 -1.0]])

    barefbU = Propagators.calculate_bare_propagators(; Hamiltonian=H, dt, ntimes)
    t, ρsQuAPI = QuAPI.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, kmax=9, svec=[1.0 -1.0])

    new_figure("full")
    plt.plot(t, real.(ρsQuAPI[:, 1, 1]), lw=0.75, label="Exact QuAPI results")
    plt.plot(time, real.(ρs[:, 1, 1]), lw=0.75, label="BRME")
    plt.legend()
    plt.xlabel(L"\Omega t")
    plt.ylabel(L"P(t)")
    plt.savefig("brme_example.pdf"; bbox_inches="tight")
end

brme_example()