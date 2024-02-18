using QuantumDynamics
using LaTeXStrings
import PyPlot;
const plt = PyPlot;

invcm2au = 4.55633e-6

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

H = Matrix{ComplexF64}([
    0.0 0.000525
    0.000525 0.0
]) * invcm2au

β = 87.7899537566
svec = [1.0 -1.0]

state1 = [1.0 0.0; 0.0 0.0]
A = 1im * Utilities.commutator(H, state1)
B = copy(state1)
idmat = [1.0 0.0; 0.0 1.0]

Jw = [SpectralDensities.ExponentialCutoff(; ξ=0.1, ωc=500invcm2au)]
At, avgbond = ComplexTimePI.A_of_t(; Hamiltonian=H, β, t=0.0, N=50, Jw, svec, A=idmat)
Q = real(tr(At*state1))
time_xi01, corr_xi01, avg_bond_dim_xi01 = ComplexTimePI.correlation_function_tnpi(; Hamiltonian=H, β, tfinal=8000.0, dt=100.0, N=50, Jw, svec, A, B=[B], Z=Q, verbose=true)

Jw = [SpectralDensities.ExponentialCutoff(; ξ=0.5, ωc=500invcm2au)]
At, avgbond = ComplexTimePI.A_of_t(; Hamiltonian=H, β, t=0.0, N=50, Jw, svec, A=idmat)
Q = real(tr(At*state1))
time_xi05, corr_xi05, avg_bond_dim_xi05 = ComplexTimePI.correlation_function_tnpi(; Hamiltonian=H, β, tfinal=8000.0, dt=100.0, N=50, Jw, svec, A, B=[B], Z=Q, verbose=true)

new_figure()
plt.plot(0.002278165 * time_xi01, real.(corr_xi01) ./ H[1,2], lw=0.75, label=L"\xi=0.1")
plt.plot(0.002278165 * time_xi05, real.(corr_xi05) ./ H[1,2], lw=0.75, label=L"\xi=0.5")
plt.legend()
plt.xlabel(L"\omega_c t")
plt.ylabel(L"C_{fs}(t) / \Omega")
plt.savefig("rate_theory.pdf", bbox_inches="tight")

new_figure()
plt.plot(0.002278165 * time_xi01, avg_bond_dim_xi01, lw=0.75, label=L"\xi=0.1")
plt.plot(0.002278165 * time_xi05, avg_bond_dim_xi05, lw=0.75, label=L"\xi=0.5")
plt.legend()
plt.xlabel(L"\omega_c t")
plt.ylabel(L"\bar{\beta}(t)")
plt.savefig("bond_dim_rate_theory.pdf", bbox_inches="tight")