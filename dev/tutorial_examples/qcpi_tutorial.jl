using Revise
using QuantumDynamics

H0 = Utilities.create_tls_hamiltonian(; ϵ=2.0, Δ=2.0)
ρ0 = [1.0+0.0im 0.0; 0.0 0.0]
β = 5.0
dt = 0.25
ntimes = 100

Jw = SpectralDensities.ExponentialCutoff(; ξ=0.1, ωc=7.5)    # 1.2 Define the spectral density
ω, c = SpectralDensities.discretize(Jw, 100)
svec = [1.0, -1.0]
hb = Solvents.HarmonicBath(β, ω, c, svec, 10000);

barefbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt, ntimes)
times_QuAPI, ρs_QuAPI = TTM.propagate(; fbU=barefbU, Jw=[Jw], svec=[1.0 -1.0], β, ρ0, dt, ntimes, rmax=10, path_integral_routine=QuAPI.build_augmented_propagator, extraargs=QuAPI.QuAPIArgs())

EACP_fbU = Propagators.calculate_average_reference_propagators(; Hamiltonian=H0, solvent=hb, classical_dt=dt / 100, dt, ntimes);
times_EACP, ρs_EACP = Utilities.apply_propagator(; propagators=EACP_fbU, ρ0, ntimes, dt);
times_QCPI, ρs_QCPI = QCPI.propagate(; Hamiltonian=H0, Jw, solvent=hb, ρ0, classical_dt=dt / 100, dt, ntimes, kmax=3, extraargs=QuAPI.QuAPIArgs(), path_integral_routine=QuAPI.propagate)

new_figure("double")
plt.plot(times_EACP, real.(ρs_EACP[:, 1, 1] .- ρs_EACP[:, 2, 2]), label="EACP")
plt.plot(times_QCPI, real.(ρs_QCPI[:, 1, 1] .- ρs_QCPI[:, 2, 2]), label="QCPI")
plt.plot(times_QuAPI, real.(ρs_QuAPI[:, 1, 1] .- ρs_QuAPI[:, 2, 2]), label="QuAPI")
plt.legend()
plt.xlabel(L"t")
plt.ylabel(L"\langle\sigma_z(t)\rangle")
plt.savefig("QCPI.png"; bbox_inches="tight")