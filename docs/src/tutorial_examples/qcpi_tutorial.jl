using Revise
using QuantumDynamics

H0 = Utilities.create_tls_hamiltonian(; ϵ=2.0, Δ=2.0)
ρ0 = [1.0+0.0im 0.0; 0.0 0.0]
β = 5.0
dt = 0.25
ntimes = 100

Jw = SpectralDensities.ExponentialCutoff(; ξ=0.1, ωc=7.5)    # 1.2 Define the spectral density
@info "Discretizing"
ω, c = SpectralDensities.discretize(Jw, 100)
plt.stem(ω, c)
plt.show()
@info "Done discretizing"
svec = [1.0 -1.0]
hb = Solvents.HarmonicBath(β, [ω], [c], svec, 1000);

barefbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt, ntimes)
times_TEMPO, ρs_TEMPO = TTM.propagate(; fbU=barefbU, Jw=[Jw], svec=[1.0 -1.0], β, ρ0, dt, ntimes, rmax=10, path_integral_routine=TEMPO.build_augmented_propagator, extraargs=TEMPO.TEMPOArgs(), verbose=true)

@time EACP_fbU = Propagators.calculate_average_reference_propagators(; Hamiltonian=H0, solvent=hb, classical_dt=dt / 10, ρ0, dt, ntimes, verbose=true, reference_choice="fixed")
@info "Sampling done"
times_EACP, ρs_EACP = Utilities.apply_propagator(; propagators=EACP_fbU, ρ0, ntimes, dt);

@time ehrenfest_surface_fbU = Propagators.calculate_average_reference_propagators(; Hamiltonian=H0, solvent=hb, ρ0, classical_dt=dt / 10, dt, ntimes, verbose=true, reference_choice="ehrenfest_surface")
@info "Sampling done"
times_ehrenfest_surface, ρs_ehrenfest_surface = Utilities.apply_propagator(; propagators=ehrenfest_surface_fbU, ρ0, ntimes, dt);

@time ehrenfest_fbU = Propagators.calculate_average_reference_propagators(; Hamiltonian=H0, solvent=hb, ρ0, classical_dt=dt / 10, dt, ntimes, verbose=true, reference_choice="ehrenfest")
@info "Sampling done"
times_ehrenfest, ρs_ehrenfest = Utilities.apply_propagator(; propagators=ehrenfest_fbU, ρ0, ntimes, dt);

@time max_prob_fbU = Propagators.calculate_average_reference_propagators(; Hamiltonian=H0, solvent=hb, ρ0, classical_dt=dt / 10, dt, ntimes, verbose=true, reference_choice="max_prob")
@info "Sampling done"
times_max_prob, ρs_max_prob = Utilities.apply_propagator(; propagators=max_prob_fbU, ρ0, ntimes, dt);

@time dcsh_fbU = Propagators.calculate_average_reference_propagators(; Hamiltonian=H0, solvent=hb, ρ0, classical_dt=dt / 10, dt, ntimes, verbose=true, reference_choice="dcsh")
@info "Sampling done"
times_dcsh, ρs_dcsh = Utilities.apply_propagator(; propagators=dcsh_fbU, ρ0, ntimes, dt);

@time times_QCPI_fixed, ρs_QCPI_fixed = QCPI.propagate(; Hamiltonian=H0, Jw=[Jw], solvent=hb, ρ0, classical_dt=dt / 10, dt, ntimes, kmax=3, reference_choice="fixed", extraargs=QuAPI.QuAPIArgs(), path_integral_routine=QuAPI.propagate, verbose=true)
@time times_QCPI_ehrenfest_surface, ρs_QCPI_ehrenfest_surface = QCPI.propagate(; Hamiltonian=H0, Jw=[Jw], solvent=hb, ρ0, classical_dt=dt / 10, dt, ntimes, kmax=3, reference_choice="ehrenfest_surface", extraargs=QuAPI.QuAPIArgs(), path_integral_routine=QuAPI.propagate, verbose=true)
@time times_QCPI_dcsh, ρs_QCPI_dcsh = QCPI.propagate(; Hamiltonian=H0, Jw=[Jw], solvent=hb, ρ0, classical_dt=dt / 10, dt, ntimes, kmax=5, reference_choice="dcsh", extraargs=QuAPI.QuAPIArgs(), path_integral_routine=QuAPI.propagate, verbose=true)

new_figure("double")
plt.plot(times_TEMPO, real.(ρs_TEMPO[:, 1, 1] .- ρs_TEMPO[:, 2, 2]), label="TEMPO")
plt.plot(times_EACP, real.(ρs_EACP[:, 1, 1] .- ρs_EACP[:, 2, 2]), label="EACP")
plt.plot(times_ehrenfest, real.(ρs_ehrenfest[:, 1, 1] .- ρs_ehrenfest[:, 2, 2]), label="ehrenfest")
plt.plot(times_dcsh, real.(ρs_dcsh[:, 1, 1] .- ρs_dcsh[:, 2, 2]), label="dcsh")
plt.plot(times_ehrenfest_surface, real.(ρs_ehrenfest_surface[:, 1, 1] .- ρs_ehrenfest_surface[:, 2, 2]), label="ehrenfest_surface")
plt.plot(times_max_prob, real.(ρs_max_prob[:, 1, 1] .- ρs_max_prob[:, 2, 2]), label="max_prob")
plt.plot(times_QCPI_fixed, real.(ρs_QCPI_fixed[:, 1, 1] .- ρs_QCPI_fixed[:, 2, 2]), label="QCPI fixed reference")
plt.plot(times_TEMPO_fixed, real.(ρs_TEMPO_fixed[:, 1, 1] .- ρs_TEMPO_fixed[:, 2, 2]), label="TEMPO fixed reference")
plt.plot(times_QCPI_ehrenfest_surface, real.(ρs_QCPI_ehrenfest_surface[:, 1, 1] .- ρs_QCPI_ehrenfest_surface[:, 2, 2]), label="QCPI Ehrenfest Surface")
plt.plot(times_QCPI_dcsh, real.(ρs_QCPI_dcsh[:, 1, 1] .- ρs_QCPI_dcsh[:, 2, 2]), label="QCPI DCSH")
plt.legend()
plt.xlabel(L"t")
plt.ylabel(L"\langle\sigma_z(t)\rangle")
plt.savefig("QCPI.png"; bbox_inches="tight")