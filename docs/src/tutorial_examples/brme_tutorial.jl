@time using QuantumDynamics

@time H0 = Utilities.create_tls_hamiltonian(; ϵ=0.0, Δ=2.0)        # 1.1 Define the system Hamiltonian
@time Jw = SpectralDensities.ExponentialCutoff(; ξ=0.1, ωc=7.5)    # 1.2 Define the spectral density
β = 5.0    # 1.3 Inverse temperature
dt = 0.25
ntimes = 100
ρ0 = [1.0+0.0im 0; 0 0]

@time tarr, ρs = BlochRedfield.propagate(; Hamiltonian=H0, Jw=[Jw], β, ρ0, dt, ntimes, sys_ops=[[1.0+0.0im 0.0; 0.0 -1.0]])

fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt=dt, ntimes=ntimes)
@time t, ρs_quapi = QuAPI.propagate(; fbU=fbU, Jw=[Jw], β=β, ρ0=ρ0, dt=dt, ntimes=ntimes, kmax=7)

new_figure()
plt.plot(t, real.(ρs_quapi[:, 1, 1] .- ρs_quapi[:, 2, 2]), lw=0.75, label="QuAPI")
plt.plot(tarr, real.(ρs[:, 1, 1] .- ρs[:, 2, 2]), lw=0.75, label="BRME")
plt.legend()
plt.xlabel(L"t")
plt.ylabel(L"\langle\sigma_z(t)\rangle")
plt.savefig("BRME.png"; bbox_inches="tight")
