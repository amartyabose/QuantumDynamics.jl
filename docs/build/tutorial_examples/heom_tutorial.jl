using MKL
using QuantumDynamics

H0 = Utilities.create_tls_hamiltonian(; ϵ=0.0, Δ=2.0)        # 1.1 Define the system Hamiltonian
ρ0 = [1.0+0.0im 0; 0 0]
β = 0.5
Jw = SpectralDensities.DrudeLorentz(; λ=1.5, γ=7.5)
dt = 0.125
ntimes = 200

barefbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt, ntimes);
times, ρs = @time TEMPO.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, kmax=10)
times_HEOM, ρs_HEOM = @time HEOM.propagate(; Hamiltonian=H0, ρ0, β, dt, ntimes, Jw=[Jw], sys_ops=[[1.0+0.0im 0.0; 0.0 -1.0]], num_modes=1, Lmax=2)

new_figure("double")
plt.plot(times, real.(ρs[:, 1, 1] .- ρs[:, 2, 2]), "o", label="QuAPI")
plt.plot(times_HEOM, real.(ρs_HEOM[:, 1, 1] .- ρs_HEOM[:, 2, 2]), label="HEOM")
plt.legend()
plt.xlabel(L"t")
plt.ylabel(L"\expval{\sigma_z(t)}")
plt.savefig("HEOM.png"; bbox_inches="tight")