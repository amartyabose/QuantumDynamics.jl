using QuantumDynamics

H0 = Utilities.create_tls_hamiltonian(; ϵ=0.0, Δ=2.0)        # 1.1 Define the system Hamiltonian
V(t) = 11.96575 * cos(10.0 * t)   # This is the monochromatic light
EF = Utilities.ExternalField(V, [1.0+0.0im 0.0; 0.0 -1.0])
Jw = SpectralDensities.ExponentialCutoff(; ξ=0.16, ωc=7.5)    # 1.2 Define the spectral density
β = 0.5    # 1.3 Inverse temperature

dt = 0.125
ntimes = 100
fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt=dt, ntimes=ntimes, external_fields=[EF])
nofield_fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt=dt, ntimes=ntimes)

ρ0 = [1.0+0.0im 0; 0 0]
tlight, ρlight = QuAPI.propagate(; fbU=fbU, Jw=[Jw], β=β, ρ0=ρ0, dt=dt, ntimes=ntimes, kmax=9)
t, ρ = TTM.propagate(; fbU=nofield_fbU, Jw=[Jw], β=β, ρ0=ρ0, dt=dt, ntimes=ntimes, rmax=9, extraargs=QuAPI.QuAPIArgs(), path_integral_routine=QuAPI.build_augmented_propagator)
t_nodissip, ρ_nodissip = Utilities.apply_propagator(; propagators=fbU, ρ0=ρ0, ntimes=ntimes, dt=dt)

new_figure("double")
plt.plot(t, real.(ρ[:, 1, 1] .- ρ[:, 2, 2]), label="No light")
plt.plot(tlight, real.(ρlight[:, 1, 1] .- ρlight[:, 2, 2]), label="With light")
plt.plot(t_nodissip, real.(ρ_nodissip[:, 1, 1] .- ρ_nodissip[:, 2, 2]), label="Without dissipation")
# plt.grid(true)
plt.legend()
plt.xlabel(L"t")
plt.ylabel(L"\langle\sigma_z(t)\rangle")
plt.savefig("light.png"; bbox_inches="tight")

V1(t) = 11.96575 * cos(10.0 * t) * exp(-t^2 / 8)   # This is the light pulse
EF1 = Utilities.ExternalField(V1, [1.0+0.0im 0.0; 0.0 -1.0])
fbU_pulse = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt=dt, ntimes=ntimes, external_fields=[EF1])
@time time_pulse, ρs = QuAPI.propagate(; fbU=fbU_pulse, Jw=[Jw], β=β, ρ0=ρ0, dt=dt, ntimes=ntimes, kmax=9)
sigma_z_pulse = real.(ρs[:, 1, 1] .- ρs[:, 2, 2])

new_figure("double")
plt.plot(t, real.(ρ[:, 1, 1] .- ρ[:, 2, 2]), label="No light")
plt.plot(tlight, real.(ρlight[:, 1, 1] .- ρlight[:, 2, 2]), label="CW")
plt.plot(time_pulse, sigma_z_pulse, label="Pulse")
plt.legend()
plt.xlabel(L"t")
plt.ylabel(L"\langle\sigma_z(t)\rangle")
plt.savefig("light_pulse.png"; bbox_inches="tight")