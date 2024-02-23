using Revise
using QuantumDynamics
using FLoops

H0 = Utilities.create_tls_hamiltonian(; ϵ=2.0, Δ=2.0)
ρ0 = [1.0+0.0im 0.0; 0.0 0.0]
β = 5.0
dt = 0.25
ntimes = 100

Jw = SpectralDensities.ExponentialCutoff(; ξ=0.1, ωc=7.5)    # 1.2 Define the spectral density
ω, c = SpectralDensities.discretize(Jw, 100)
svec = [1.0, -1.0]
hb = Solvents.HarmonicBath(β, ω, c, svec, 1000);

barefbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt, ntimes)
times_QuAPI, ρs_QuAPI = TTM.propagate(; fbU=barefbU, Jw=[Jw], svec=[1.0 -1.0], β, ρ0, dt, ntimes, rmax=10, path_integral_routine=QuAPI.build_augmented_propagator, extraargs=QuAPI.QuAPIArgs())

@floop for j = 1:80
    @reduce EACP_fbU = zero(barefbU) + Propagators.calculate_average_reference_propagators(; Hamiltonian=H0, solvent=hb, classical_dt=dt / 100, dt, ntimes)
end
times_EACP, ρs_EACP = Utilities.apply_propagator(; propagators=EACP_fbU, ρ0, ntimes, dt);
times_QCPI, ρs_QCPI = QCPI.propagate(; Hamiltonian=H0, Jw, solvent=hb, ρ0, classical_dt=dt / 100, dt, ntimes, kmax=3, extraargs=QuAPI.QuAPIArgs(), path_integral_routine=QuAPI.propagate)

function get_timestep_augmented_propagator(EACP_fbU, nsteps=1)
    ntimes = size(EACP_fbU, 1)
    sdim2 = size(EACP_fbU, 2)
    dtfbU = zeros(ComplexF64, ntimes ÷ nsteps, sdim2, sdim2)
    dtfbU[1, :, :] = EACP_fbU[nsteps, :, :]
    for j = 2:ntimes÷nsteps
        dtfbU[j, :, :] .= EACP_fbU[j*nsteps, :, :] * inv(EACP_fbU[(j-1)*nsteps, :, :])
    end
    dtfbU
end
dtfbU = get_timestep_augmented_propagator(EACP_fbU, 1)
times_blip, ρs_blip = TEMPO.propagate(; fbU=dtfbU, Jw=[Jw], svec=[1.0 -1.0], β, ρ0, dt, ntimes, kmax=10, reference_prop=true, verbose=true)
dtfbU = get_timestep_augmented_propagator(EACP_fbU, 2)
times_blip2, ρs_blip2 = TEMPO.propagate(; fbU=dtfbU, Jw=[Jw], svec=[1.0 -1.0], β, ρ0, dt=2dt, ntimes=20, kmax=10, reference_prop=true, verbose=true)

new_figure("double")
plt.plot(times_EACP, real.(ρs_EACP[:, 1, 1] .- ρs_EACP[:, 2, 2]), label="EACP")
plt.plot(times_QCPI, real.(ρs_QCPI[:, 1, 1] .- ρs_QCPI[:, 2, 2]), label="QCPI")
plt.plot(times_QuAPI, real.(ρs_QuAPI[:, 1, 1] .- ρs_QuAPI[:, 2, 2]), label="QuAPI")
plt.plot(times_blip, real.(ρs_blip[:, 1, 1] .- ρs_blip[:, 2, 2]), label="EACP++")
plt.plot(times_blip2, real.(ρs_blip2[:, 1, 1] .- ρs_blip2[:, 2, 2]), label="EACP++ 2dt")
plt.legend()
plt.xlabel(L"t")
plt.ylabel(L"\langle\sigma_z(t)\rangle")
plt.savefig("QCPI.png"; bbox_inches="tight")