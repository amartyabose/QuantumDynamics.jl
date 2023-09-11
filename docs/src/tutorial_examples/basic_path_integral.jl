@time using QuantumDynamics

@time H0 = Utilities.create_tls_hamiltonian(; ϵ=0.0, Δ=2.0)        # 1.1 Define the system Hamiltonian
@time Jw = SpectralDensities.ExponentialCutoff(; ξ=0.1, ωc=7.5)    # 1.2 Define the spectral density
β = 5.0    # 1.3 Inverse temperature

ω = 0:0.1:100
new_figure()
plt.plot(ω, Jw.(ω), lw=0.75, "k")
plt.xlabel(L"\omega")
plt.ylabel(L"J(\omega)")
plt.savefig("spectral_density.png"; bbox_inches="tight")

dt = 0.25
ntimes = 100
fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt=dt, ntimes=ntimes)

ρ0 = [1.0+0.0im 0; 0 0]

println("QuAPI")
sigma_z = []
kmax = [2, 5, 9]
time = Vector{Float64}()
for k in kmax
    @show k
    @time t, ρs = QuAPI.propagate(; fbU=fbU, Jw=[Jw], β=β, ρ0=ρ0, dt=dt, ntimes=ntimes, kmax=k)
    global time = t
    push!(sigma_z, real.(ρs[:, 1, 1] .- ρs[:, 2, 2]))
end

new_figure()
for (j, k) in enumerate(kmax)
    plt.plot(time, sigma_z[j], lw=0.75, label=L"k = %$k")
end
plt.legend()
plt.xlabel(L"t")
plt.ylabel(L"\expval{\sigma_z(t)}")
plt.savefig("QuAPI.png"; bbox_inches="tight")

println("TEMPO")
sigma_z_TEMPO = []
kmax = [2, 5, 9]
time = Vector{Float64}()
for k in kmax
    @show k
    @time t, ρs = TEMPO.propagate(; fbU=fbU, Jw=[Jw], β=β, ρ0=ρ0, dt=dt, ntimes=ntimes, kmax=k)
    global time = t
    push!(sigma_z_TEMPO, real.(ρs[:, 1, 1] .- ρs[:, 2, 2]))
end

new_figure()
for (j, k) in enumerate(kmax)
    plt.plot(time, sigma_z_TEMPO[j], lw=0.75, label=L"k = %$k")
end
plt.legend()
plt.xlabel(L"t")
plt.ylabel(L"\expval{\sigma_z(t)}")
plt.savefig("TEMPO.png"; bbox_inches="tight")

println("QuAPI TTM")
sigma_z = []
rmax = [2, 5, 9]
time = Vector{Float64}()
for r in rmax
    @show r
    @time t, ρs = TTM.propagate(; fbU=fbU, Jw=[Jw], β=β, ρ0=ρ0, dt=dt, ntimes=ntimes, rmax=r, path_integral_routine=QuAPI.build_augmented_propagator, extraargs=QuAPI.QuAPIArgs())
    global time = t
    push!(sigma_z, real.(ρs[:, 1, 1] .- ρs[:, 2, 2]))
end

new_figure()
for (j, k) in enumerate(kmax)
    plt.plot(time, sigma_z[j], lw=0.75, label=L"k = %$k")
end
plt.legend()
plt.xlabel(L"t")
plt.ylabel(L"\expval{\sigma_z(t)}")
plt.savefig("QuAPI-TTM.png"; bbox_inches="tight")