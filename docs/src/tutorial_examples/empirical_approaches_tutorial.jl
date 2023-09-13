using QuantumDynamics

dt = 0.125
ntimes = 200
ρ0 = [1.0+0.0im 0.0; 0.0 0.0]

H = Utilities.create_tls_hamiltonian(; ϵ=0.0, Δ=2.0)

times, ρs = Bare.propagate(; Hamiltonian=H, ρ0, dt, ntimes)
new_figure("double")
plt.plot(times, real.(ρs[:, 1, 1]), label=L"P_1(t)")
plt.plot(times, real.(ρs[:, 2, 2]), label=L"P_2(t)")
plt.legend()
plt.xlabel(L"t")
plt.ylabel("Population")
plt.savefig("bare_system.png"; bbox_inches="tight")

L = [[0.0+0.0im 0; 0.75 0]]
times, ρs = Bare.propagate(; Hamiltonian=H, ρ0, dt, ntimes, L);
new_figure("double")
plt.plot(times, real.(ρs[:, 1, 1]), label=L"P_1(t)")
plt.plot(times, real.(ρs[:, 2, 2]), label=L"P_2(t)")
plt.legend()
plt.xlabel(L"t")
plt.ylabel("Population")
plt.savefig("lindblad.png"; bbox_inches="tight")

H = [exp(-1.5im) -0.75; -0.75 exp(0.75im)]
datum = (H[1, 1] + H[2, 2]) / 2
H[1, 1] -= datum
H[2, 2] -= datum

times, ρs = Bare.propagate(; Hamiltonian=H, ρ0, dt, ntimes);
new_figure("double")
plt.plot(times, real.(ρs[:, 1, 1]), label=L"P_1(t)")
plt.plot(times, real.(ρs[:, 2, 2]), label=L"P_2(t)")
plt.legend()
plt.xlabel(L"t")
plt.ylabel("Population")
plt.savefig("nonhermitian.png"; bbox_inches="tight")
