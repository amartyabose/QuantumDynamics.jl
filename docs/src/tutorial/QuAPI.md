# Quasi-Adiabatic Propagator Path Integral

Quasi-Adiabatic Propagator Path Integral (QuAPI) is one of the foundational numerically exact techniques for simulating a quantum system interacting with a harmonic environment. It simulates the reduced density matrix of an `n`-level quantum system using path integrals and the harmonic bath is incorporated through the Feynman-Vernon influence functional. The tracing out of the harmonic bath leads to a non-Markovian memory, which is used as a convergence parameter.

# Using the QuAPI module of QuantumDynamics
To use the QuAPI model, one needs to first set up the system and the harmonic bath. Suppose we are simulating a classical spin-boson Hamiltonian, then the system Hamiltonian is given as
```@example quapi_eg1
H0 = [0.0 -1.0; -1.0 0.0]
```
and the solvent Hamiltonian is described by a spectral density chosen here to be the Ohmic spectral density with an exponential cutoff held at an inverse temperature `\beta = 5.0`:
```@example quapi_eg1
using QuantumDynamics
Jw = SpectralDensities.ExponentialCutoff(ξ=0.1, ωc=7.5)
β = 5.0
```

With the main system having been specified, we need to pick a proper time-step, which is a convergence parameter.
```@example quapi_eg1
dt = 0.25
ntimes = 100
time = 0:dt:dt*ntimes;
```

We simulate the dynamics for a sequence of increasing non-Markovian memory lengths, `kmax`:
```@example quapi_eg1
ρ0 = [1.0 0; 0 0]
sigma_z = []
kmax = 1:2:9
for k in kmax
    ρs = QuAPI.propagate(; Hamiltonian=H0, Jw, β, ρ0, dt, ntimes, kmax=k)
    push!(sigma_z, real(ρs[1,1,:] - ρs[2,2,:]))
end
```

The dynamics can now be simply plotted.
```@example quapi_eg1
using Plots, LaTeXStrings
colors = ["red" "green" "blue" "teal" "magenta"]
plot(time, sigma_z, lw = 2, label=permutedims([L"k = %$k" for k in kmax]), size=(800, 600), seriescolor=colors)
xlabel!(L"t")
ylabel!(L"\langle\sigma_z(t)\rangle")
```