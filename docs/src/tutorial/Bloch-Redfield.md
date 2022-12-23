# Bloch-Redfield Master Equation

QuantumDynamics also offers the option of simulating the dynamics of an open quantum system using the Bloch-Redfield equations. The main interface is similar to that of the path integral-based methods except for the crucial difference that instead of building on the forward-backward propagator, the Bloch-Redfield Master Equations (BRME) uses the Hamiltonian.

First, we define the system and the spectral density describing the solvent
```@example brme_eg1
using QuantumDynamics
using Plots, LaTeXStrings

H0 = [0.0+0.0im -1.0; -1.0 0.0]   # 1.1 Define the system Hamiltonian
Jw = SpectralDensities.ExponentialCutoff(; ξ=0.1, ωc=7.5)    # 1.2 Define the spectral density
β = 5.0    # 1.3 Inverse temperature
dt = 0.25
ntimes = 100
ρ0 = [1.0+0.0im 0; 0 0]
nothing
```

The interface to BRME is provided in the Bloch-Redfield module as the propagate function.
```@example brme_eg1
time, ρs = BlochRedfield.propagate(; Hamiltonian=H0, Jw=[Jw], β, ρ0, dt, ntimes, svec=[[1.0 0.0; 0.0 -1.0]]);
plot(time, real.(ρs[:,1,1] .- ρs[:,2,2]))
xlabel!(L"t")
ylabel!(L"\langle\sigma_z(t)\rangle")
```