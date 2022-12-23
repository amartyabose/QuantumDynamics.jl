# Numerically Exact Path Integral Approaches

The family of methods based on Quasi-Adiabatic Propagator Path Integral (QuAPI) is a family of numerically exact non-perturbative techniques for simulating a quantum system interacting with a harmonic environment. It simulates the reduced density matrix of an `n`-level quantum system using path integrals and the harmonic bath is incorporated through the Feynman-Vernon influence functional. The tracing out of the harmonic bath leads to a non-Markovian memory, which is used as a convergence parameter.

The most common prototypical model problem of open quantum systems is the spin-boson problem. We will illustrate the approach taken by QuantumDynamics to make the various methods compatible with each other by demonstrating how the same basic setup works for all the basic methods.

The basic steps involved for these simulations are

1. Define the system
    1. Define the Hamiltonian
    2. Define the spectral density corresponding to the solvent
    3. Specify the temperature
2. Obtain the short-time propagators that are used to construct the path integral
3. Build on top of the short-time propagators using the Feynman-Vernon influence functional.

## Example
```@example quapi_eg1
using QuantumDynamics
using Plots, LaTeXStrings

H0 = [0.0+0.0im -1.0; -1.0 0.0]   # 1.1 Define the system Hamiltonian
Jw = SpectralDensities.ExponentialCutoff(; ξ=0.1, ωc=7.5)    # 1.2 Define the spectral density
β = 5.0    # 1.3 Inverse temperature
nothing
```

Let us plot the spectral density:
```@example quapi_eg1
ω = 0:0.1:100
plot(ω, Jw.(ω), lw=2, label="")
xlabel!(L"\omega")
ylabel!(L"J(\omega)")
```

Next, we calculate the short-time forward-backward propagators, which require us to define the time-step and number of steps of simulation.
```@example quapi_eg1
dt = 0.25
ntimes = 100
fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt=dt, ntimes=ntimes)
nothing # suppress output
```

Finally, the methods incorporate the influence functional on top of the propagator. Here, we demonstrate the basic QuAPI algorithm at different memory lengths, ``kmax``.
```@example quapi_eg1
ρ0 = [1.0+0.0im 0; 0 0]
sigma_z = []
kmax = 1:2:9
time = Vector{Float64}()
for k in kmax
    @time t, ρs = QuAPI.propagate(; fbU=fbU, Jw=[Jw], β=β, ρ0=ρ0, dt=dt, ntimes=ntimes, kmax=k)
    global time = t
    push!(sigma_z, real.(ρs[:,1,1] .- ρs[:,2,2]))
end
```
```@example quapi_eg1
colors = ["red" "green" "blue" "teal" "magenta"]
plot(time, sigma_z, lw=2, label=permutedims([L"k = %$k" for k in kmax]), seriescolor=colors)
xlabel!(L"t")
ylabel!(L"\langle\sigma_z(t)\rangle")
```

Since the iteration regime can be quite costly, we have implemented an extension to the non-Markovian transfer tensor method (TTM) which is compatible with the QuAPI scheme. This is invoked in the following manner:
```@example quapi_eg1
ρ0 = [1.0+0.0im 0; 0 0]
sigma_z = []
rmax = 1:2:9
time = Vector{Float64}()
for r in rmax
    @time t, ρs = TTM.propagate(; fbU=fbU, Jw=[Jw], β=β, ρ0=ρ0, dt=dt, ntimes=ntimes, rmax=r, build_propagator=QuAPI.build_augmented_propagator)
    global time = t
    push!(sigma_z, real.(ρs[:,1,1] .- ρs[:,2,2]))
end
```
The `TTM.propagate` method, in addition to the usual arguments, takes a function to build the initial propagators for the full-path regime of the simulation. In this case, we are using QuAPI to build the propagators in the full-path segment, as indicated by `build_propagator=QuAPI.build_augmented_propagator`.
```@example quapi_eg1
colors = ["red" "green" "blue" "teal" "magenta"]
plot(time, sigma_z, lw=2, label=permutedims([L"k = %$r" for r in rmax]), seriescolor=colors)
xlabel!(L"t")
ylabel!(L"\langle\sigma_z(t)\rangle")
```

TTM can also use the so-called blip-decomposed propagators where the augmented propagators are calculated using blip-decomposed path integrals. The code remains practically identical, except the `build_propagator` argument changes from `QuAPI.build_augmented_propagator` to `Blip.build_augmented_propagator`.
```@example quapi_eg1
ρ0 = [1.0+0.0im 0; 0 0]
sigma_z = []
rmax = 1:2:9
time = Vector{Float64}()
for r in rmax
    @time t, ρs = TTM.propagate(; fbU=fbU, Jw=[Jw], β=β, ρ0=ρ0, dt=dt, ntimes=ntimes, rmax=r, build_propagator=Blip.build_augmented_propagator)
    global time = t
    push!(sigma_z, real.(ρs[:,1,1] .- ρs[:,2,2]))
end
```
```@example quapi_eg1
colors = ["red" "green" "blue" "teal" "magenta"]
plot(time, sigma_z, lw=2, label=permutedims([L"k = %$r" for r in rmax]), seriescolor=colors)
xlabel!(L"t")
ylabel!(L"\langle\sigma_z(t)\rangle")
```