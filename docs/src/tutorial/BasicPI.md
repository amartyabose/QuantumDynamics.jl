# Numerically Exact Path Integral Approaches

The family of methods based on Quasi-Adiabatic Propagator Path Integral (QuAPI) is a family of numerically exact non-perturbative techniques for simulating a quantum system interacting with a harmonic environment. It simulates the reduced density matrix of an `n`-level quantum system using path integrals and the harmonic bath is incorporated through the Feynman-Vernon influence functional. The tracing out of the harmonic bath leads to a non-Markovian memory, which is used as a convergence parameter.

While at a first glance, the restriction to harmonic environments may seem arbitrarily limiting, it is actually quite general. Under the Gaussian response theory, when the environment is large and has enough "independent" degrees of freedom, the impact of an atomistically-described environment can be mapped onto a bath of harmonic oscillators with given frequencies and coupling strengths. Together these frequencies and couplings are described through the spectral density of the solvent which is given by
```math
J(\omega) = \frac{\pi}{2} \sum_j \frac{c_j^2}{\omega_j}\delta\left(\omega-\omega_j\right).
```
The full system-harmonic environment Hamiltonian is then given by
```math
\hat{H} = \hat{H}_0 + \hat{H}_\text{env}\\
\hat{H}_\text{env} = \sum_j \frac{p_j^2}{2} + \frac{1}{2}\omega_j^2\left(x_j - \frac{c_j}{\omega_j^2}\hat{s}\right)^2,
```
where ``\hat{s}`` is the system operator that interacts with the environment.

The most common prototypical model problem of open quantum systems is the spin-boson problem. We will illustrate the approach taken by QuantumDynamics to make the various methods compatible with each other by demonstrating how the same basic setup works for all the basic methods.

The basic steps involved for these simulations are

1. Define the system
    1. Define the Hamiltonian
    2. Define the spectral density corresponding to the solvent
    3. Specify the temperature
2. Obtain the short-time propagators that are used to construct the path integral
3. Build on top of the short-time propagators using the Feynman-Vernon influence functional.

In this tutorial, we will show how to use the different methods of the QuAPI family to obtain results for a single parameter. This side-by-side use of all the algorithms serve to emphasize the similarity of the APIs involved.

```@example quapi_eg1
using QuantumDynamics
using Plots, LaTeXStrings

H0 = Utilities.create_tls_hamiltonian(; ϵ=0.0, Δ=2.0)        # 1.1 Define the system Hamiltonian
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

## Iterative Quasi-Adiabatic propagator Path Integral (QuAPI)
Finally, the methods incorporate the influence functional on top of the propagator. First, we demonstrate the basic QuAPI algorithm ([QuAPI review](https://doi.org/10.1063/1.531046)) at different memory lengths, `kmax`. The exact method can also be used with filtering if the optional argument of `extraargs` of type `QuAPI.QuAPIArgs` is provided.
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

## Time-Evolved Matrix Product Operator (TEMPO)
Recently ideas of tensor network have been used to make path integral calculations more efficient. The correlation between the time-points decrease with the temporal separation between them. This allows for significantly compressed matrix product state (MPS) representation of the so-called path-amplitude tensor. The influence functional is represented as a matrix product operator and applied to this path-amplitude MPS to incorporate the effect of the baths. The interface is kept consistent with the other path integral methods like QuAPI. The MPO-MPS applications is controlled through a `cutoff` threshold and a `maxdim` threshold. The method used for applying an MPO to an MPS can be chosen to be one of `naive` and `densitymatrix`. These settings are passed as `extraargs`, which is an object of `TNPI.TNPIArgs`. By default, `cutoff=1e-8`, `maxdim=50` and `method=naive`. These ideas have been outlined in [TEMPO](https://dx.doi.org/10.1038/s41467-018-05617-3). The implementation follows the details of [TNPI](https://arxiv.org/abs/2106.12523) incorporating multiple baths and the QuAPI splitting.

```@example quapi_eg1
ρ0 = [1.0+0.0im 0; 0 0]
sigma_z_TEMPO = []
kmax = 2:2:9
time = Vector{Float64}()
for k in kmax
    @time t, ρs = TEMPO.propagate(; fbU=fbU, Jw=[Jw], β=β, ρ0=ρ0, dt=dt, ntimes=ntimes, kmax=k)
    global time = t
    push!(sigma_z_TNPI, real.(ρs[:,1,1] .- ρs[:,2,2]))
end
```
```@example quapi_eg1
colors = ["red" "green" "blue" "teal" "magenta"]
plot(time, sigma_z_TEMPO, lw=2, label=permutedims([L"k = %$k" for k in kmax]), seriescolor=colors[2:end])
xlabel!(L"t")
ylabel!(L"\langle\sigma_z(t)\rangle")
```

## Transfer Tensor Method with QuAPI and Blips
Since the iteration regime can be quite costly, we have implemented an extension to the non-Markovian transfer tensor method (TTM) ([TTM](https://link.aps.org/doi/10.1103/PhysRevLett.112.110401)) which is compatible with the QuAPI scheme. This is invoked in the following manner:
```@example quapi_eg1
ρ0 = [1.0+0.0im 0; 0 0]
sigma_z = []
rmax = 1:2:9
time = Vector{Float64}()
for r in rmax
    @time t, ρs = TTM.propagate(; fbU=fbU, Jw=[Jw], β=β, ρ0=ρ0, dt=dt, ntimes=ntimes, rmax=r, extraargs=QuAPI.QuAPIArgs(), path_integral_routine=QuAPI.build_augmented_propagator)
    global time = t
    push!(sigma_z, real.(ρs[:,1,1] .- ρs[:,2,2]))
end
```
The `TTM.propagate` method, in addition to the usual arguments, takes a function to build the initial propagators for the full-path regime of the simulation. In this case, we are using QuAPI to build the propagators in the full-path segment, as indicated by `path_integral_routine=QuAPI.build_augmented_propagator`. Other possible choices are `path_integral_routine=Blip.build_augmented_propagator` and `path_integral_routine=TNPI.build_augmented_propagator`. Also notice that because each of these `path_integral_routine`s take different `extraargs`, it is not possible to provide a default. Here, it is necessary for the `extraargs` to be provided and it needs to be consistent with the `path_integral_routine`.
```@example quapi_eg1
colors = ["red" "green" "blue" "teal" "magenta"]
plot(time, sigma_z, lw=2, label=permutedims([L"k = %$r" for r in rmax]), seriescolor=colors)
xlabel!(L"t")
ylabel!(L"\langle\sigma_z(t)\rangle")
```

TTM can also use the so-called blip-decomposed propagators where the augmented propagators are calculated using blip-decomposed path integrals. The code remains practically identical, except the `path_integral_routine` argument changes from `QuAPI.build_augmented_propagator` to `Blip.build_augmented_propagator`.
```@example quapi_eg1
ρ0 = [1.0+0.0im 0; 0 0]
sigma_z = []
rmax = 1:2:9
time = Vector{Float64}()
for r in rmax
    @time t, ρs = TTM.propagate(; fbU=fbU, Jw=[Jw], β=β, ρ0=ρ0, dt=dt, ntimes=ntimes, rmax=r, extraargs=Blip.BlipArgs(), path_integral_routine=Blip.build_augmented_propagator)
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