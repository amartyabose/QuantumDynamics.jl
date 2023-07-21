# Dynamics in presence of an external light

It has been shown that dissipative tunneling dynamics can be controlled by continuous wave light. We replicate some of the results here. 

As usual, first, we set up the system:
```@example external_eg1
using QuantumDynamics
using Plots, LaTeXStrings

H0 = Utilities.create_tls_hamiltonian(; ϵ=0.0, Δ=2.0)        # 1.1 Define the system Hamiltonian
V(t) = 11.96575 * cos(10.0 * t)   # This is the monochromatic light
EF = Utilities.ExternalField(V, [1.0+0.0im 0.0; 0.0 -1.0])
Jw = SpectralDensities.ExponentialCutoff(; ξ=0.16, ωc=7.5)    # 1.2 Define the spectral density
β = 0.5    # 1.3 Inverse temperature
nothing
```

Calculate the forward-backward propagators. For the case with the external field, we use the `Propagators.calculate_bare_propagators_external_field` function.
```@example external_eg1
dt = 0.125
ntimes = 100
fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt=dt, ntimes=ntimes, external_fields=[EF])
nofield_fbU = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt=dt, ntimes=ntimes)
nothing # suppress output
```

Simulate the system with the field. TTM does not yet work with time-dependent Hamiltonians. So, we resort to plain QuAPI.
```@example external_eg1
ρ0 = [1.0+0.0im 0; 0 0]
sigma_z = []
kmax = [2,5,9]
time = Vector{Float64}()
for k in kmax
    @time t, ρs = QuAPI.propagate(; fbU=fbU, Jw=[Jw], β=β, ρ0=ρ0, dt=dt, ntimes=ntimes, kmax=k)
    global time = t
    push!(sigma_z, real.(ρs[:,1,1] .- ρs[:,2,2]))
end
```

Use TTM to simulate the case without the external field.
```@example external_eg1
sigma_z_nofield = []
kmax = [2,5,9]
time = Vector{Float64}()
for k in kmax
    @time t, ρs = TTM.propagate(; fbU=nofield_fbU, Jw=[Jw], β=β, ρ0=ρ0, dt=dt, ntimes=ntimes, rmax=k, extraargs=QuAPI.QuAPIArgs(), path_integral_routine=QuAPI.build_augmented_propagator_QuAPI_TTM, QuAPI=true)
    global time = t
    push!(sigma_z_nofield, real.(ρs[:,1,1] .- ρs[:,2,2]))
end
```

Obtain the Markovian dynamics in presence of light but in absence of the dissipative medium. 
```@example external_eg1
time, ρs_nodissip = Utilities.apply_propagator(; propagators=fbU, ρ0=ρ0, ntimes=ntimes, dt=dt)
nothing
```

Plot the results.
```@example external_eg1
colors = ["red" "green" "blue" "teal" "magenta"]
plot()
for (j, k) in enumerate(kmax)
    plot!(time, sigma_z[j], lw=2, label="Light " * L"k = %$k", seriescolor=colors[j])
    plot!(time, sigma_z_nofield[j], lw=2, ls=:dash, label="No light " * L"k = %$k", seriescolor=colors[j])
end
# plot(time, sigma_z, lw=2, label=permutedims(["Light k = $k" for k in kmax]), seriescolor=colors)
# plot!(time, sigma_z_nofield, lw=2, ls=:dash, label=permutedims(["No light r = $k" for k in kmax]), seriescolor=colors)
plot!(time, real.(ρs_nodissip[:,1,1] .- ρs_nodissip[:,2,2]), lw=2, label="No dissipation")
xlabel!(L"t")
ylabel!(L"\langle\sigma_z(t)\rangle")
```

The localization phenomenon, though not as pronounced as in absence of dissipative media, is still clearly visible. As a comparison, we also simulate the dynamics in presence of a light pulse
```@example external_eg1
V1(t) = 11.96575 * cos(10.0 * t) * exp(-t^2 / 8)   # This is the light pulse
EF1 = Utilities.ExternalField(V1, [1.0+0.0im 0.0; 0.0 -1.0])
fbU_pulse = Propagators.calculate_bare_propagators(; Hamiltonian=H0, dt=dt, ntimes=ntimes, external_fields=[EF1])
kmax = [2,5,9]
@time time, ρs = QuAPI.propagate(; fbU=fbU_pulse, Jw=[Jw], β=β, ρ0=ρ0, dt=dt, ntimes=ntimes, kmax=9)
sigma_z_pulse = real.(ρs[:,1,1] .- ρs[:,2,2])
nothing
```

Plot the results.
```@example external_eg1
colors = ["black"]
plot(time, sigma_z_nofield[end], lw=2, ls=:dashdotdot, label="No light", seriescolor=colors[1])
plot!(time, sigma_z[end], lw=2, label="CW k = $(kmax[end])", seriescolor=colors[1])
plot!(time, sigma_z_pulse, lw=2, ls=:dash, label="Pulse k = 9", seriescolor=colors[1])
xlabel!(L"t")
ylabel!(L"\langle\sigma_z(t)\rangle")
```
