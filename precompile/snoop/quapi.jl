using QuantumDynamics

H = Utilities.create_tls_hamiltonian(; ϵ=0.0, Δ=2.0)
ρ0 = Matrix{ComplexF64}([
    1.0 0.0
    0.0 0.0
])

Jw = SpectralDensities.ExponentialCutoff(; ξ=0.1, ωc=7.5)
SpectralDensities.reorganization_energy(Jw)

β = 1.0

dt = 0.25
ntimes = 100
barefbU = Propagators.calculate_bare_propagators(; Hamiltonian=H, dt, ntimes)

t, ρs = BlochRedfield.propagate(; Hamiltonian=H, Jw=[Jw], β, ρ0, dt, ntimes, sys_ops=[[1.0+0.0im 0.0; 0.0 -1.0]])
t, ρs = QuAPI.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, kmax=5, svec=[1.0 -1.0])
t, ρs = TEMPO.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, kmax=5, svec=[1.0 -1.0], extraargs=TEMPO.TEMPOArgs(; algorithm="naive"))
t, ρs = TEMPO.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, kmax=5, svec=[1.0 -1.0], extraargs=TEMPO.TEMPOArgs(; algorithm="densitymatrix"))
t, ρs = TEMPO.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, kmax=5, svec=[1.0 -1.0], extraargs=TEMPO.TEMPOArgs(; algorithm="fit"))
U = PCTNPI.build_augmented_propagator(; fbU=barefbU, Jw=[Jw], β, dt, ntimes=5, svec=[1.0 -1.0])
Utilities.apply_propagator(; propagators=U, ρ0, ntimes, dt)
U = Blip.build_augmented_propagator(; fbU=barefbU, Jw=[Jw], β, dt, ntimes=5, svec=[1.0 -1.0])
Utilities.apply_propagator(; propagators=U, ρ0, ntimes, dt)
t, ρs = TTM.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, rmax=5, path_integral_routine=TEMPO.build_augmented_propagator, extraargs=TEMPO.TEMPOArgs(; algorithm="naive"))
t, ρs = TTM.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, rmax=5, path_integral_routine=TEMPO.build_augmented_propagator, extraargs=TEMPO.TEMPOArgs(; algorithm="densitymatrix"))
t, ρs = TTM.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, rmax=5, path_integral_routine=TEMPO.build_augmented_propagator, extraargs=TEMPO.TEMPOArgs(; algorithm="fit"))
t, ρs = TTM.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, rmax=5, path_integral_routine=QuAPI.build_augmented_propagator, extraargs=QuAPI.QuAPIArgs())
t, ρs = TTM.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, rmax=5, path_integral_routine=QuAPI.build_augmented_propagator_parallel, extraargs=QuAPI.QuAPIArgs())
t, ρs = TTM.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, rmax=5, path_integral_routine=Blip.build_augmented_propagator, extraargs=Blip.BlipArgs())
t, ρs = TTM.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, rmax=5, path_integral_routine=Blip.build_augmented_propagator_parallel, extraargs=Blip.BlipArgs())