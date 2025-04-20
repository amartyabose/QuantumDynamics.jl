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
bareU = Propagators.calculate_bare_propagators(; Hamiltonian=H, dt, ntimes, forward_backward=false)

N = 5
mstnpi = MSTNPI.Setup(N, 1, "Exciton")

hops = OpSum()
for (j, jnext) in zip(mstnpi.hamiltonian_indices, mstnpi.hamiltonian_indices[2:end])
    hops += 0.5, "G->E", j, "E->G", jnext
    hops += 0.5, "E->G", j, "G->E", jnext
end
fbprop = Propagators.calculate_bare_propagators(; Hamiltonian=hops, dt=0.25, ntimes=5, mstnpi)

t, ρs = BlochRedfield.propagate(; Hamiltonian=H, Jw=[Jw], β, ρ0, dt, ntimes, sys_ops=[[1.0+0.0im 0.0; 0.0 -1.0]])
t, ρs = QuAPI.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, kmax=5, svec=[1.0 -1.0])
@info "Naive algorithm"
t, ρs = TEMPO.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, kmax=5, svec=[1.0 -1.0], extraargs=TEMPO.TEMPOArgs(; algorithm="naive"))
@info "DensityMatrix algorithm"
t, ρs = TEMPO.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, kmax=5, svec=[1.0 -1.0], extraargs=TEMPO.TEMPOArgs(; algorithm="densitymatrix"))
@info "Fit algorithm"
t, ρs = TEMPO.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, kmax=5, svec=[1.0 -1.0], extraargs=TEMPO.TEMPOArgs(; algorithm="fit"))
U = PCTNPI.build_augmented_propagator(; fbU=barefbU, Jw=[Jw], β, dt, ntimes=5, svec=[1.0 -1.0])
Utilities.apply_propagator(; propagators=U, ρ0, ntimes, dt)
U = Blip.build_augmented_propagator(; fbU=barefbU, Jw=[Jw], β, dt, ntimes=5, svec=[1.0 -1.0], exec=FLoops.ThreadedEx())
U = Blip.build_augmented_propagator(; fbU=barefbU, Jw=[Jw], β, dt, ntimes=5, svec=[1.0 -1.0], exec=FLoops.SequentialEx())
U = QuAPI.build_augmented_propagator(; fbU=barefbU, Jw=[Jw], β, dt, ntimes=5, svec=[1.0 -1.0], exec=FLoops.ThreadedEx())
U = QuAPI.build_augmented_propagator(; fbU=barefbU, Jw=[Jw], β, dt, ntimes=5, svec=[1.0 -1.0], exec=FLoops.SequentialEx())
Utilities.apply_propagator(; propagators=U, ρ0, ntimes, dt)
t, ρs = TTM.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, rmax=5, path_integral_routine=TEMPO.build_augmented_propagator, extraargs=TEMPO.TEMPOArgs(; algorithm="naive"))
t, ρs = TTM.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, rmax=5, path_integral_routine=TEMPO.build_augmented_propagator, extraargs=TEMPO.TEMPOArgs(; algorithm="densitymatrix"))
t, ρs = TTM.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, rmax=5, path_integral_routine=TEMPO.build_augmented_propagator, extraargs=TEMPO.TEMPOArgs(; algorithm="fit"))
t, ρs = TTM.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, rmax=5, path_integral_routine=QuAPI.build_augmented_propagator, extraargs=QuAPI.QuAPIArgs(), exec=FLoops.SequentialEx())
t, ρs = TTM.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, rmax=5, path_integral_routine=QuAPI.build_augmented_propagator, extraargs=QuAPI.QuAPIArgs(), exec=FLoops.ThreadedEx())
t, ρs = TTM.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, rmax=5, path_integral_routine=Blip.build_augmented_propagator, extraargs=Blip.BlipArgs(), exec=FLoops.SequentialEx())
t, ρs = TTM.propagate(; fbU=barefbU, Jw=[Jw], β, ρ0, dt, ntimes, rmax=5, path_integral_routine=Blip.build_augmented_propagator, extraargs=Blip.BlipArgs(), exec=FLoops.ThreadedEx())
t, ρs = TTM.propagate(; fbU=bareU, Jw=[Jw], β, ρ0, dt, ntimes, rmax=5, path_integral_routine=Blip.build_augmented_propagator, extraargs=Blip.BlipArgs(), exec=FLoops.ThreadedEx(), forward_backward=false)

invcm2au = 4.55633e-6

H = Matrix{ComplexF64}([
    0.0 0.000525
    0.000525 0.0
]) * invcm2au

β = 87.7899537566
svec = [1.0 -1.0]

state1 = [1.0 0.0; 0.0 0.0]
A = 1im * Utilities.commutator(H, state1)
B = copy(state1)
idmat = [1.0 0.0; 0.0 1.0]

Jw = [SpectralDensities.ExponentialCutoff(; ξ=0.1, ωc=500invcm2au)]
At, avgbond = ComplexTNPI.A_of_t(; Hamiltonian=H, β, t=0.0, N=50, Jw, svec, A=idmat)
Q = real(tr(At*state1))
time_xi01, corr_xi01, avg_bond_dim_xi01 = ComplexTNPI.complex_correlation_function(; Hamiltonian=H, β, tfinal=8000.0, dt=100.0, N=50, Jw, svec, A, B=[B], Z=Q, verbose=true, extraargs=Utilities.TensorNetworkArgs())

At, num_paths = ComplexQuAPI.A_of_t(; Hamiltonian=H, β, t=0.0, N=10, Jw, svec, A=idmat)
Q = real(tr(At*state1))
time_xi01, corr_xi01, num_paths_xi01 = ComplexQuAPI.complex_correlation_function(; Hamiltonian=H, β, tfinal=8000.0, dt=100.0, N=10, Jw, svec, A, B=[B], Z=Q, verbose=true, extraargs=QuAPI.QuAPIArgs(; cutoff=1e-20), exec=ThreadedEx())