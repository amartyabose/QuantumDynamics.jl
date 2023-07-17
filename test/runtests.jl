using QuantumDynamics
using Test
using Plots

@testset "Isolated Hamiltonian" begin
    H = [-0.1im -1.0; -1.0 -0.5im]
    V(t) = 12*cos(10.0*t)
    
    EF = Utilities.ExternalField(V, [1.0+0.0im 0.0; 0.0 -1.0])

    ρ0 = [1.0+1.0im 0.0; 0.0 0.0]
    dt = 0.125
    ntimes = 100

    times, ρs = Bare.propagate(; Hamiltonian=H, ρ0, dt, ntimes, external_fields=[EF])
    
    plot(times, real.(ρs[:,1,1]))
end


