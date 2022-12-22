module Bare

using DifferentialEquations
using LinearAlgebra

"""
    propagate(; Hamiltonian::Matrix{ComplexF64}, ρ0::Matrix{ComplexF64}, dt::Real, ntimes::Int)
Given a potentially non-Hermitian Hamiltonian, this solves the equation of motion to propagate the input initial reduced density matrix, ρ0, with a time-step of `dt` for `ntimes` time steps.
"""
function propagate(; Hamiltonian::Matrix{ComplexF64}, ρ0::Matrix{ComplexF64}, dt::Real, ntimes::Int, cutoff=1e-10, solver=Tsit5())
    f(ρ, H, t) = -1im * (H * ρ - ρ * H')
    tspan = (0.0, ntimes * dt)
    prob = ODEProblem(f, ρ0, tspan, Hamiltonian)
    sol = solve(prob, solver, reltol=cutoff, abstol=cutoff, saveat=dt)
    sdim = size(Hamiltonian, 1)
    ρs = zeros(ComplexF64, length(sol.t), sdim, sdim)
    for j = 1:length(sol.t)
        @inbounds ρs[j, :, :] .= sol.u[j]
    end
    if !ishermitian(Hamiltonian)
        for j = 1:length(sol.t)
            @inbounds ρs[j, :, :] ./= tr(ρs[j, :, :])
        end
    end
    sol.t, ρs
end

end