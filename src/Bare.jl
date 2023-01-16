module Bare

using OrdinaryDiffEq
using LinearAlgebra
using ..Utilities

prop_RHS(ρ, H, t) = -1im * (H * ρ - ρ * H')

struct Hamiltonian_Linblad
    H :: Matrix{ComplexF64}
    L :: Vector{Matrix{ComplexF64}}
end
function prop_RHS_linblad(ρ, HL, t)
    dρ = -1im * Utilities.commutator(HL.H, ρ)
    for L in HL.L
        dρ .+= L * ρ * L' .- 0.5 .* L' * L * ρ .- 0.5 .* ρ * L' * L
    end
    dρ
end

"""
    propagate(; Hamiltonian::Matrix{ComplexF64}, ρ0::Matrix{ComplexF64}, dt::Real, ntimes::Int, L::Union{Nothing, Vector{Matrix{ComplexF64}}}=nothing, extraargs::Utilities.DiffEqArgs=Utilities.DiffEqArgs())
Given a potentially non-Hermitian Hamiltonian, this solves the equation of motion to propagate the input initial reduced density matrix, ρ0, with a time-step of `dt` for `ntimes` time steps. If a solution to the Linblad Master Equation is desired, make the Hamiltonian Hermitian, and keep all the non-Hermitian dissipative operators in `L`.
"""
function propagate(; Hamiltonian::Matrix{ComplexF64}, ρ0::Matrix{ComplexF64}, dt::Real, ntimes::Int, L::Union{Nothing, Vector{Matrix{ComplexF64}}}=nothing, extraargs::Utilities.DiffEqArgs=Utilities.DiffEqArgs())
    tspan = (0.0, ntimes * dt)
    prob = ODEProblem(prop_RHS, ρ0, tspan, Hamiltonian)
    if !isnothing(L)
        HL = Hamiltonian_Linblad(Hamiltonian, L)
        prob = ODEProblem(prop_RHS_linblad, ρ0, tspan, HL)
    end
    sol = solve(prob, extraargs.solver, reltol=extraargs.reltol, abstol=extraargs.abstol, saveat=dt)
    sdim = size(Hamiltonian, 1)
    ρs = zeros(ComplexF64, length(sol.t), sdim, sdim)
    for j = 1:length(sol.t)
        @inbounds ρs[j, :, :] .= sol.u[j]
    end
    # if !ishermitian(Hamiltonian)
    #     for j = 1:length(sol.t)
    #         @inbounds ρs[j, :, :] ./= tr(ρs[j, :, :])
    #     end
    # end
    sol.t, ρs
end

end