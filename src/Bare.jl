module Bare

using OrdinaryDiffEq
using LinearAlgebra
using ..Utilities

struct SysParams
    H::AbstractMatrix{ComplexF64}
    L::Union{Nothing,Vector{Matrix{ComplexF64}}}
    EF::Union{Nothing,Vector{Utilities.ExternalField}}
end
function prop_RHS(ρ, HL, t)
    H = deepcopy(HL.H)
    if !isnothing(HL.EF)
        for ef in HL.EF
            H .+= ef.V(t) * ef.coupling_op
        end
    end
    dρ = -1im * (H * ρ - ρ * H')
    if !isnothing(HL.L)
        for L in HL.L
            dρ .+= L * ρ * L' .- 0.5 .* L' * L * ρ .- 0.5 .* ρ * L' * L
        end
    end
    dρ
end

"""
    propagate(; Hamiltonian::AbstractMatrix{ComplexF64}, ρ0::AbstractMatrix{ComplexF64}, dt::Real, ntimes::Int, L::Union{Nothing,Vector{Matrix{ComplexF64}}}=nothing, external_fields::Union{Nothing,Vector{Utilities.ExternalField}}=nothing, extraargs::Utilities.DiffEqArgs=Utilities.DiffEqArgs())
Given a potentially non-Hermitian Hamiltonian, this solves the equation of motion to propagate the input initial reduced density matrix, ρ0, with a time-step of `dt` for `ntimes` time steps. If a solution to the Lindblad Master Equation is desired, make the Hamiltonian Hermitian, and keep all the non-Hermitian dissipative operators in `L`.
"""
function propagate(; Hamiltonian::AbstractMatrix{ComplexF64}, ρ0::AbstractMatrix{ComplexF64}, dt::Real, ntimes::Int, L::Union{Nothing,Vector{Matrix{ComplexF64}}}=nothing, external_fields::Union{Nothing,Vector{Utilities.ExternalField}}=nothing, extraargs::Utilities.DiffEqArgs=Utilities.DiffEqArgs())
    tspan = (0.0, ntimes * dt)
    params = SysParams(Hamiltonian, L, external_fields)
    prob = ODEProblem(prop_RHS, ρ0, tspan, params)
    sol = solve(prob, extraargs.solver, reltol=extraargs.reltol, abstol=extraargs.abstol, saveat=dt)
    sdim = size(Hamiltonian, 1)
    ρs = zeros(ComplexF64, length(sol.t), sdim, sdim)
    for j = 1:length(sol.t)
        @inbounds ρs[j, :, :] .= sol.u[j]
    end
    sol.t, ρs
end

end