module BlochRedfield

using OrdinaryDiffEq
using LinearAlgebra
using ..SpectralDensities, ..Utilities

"""
    get_Rtensor(eigvals, eigvecs, Jw::AbstractVector{SpectralDensities.SpectralDensity}, svec::AbstractVector{Matrix{ComplexF64}}, β::Real)
Calculates the Bloch-Redfield R tensor given the eigenvalues, `eigvals`, and eigenvectors, `eigvecs`, of the system Hamiltonian, an inverse temperature `β`, and a number of baths specified by their spectral densities, `Jw`, and the operator through which they interact, `svec`.
"""
function get_Rtensor(eigvals, eigvecs, Jw::AbstractVector{SpectralDensities.SpectralDensity}, svec::AbstractVector{Matrix{ComplexF64}}, β::Real)
    sdim = size(eigvecs, 1)
    R = zeros(sdim, sdim, sdim, sdim)
    ω = eigvals .- eigvals'
    Γ = zeros(ComplexF64, sdim, sdim, sdim, sdim)
    @inbounds begin
        for (A, J) in zip(svec, Jw)
            A_eig = eigvecs' * A * eigvecs
            for a in axes(Γ, 1), b in axes(Γ, 2), c in axes(Γ, 3), d in axes(Γ, 4)
                Γ[a, b, c, d] += A_eig[a, b] * A_eig[c, d] * SpectralDensities.eval_spectrum(J, ω[d, c], β)
            end
        end
        for a = 1:sdim, b = 1:sdim, c = 1:sdim, d = 1:sdim
            R[a,b,c,d] -= Γ[c,b,a,d] + conj(Γ[d,a,c,b])
            for e = 1:sdim
                if b==d
                    R[a,b,c,d] += Γ[a,e,c,e]
                end
                if a==c
                    R[a,b,c,d] += conj(Γ[b,e,d,e])
                end
            end
        end
    end
    R
end

struct Params
    H::AbstractMatrix{ComplexF64}
    R
    eigvals::AbstractVector{Float64}
    EF::Union{Nothing,Vector{Utilities.ExternalField}}
end

function func_BRME(ρ, params, t)
    sdim = size(ρ, 1)
    H = deepcopy(params.H)
    if !isnothing(params.EF)
        for ef in params.EF
            H .+= ef.V(t) * ef.coupling_op
        end
    end
    dρ = -1im * Utilities.nh_commutator(H, ρ)
    for a = 1:sdim
        for b = 1:sdim
            for c = 1:sdim
                for d = 1:sdim
                    @inbounds dρ[a, b] += params.R[a, b, c, d] * ρ[c, d]
                end
            end
        end
    end
    return dρ
end

"""
    propagate(; Hamiltonian::AbstractMatrix{ComplexF64}, Jw::AbstractVector{T}, β::Real, ρ0::AbstractMatrix{ComplexF64}, dt::Real, ntimes::Int, sys_ops::Vector{Matrix{ComplexF64}}, extraargs::Utilities.DiffEqArgs=Utilities.DiffEqArgs()) where {T<:SpectralDensities.AnalyticalSpectralDensity}
Given a system Hamiltonian, the spectral densities describing the solvent, `Jw`, and an inverse temperature, this uses Bloch-Redfield Master Equations to propagate the input initial reduced density matrix, ρ0, with a time-step of `dt` for `ntimes` time steps. The ith bath, described by `Jw[i]`, interacts with the system through the operator with the values of `svec[j]`. The default solver used here is Tsit5 with a relative and absolute error cutoffs of 1e-10.
"""
function propagate(; Hamiltonian::AbstractMatrix{ComplexF64}, Jw::AbstractVector{SpectralDensities.SpectralDensity}, β::Real, ρ0::AbstractMatrix{ComplexF64}, dt::Real, ntimes::Int, sys_ops::Vector{Matrix{ComplexF64}}, external_fields::Union{Nothing,Vector{Utilities.ExternalField}}=nothing, extraargs::Utilities.DiffEqArgs=Utilities.DiffEqArgs())
    @assert ishermitian(Hamiltonian)
    for j in Jw
        @assert j isa SpectralDensities.AnalyticalSpectralDensity "Spectral density needs to be an AnalyticalSpectralDensity."
    end
    eigvals, eigvecs = eigen(Hamiltonian)
    eigvals = real.(eigvals)
    R = get_Rtensor(eigvals, eigvecs, Jw, sys_ops, β)
    H_diag = diagm(eigvals)
    params = Params(H_diag, R, eigvals, external_fields)
    ρinit = eigvecs' * ρ0 * eigvecs
    tspan = (0.0, ntimes * dt)
    prob = ODEProblem(func_BRME, ρinit, tspan, params)
    sol = solve(prob, extraargs.solver, reltol=extraargs.reltol, abstol=extraargs.abstol, saveat=dt)
    sdim = size(ρ0, 1)
    ρs = zeros(ComplexF64, length(sol.t), sdim, sdim)
    for j = 1:length(sol.t)
        @inbounds ρs[j, :, :] .= eigvecs * sol.u[j] * eigvecs'
    end
    sol.t, ρs
end

end
