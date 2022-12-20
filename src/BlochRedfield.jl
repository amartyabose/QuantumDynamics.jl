module BlochRedfield

using DifferentialEquations
using LinearAlgebra, Tullio
using ..SpectralDensities

function get_Rtensor(eigvals, eigvecs, Jw::Vector{T}, svec::Vector{Matrix{Float64}}, β::Real) where {T<:SpectralDensities.AnalyticalSpectralDensity}
    sdim = size(eigvecs, 1)
    R = zeros(4, sdim, sdim, sdim, sdim)
    for (A, J) in zip(svec, Jw)
        A_eig = eigvecs * A * inv(eigvecs)
        for a = 1:sdim
            for b = 1:sdim
                for c = 1:sdim
                    for d = 1:sdim
                        @inbounds R[1, a, b, c, b] -= 0.5 * ( A_eig[a, d] * A_eig[d, c] * SpectralDensities.eval_spectrum(J, eigvals[c]-eigvals[d], β) )
                        @inbounds R[2, a, b, c, d] += 0.5 * ( A_eig[a, c] * A_eig[d, b] * SpectralDensities.eval_spectrum(J, eigvals[c]-eigvals[a], β) )
                        @inbounds R[3, a, b, a, d] -= 0.5 * ( A_eig[d, c] * A_eig[c, b] * SpectralDensities.eval_spectrum(J, eigvals[d]-eigvals[c], β) )
                        @inbounds R[4, a, b, c, d] += 0.5 * ( A_eig[a, c] * A_eig[d, b] * SpectralDensities.eval_spectrum(J, eigvals[d]-eigvals[b], β) )
                    end
                end
            end
        end
    end
    R
end

struct Params
    H :: Matrix{ComplexF64}
    R
    eigvals :: Vector{Float64}
end

function func_nosec(ρ, params, t) 
    sdim = size(ρ, 1)
    dρ = -1im * (params.H * ρ - ρ * params.H)
    for a = 1:sdim
        for b = 1:sdim
            for c = 1:sdim
                for d = 1:sdim
                    @inbounds dρ[a,b] += params.R[1,a,b,c,d] * exp(1im * (params.eigvals[a] - params.eigvals[c]) * t) * ρ[c,d] + params.R[2,a,b,c,d] * exp(1im * (params.eigvals[a] - params.eigvals[c] + params.eigvals[d] - params.eigvals[b]) * t) * ρ[c,d] + params.R[3,a,b,c,d] * exp(1im * (params.eigvals[a] - params.eigvals[c] + params.eigvals[d] - params.eigvals[b]) * t) * ρ[c,d] + params.R[4,a,b,c,d] * exp(1im * (params.eigvals[d] - params.eigvals[b]) * t) * ρ[c,d]
                end
            end
        end
    end
    return dρ
end

function func_sec(ρ, params, t) 
    sdim = size(ρ, 1)
    dρ = -1im * (params.H * ρ - ρ * params.H)
    for a = 1:sdim
        for b = 1:sdim
            for c = 1:sdim
                for d = 1:sdim
                    @inbounds dρ[a,b] += params.R[1,a,b,c,d] * ρ[c,d] + params.R[2,a,b,c,d] * ρ[c,d] + params.R[3,a,b,c,d] * ρ[c,d] + params.R[4,a,b,c,d] * ρ[c,d]
                end
            end
        end
    end
    return dρ
end

function propagate(; Hamiltonian::Matrix{ComplexF64}, Jw::Vector{T}, β::Real, ρ0::Matrix{ComplexF64}, dt::Real, ntimes::Int, svec::Vector{Matrix{Float64}}, solver=Tsit5(), secular::Bool=false) where {T<:SpectralDensities.AnalyticalSpectralDensity}
    eigvals, eigvecs = eigen(Hamiltonian)
    R = get_Rtensor(eigvals, eigvecs, Jw, svec, β)
    H_diag = diagm(eigvals)
    params = Params(H_diag, R, eigvals)
    ρinit = inv(eigvecs) * ρ0 * eigvecs
    tspan = (0.0, ntimes * dt)
    prob = ODEProblem(func_nosec, ρinit, tspan, params)
    if secular
        prob = ODEProblem(func_sec, ρinit, tspan, params)
    end
    sol = solve(prob, solver, saveat=dt)
    sdim = size(ρ0, 1)
    ρs = zeros(ComplexF64, length(sol.t), sdim, sdim)
    for j = 1:length(sol.t)
        @inbounds ρs[j, :, :] .= eigvecs * sol.u[j] * inv(eigvecs)
    end
    sol.t, ρs
end

end