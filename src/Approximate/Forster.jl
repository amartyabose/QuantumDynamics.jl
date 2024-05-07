module Forster

using LinearAlgebra
using ..SpectralDensities, ..Utilities

function get_F_A(ω, jw, λ, ϵ, β)
    dω = ω[2] - ω[1]
    tmax = π / dω
    dt = π / (10 * ω[end])
    t = 0:dt:tmax
    α = zeros(ComplexF64, length(t))
    for (i, time) in enumerate(t)
        α[i] = 1.0 / π * Utilities.trapezoid(ω, jw .* (coth.(ω * β / 2) .* cos.(ω * time) - 1im * sin.(ω * time)))
    end
    inner_integral = zeros(ComplexF64, length(t))
    for i in 2:length(t)
        inner_integral[i] = Utilities.trapezoid(t[1:i], α[1:i])
    end
    g = zeros(ComplexF64, length(t))
    for i in 2:length(t)
        g[i] = Utilities.trapezoid(t[1:i], inner_integral[1:i])
    end
    F = exp.(-1im * (ϵ - λ) * t .- conj(g))
    A = exp.(-1im * (ϵ + λ) * t .- g)
    t, F, A
end

"""
    build_incoherent_propagator(Jws::Vector{SpectralDensities.SpectralDensityTable}, H::Matrix{ComplexF64}, dt::Float64, β::Float64; verbose::Bool=false)

Calculate the incoherent propagator and rate matrix under the approximation Forster theory.
"""
function build_incoherent_propagator(Jws::Vector{SpectralDensities.SpectralDensityTable}, H::Matrix{ComplexF64}, dt::Float64, β::Float64; verbose::Bool=false)
    nsites = length(Jws)
    F = Vector{Vector{ComplexF64}}()
    A = Vector{Vector{ComplexF64}}()
    times = Vector{Real}
    if verbose
        @info "Calculating F_j and A_j"
    end
    for i = 1:nsites
        if verbose
            @info "Calculating for unit $(i)"
        end
        ω, jw = SpectralDensities.tabulate(Jws[i], false)
        λ = SpectralDensities.reorganization_energy(Jws[i])
        t, Fi, Ai = get_F_A(ω, jw, λ, H[i, i], β)
        push!(F, copy(Fi))
        push!(A, copy(Ai))
        times = t
    end
    if verbose
        @info "Finished calculating F_j and A_j"
        @info "Calculating rate-constant matrix"
    end

    # dP / dt = K * P => dP_j / dt = \sum_{i} K_{j,i} * P_i
    k = zeros(nsites, nsites)
    for i = 1:nsites
        for j = 1:nsites
            k[j, i] = 2 * (H[i, j])^2 * real(Utilities.trapezoid(times, conj.(F[i]) .* A[j]))
        end
        k[i, i] = 0.0
    end
    for i = 1:nsites
        k[i, i] -= sum(k[:, i])
    end
    @info "Diagonalizing rate-constant matrix"
    vals, vecs = eigen(k)
    k, real.(vecs * diagm(exp.(vals * dt)) * inv(vecs))
end

end