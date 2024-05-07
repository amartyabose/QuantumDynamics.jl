module Spectroscopy

using LinearAlgebra
using ..GQME, ..Utilities

μp_ρ_commutator(ρ::AbstractMatrix{<:Complex}, μp::AbstractArray{<:Real}, _::AbstractArray{<:Real}) = Utilities.commutator(μp, ρ)
μm_ρ_commutator(ρ::AbstractMatrix{<:Complex}, _::AbstractArray{<:Real}, μm::AbstractArray{<:Real}) = Utilities.commutator(μm, ρ)
left_multiply_by_μm(ρ::AbstractMatrix{<:Complex}, _::AbstractArray{<:Real}, μm::AbstractArray{<:Real}) = μm * ρ
left_multiply_by_μp(ρ::AbstractMatrix{<:Complex}, μp::AbstractArray{<:Real}, _::AbstractArray{<:Real}) = μp * ρ
right_multiply_by_μm(ρ::AbstractMatrix{<:Complex}, _::AbstractArray{<:Real}, μm::AbstractArray{<:Real}) = (-1) * ρ * μm
right_multiply_by_μp(ρ::AbstractMatrix{<:Complex}, μp::AbstractArray{<:Real}, _::AbstractArray{<:Real}) = (-1) * ρ * μp

const rephasing_2D = [μm_ρ_commutator, μp_ρ_commutator, μp_ρ_commutator, left_multiply_by_μm]
const nonrephasing_2D = [μp_ρ_commutator, μm_ρ_commutator, μp_ρ_commutator, left_multiply_by_μm]

const rephasing_gsb = [right_multiply_by_μm, right_multiply_by_μp, left_multiply_by_μp, left_multiply_by_μm]
const nonrephasing_gsb = [left_multiply_by_μp, left_multiply_by_μm, left_multiply_by_μp, left_multiply_by_μm]

const rephasing_se = [right_multiply_by_μm, left_multiply_by_μp, right_multiply_by_μp, left_multiply_by_μm]
const nonrephasing_se = [left_multiply_by_μp, right_multiply_by_μm, right_multiply_by_μp, left_multiply_by_μm]

const rephasing_esa = [right_multiply_by_μm, left_multiply_by_μp, left_multiply_by_μp, left_multiply_by_μm]
const nonrephasing_esa = [left_multiply_by_μp, right_multiply_by_μm, left_multiply_by_μp, left_multiply_by_μm]

# rephasing2D = [(ρ, μp, μm) -> Utilities.commutator(μm, ρ), (ρ, μp, μm) -> Utilities.commutator(μp, ρ), (ρ, μp, μm) -> Utilities.commutator(μp, ρ), (ρ, μp, μm) -> μm * ρ]
# nonrephasing2D = [(ρ, μp, μm) -> Utilities.commutator(μp, ρ), (ρ, μp, μm) -> Utilities.commutator(μm, ρ), (ρ, μp, μm) -> Utilities.commutator(μp, ρ), (ρ, μp, μm) -> μm * ρ]
# 
# rephasing_gsb = [(ρ, μp, μm) -> ρ * μm, (ρ, μp, μm) -> ρ * μp, (ρ, μp, μm) -> μp * ρ, (ρ, μp, μm) -> μm * ρ]
# nonrephasing_gsb = [(ρ, μp, μm) -> μp * ρ, (ρ, μp, μm) -> μm * ρ, (ρ, μp, μm) -> μp * ρ, (ρ, μp, μm) -> μm * ρ]
# 
# rephasing_se = [(ρ, μp, μm) -> ρ * μm, (ρ, μp, μm) -> μp * ρ, (ρ, μp, μm) -> ρ * μp, (ρ, μp, μm) -> μm * ρ]
# nonrephasing_se = [(ρ, μp, μm) -> μp * ρ, (ρ, μp, μm) -> ρ * μm, (ρ, μp, μm) -> ρ * μp, (ρ, μp, μm) -> μm * ρ]
# 
# rephasing_esa = [(ρ, μp, μm) -> ρ * μm, (ρ, μp, μm) -> μp * ρ, (ρ, μp, μm) -> μp * ρ, (ρ, μp, μm) -> μm * ρ]
# nonrephasing_esa = [(ρ, μp, μm) -> μp * ρ, (ρ, μp, μm) -> ρ * μm, (ρ, μp, μm) -> μp * ρ, (ρ, μp, μm) -> μm * ρ]

function simulate_absorption(; K::AbstractArray{<:Complex,3}, fbU::AbstractMatrix{<:Complex}, ρ0::AbstractMatrix{<:Complex}, μ::AbstractMatrix{<:Real}, dt::Real, ntimes::Int, L::Union{Nothing,Vector{Matrix{ComplexF64}}}=nothing, ωmin::Union{Nothing,Float64}=nothing, ωmax::Union{Nothing,Float64}=nothing)
    ts, ρs = GQME.propagate_with_memory_kernel(; K, fbU, ρ0=μ * ρ0, dt, ntimes, L)
    corr = [tr(μ * ρs[j, :, :]) for j = 1:ntimes+1]
    ω, spect = Utilities.fourier_transform(ts, corr; full=false, ωmin, ωmax)
    ω, real.(spect)
end

function corr2d(; K::AbstractArray{<:Complex,3}, fbU::AbstractMatrix{<:Complex}, ρ0::AbstractMatrix{<:Complex}, μ::AbstractMatrix{<:Real}, dt::Real, n2::Int, ntimes::Int, L::Union{Nothing,Vector{Matrix{ComplexF64}}}=nothing, funcs)
    # corr[n1, n3]: t1 => n1*dt; t2 => n2*dt; t3 => n3*dt
    @inbounds begin
        μm = UpperTriangular(μ)
        μp = LowerTriangular(μ)

        corr = zeros(ComplexF64, ntimes + 1, ntimes + 1)
        _, ρs = GQME.propagate_with_memory_kernel(; K, fbU, ρ0=funcs[1](ρ0, μp, μm), dt, ntimes, L)
        for n1 = 0:ntimes
            tmpρ0 = deepcopy(ρs[1:n1+1, :, :])
            tmpρ0[end, :, :] .= funcs[2](tmpρ0[end, :, :], μp, μm)
            tmpρs = GQME.propagate_from_middle_with_memory_kernel(; K, fbU, ρ0=tmpρ0, dt, ntimes=n2, L)
            tmpρs[end, :, :] .= funcs[3](tmpρs[end, :, :], μp, μm)
            finalρs = GQME.propagate_from_middle_with_memory_kernel(; K, fbU, ρ0=tmpρs, dt, ntimes, L)
            for n3 = 0:ntimes
                corr[n1+1, n3+1] = tr(funcs[4](finalρs[n1+n2+n3+1, :, :], μp, μm))
            end
        end
    end
    corr
end

function spectrum2D(ts, corr, rephasing, ωmin, ωmax, verbose)
    @inbounds begin
        if verbose
            @info "First set of Fourier transforms"
        end
        ntimes = size(corr, 1)
        flip_sign = rephasing
        w, ft = Utilities.fourier_transform(ts, corr[:, 1]; full=false, ωmin, ωmax, flip_sign)
        ω1 = w
        fourier1 = zeros(ComplexF64, length(ω1), ntimes)
        fourier1[:, 1] = ft
        for j = 2:ntimes
            _, ft = Utilities.fourier_transform(ts, corr[:, j]; full=false, ωmin, ωmax, flip_sign)
            fourier1[:, j] .= ft
        end

        if verbose
            @info "Second set of Fourier transforms"
        end
        w, ft = Utilities.fourier_transform(ts, fourier1[1, :]; full=false, ωmin, ωmax)
        ω3 = w
        fourier2 = zeros(ComplexF64, length(ω1), length(ω3))
        fourier2[1, :] = ft
        for j = 2:length(ω1)
            _, ft = Utilities.fourier_transform(ts, fourier1[j, :]; full=false, ωmin, ωmax)
            fourier2[j, :] .= ft
        end
    end

    ω1, ω3, transpose(fourier2)
end

function simulate_2Dspectrum(; K::AbstractArray{<:Complex,3}, fbU::AbstractMatrix{<:Complex}, ρ0::AbstractMatrix{<:Complex}, μ::AbstractMatrix{<:Real}, dt::Real, n2::Int, ntimes::Int, L::Union{Nothing,Vector{Matrix{ComplexF64}}}=nothing, ωmin::Union{Nothing,Float64}=nothing, ωmax::Union{Nothing,Float64}=nothing, rephasing, non_rephasing, verbose=true)
    ts = 0:dt:ntimes*dt
    if verbose
        @info "Non-Rephasing"
    end
    corr_nr = zeros(ComplexF64, ntimes + 1, ntimes + 1)
    for nonrephas in non_rephasing
        corr_nr .+= corr2d(; K, fbU, ρ0, μ, dt, n2, ntimes, L, funcs=nonrephas)
    end
    ω1, ω3, spect2_nr = spectrum2D(ts, corr_nr, false, ωmin, ωmax, verbose)
    if verbose
        @info "Rephasing"
    end
    corr_r = zeros(ComplexF64, ntimes + 1, ntimes + 1)
    for rephas in rephasing
        corr_r .+= corr2d(; K, fbU, ρ0, μ, dt, n2, ntimes, L, funcs=rephas)
    end
    ω1, ω3, spect2_r = spectrum2D(ts, corr_r, true, ωmin, ωmax, verbose)

    ω1, ω3, real.(spect2_nr + spect2_r)
end


end