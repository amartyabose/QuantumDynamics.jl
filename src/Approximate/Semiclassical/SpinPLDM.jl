module SpinPLDM

using HDF5
using ..Utilities
using ..SolventsX, ..Systems, ..SpectralDensities
using LinearAlgebra: diagm

const references = """
- Mannouch, J. R.; Richardsion, J. O. A partially linearised spin-mapping approach for non-adiabatic dynamics. I. Derivation of the theory. J. Chem. Phys. 2020 153, 194109."""

struct SpinPLDMSysPhaseSpace <: SolventsX.PhaseSpace
    Xf::Vector{Float64}
    Pf::Vector{Float64}
    Xb::Vector{Float64}
    Pb::Vector{Float64}
end

struct SpinPLDMSys <: Systems.SpinMappedSystem
    transform::Type{<:Systems.SWTransform}
    h::AbstractMatrix
    ρ₀::AbstractMatrix{<:Complex}
    R²::Float64
    γₛ::Float64
    d::Integer
    bath::SolventsX.HarmonicBathX
    nsamples::Integer
end
function SpinPLDMSys(; transform::Type{<:Systems.SWTransform},
                     Hamiltonian::AbstractMatrix,
                     ρ₀::AbstractMatrix{<:Complex},
                     bath::SolventsX.HarmonicBathX, nsamples::Integer)
    @assert nsamples == bath.nsamples
    d = size(Hamiltonian, 1)
    SpinPLDMSys(transform, Hamiltonian, ρ₀,
                Systems.R²(transform, d), Systems.γ(transform, d),
                d, bath, nsamples)
end

function Base.iterate(sys::SpinPLDMSys, state=1)
    state > sys.nsamples && return nothing

    bathps, _ = iterate(sys.bath, state)

    Xf, Pf = Systems.sample_XP(sys)
    Xb, Pb = Systems.sample_XP(sys)

    (SpinPLDMSysPhaseSpace(Xf, Pf, Xb, Pb), bathps), state+1
end
Base.eltype(::SpinPLDMSys) = SpinPLDMSysPhaseSpace
Base.length(s::SpinPLDMSys) = s.nsamples
Base.firstindex(::SpinPLDMSys) = 1
Base.getindex(s::SpinPLDMSys, n::Integer) = iterate(s, n)[1]



transform_op_fwd(sys::SpinPLDMSys, op::Union{AbstractVector,AbstractMatrix}, ps::SpinPLDMSysPhaseSpace) =
    Systems.transform_op(sys, op, ps.Xf, ps.Pf)

transform_op_bwd(sys::SpinPLDMSys, op::Union{AbstractVector,AbstractMatrix}, ps::SpinPLDMSysPhaseSpace) =
    Systems.transform_op(sys, op, ps.Xb, ps.Pb)

function transform_kernel(sys::SpinPLDMSys,
                          X₀::Vector{<:Real}, P₀::Vector{<:Real},
                          Xₜ::Vector{<:Real}, Pₜ::Vector{<:Real},
                          U::AbstractMatrix{<:Complex})
    dual = Systems.dual(sys.transform)
    rescale = Systems.rescale_factor(sys.transform, dual, sys.d)
    γs̄ = Systems.γ(dual, sys.d)

    X̄₀ = X₀ * rescale
    P̄₀ = P₀ * rescale
    X̄ₜ = Xₜ * rescale
    P̄ₜ = Pₜ * rescale

    0.5 * ((X̄ₜ + im * P̄ₜ) * (X̄₀ + im * P̄₀)' - γs̄ * U)
end

function build_ρ!(sys::SpinPLDMSys, sps0::SpinPLDMSysPhaseSpace,
                  XPf::Vector{Float64}, XPb::Vector{Float64},
                  ρ::AbstractMatrix{<:Complex}, U::AbstractMatrix{<:Complex})
    wf = transform_kernel(sys, sps0.Xf, sps0.Pf, XPf[1:d], Xpf[d+1:2d], U)
    wb = transform_kernel(sys, sps0.Xb, sps0.Pb, XPb[1:d], Xpb[d+1:2d], U)
    ρ = sys.d^2 * (wf * sys.ρ₀ * wb' + wb * sys.ρ₀ * wf') / 2
end

function propagate_trajectory(sys::SpinPLDMSys,
                              sps0::SpinPLDMSysPhaseSpace,
                              bps0::SolventsX.PhaseSpace,
                              dt::Real, ntimes::Integer)
    XPf = [ sps0.Xf; sps0.Pf ]
    XPb = [ sps0.Xb; sps0.Pb ]
    x = bps0.q
    p = bps0.p
    d = sys.d

    ρ = zeros(ComplexF64, ntimes+1,d,d)
    U = diagm(ones(ComplexF64, sys.d))

    @views build_ρ!(sys, sps0, XPf, Xpb, ρ[1,:,:], U)

    δtₓ = dt / 100
    N½ = 50
    mω² = map(b -> -sys.bath.ω[b].^2, 1:sys.bath.nbaths)
    bs = sys.bath
    svecs = map(diagm, bs.s)
    LXP = zeros(2d,2d)
    @inbounds for t in 2:ntimes+1
        sps = SpinPLDMSysPhaseSpace(XPf[1:d], XPf[d+1:2d], XPb[1:d], XPb[d+1:2d])
        for b in 1:bs.nbaths
            s̄ₛc = bs.c[b] * (transform_op_fwd(sys, bs.s[b], sps) +
                             transform_op_bwd(sys, bs.s[b], sps)) / 2
            for _ in 1:N½
                @. p[b] = p[b] + 0.5 * (mω²[b] * x[b] + s̄ₛc) * δtₓ
                @. x[b] = x[b] + p[b] * δtₓ
                @. p[b] = p[b] + 0.5 * (mω²[b] * x[b] + s̄ₛc) * δtₓ
            end
        end

        LXP[1:d,d+1:2d] = @views sys.h - mapreduce((b, x) -> sum(bs.c[b] .* x) * svecs[b], +, 1:bs.nbaths, x)
        LXP[d+1:2d,1:d] = -LXP[1:d,d+1:2d]
        eLXP = exp(LXP * dt)
        XPf = eLXP * XPf
        XPb = eLXP * XPb
        U = exp(-im * V * dt) * U

        sps = SpinPLDMSysPhaseSpace(XPf[1:d], XPf[d+1:2d], XPb[1:d], XPb[d+1:2d])
        for b in 1:bs.nbaths
            s̄ₛc = bs.c[b] * (transform_op_fwd(sys, bs.s[b], sps) +
                             transform_op_bwd(sys, bs.s[b], sps)) / 2
            for _ in 1:N½
                @. p[b] = p[b] + 0.5 * (mω²[b] * x[b] + s̄ₛc) * δtₓ
                @. x[b] = x[b] + p[b] * δtₓ
                @. p[b] = p[b] + 0.5 * (mω²[b] * x[b] + s̄ₛc) * δtₓ
            end
        end

        @views build_ρ!(sys, sps0, XPf, Xpb, ρ[1,:,:], U)
    end

    ρ
end

function propagate_trajectories(sys::SpinPLDMSys, dt::Real, ntimes::Integer;
                                output::Union{Nothing,HDF5.Group}=nothing,
                                verbose::Bool=false, kwargs...)
    ρ = isnothing(sys.ρ₀) ? nothing : zeros(ComplexF64, ntimes+1,sys.d,sys.d)
    outputρ = if !isnothing(ρ) && haskey(kwargs, :outgroup)
        Utilities.create_and_select_group(output, kwargs[:outgroup])
    else
        nothing
    end

    mutlock = ReentrantLock()
    ndone = 0
    nthreads = Threads.nthreads()
    stats = @timed Threads.@threads for (sps0, bps0) in sys
        ρᵢ = propagate_trajectory(sys, sps0, bps0, dt, ntimes)
        lock(mutlock) do
            ndone += 1
            isnothing(ρ) || (ρ += ρᵢ)
            verbose && ndone % nthreads == 0 &&
                @info "Trajectories complete: $(100ndone / length(sys))%"
        end
    end
    @info "All trajectories complete"
    @info "Time taken = $(round(stats.time; digits=3)) sec; memory allocated = $(round(stats.bytes / 1e9; digits=3)) GB; gc time = $(round(stats.gctime; digits=3)) sec"

    if !isnothing(ρ)
        ρ /= length(sys)
        if !isnothing(outputρ)
            outputρ["rho"] = ρ
            flush(outputρ)
        end
    end

    ρ
end

"""
    propagate(; Hamiltonian::Matrix{<:Complex}, Jw::Vector{T},
              β::Real, num_osc::Vector{<:Integer}, svec::Matrix{<:Real},
              ρ0::Matrix{<:Complex}, dt::Real,
              ntimes::Real, transform::Type{<:Systems.SWTransform},
              nmc::Integer, verbose::Bool=false,
              kwargs...) where {T<:SpectralDensities.SpectralDensity}

Propagate the system using the spin-mapped PLDM method.

Arguments:
- `ρ0`: initial reduced density matrix
- `Hamiltonian`: the Hamiltonian of the sub-system
- `Jw`: list of spectral densities
- `β`: the inverse temperature of the bath
- `num_osc`: the number of oscillator for each bath
- `svec`: diagonal elements of system operators through which the
  corresponding baths interact
- `transform`: the Stratonovich–Weyl transformation to use for the
  Hamiltonian
- `dt`: the time step for the propagation
- `nmc`: the number of Monte-Carlo samples

Propagate the density matrix using a partially linearised propagator
for the system-bath problem, using the given Stratonovich–Weyl
transform for the Hamiltonian.
"""
function propagate(; Hamiltonian::Matrix{<:Complex}, Jw::Vector{T},
                   β::Real, num_osc::Vector{<:Integer}, svec::Matrix{<:Real},
                   ρ0::Matrix{<:Complex}, dt::Real,
                   ntimes::Real, transform::Type{<:Systems.SWTransform},
                   nmc::Integer, verbose::Bool=false,
                   kwargs...) where {T<:SpectralDensities.SpectralDensity}
    nbaths = length(Jw)
    c = Vector{Vector{Float64}}(undef, nbaths)
    ω = Vector{Vector{Float64}}(undef, nbaths)
    s = Vector{Vector{Float64}}(undef, nbaths)

    for n in 1:nbaths
        ω[n], c[n] = SpectralDensities.discretize(Jw[n], num_osc[n])
        s[n] = svec[n,:]
    end

    bath = SolventsX.HarmonicBathX(; β, ω, c, svecs=s, nsamples=nmc)
    sys = SpinPLDMSys(; transform, Hamiltonian, ρ₀=ρ0, bath, nsamples=nmc)

    propagate_trajectories(sys, dt, ntimes; verbose, kwargs...)
end

end
