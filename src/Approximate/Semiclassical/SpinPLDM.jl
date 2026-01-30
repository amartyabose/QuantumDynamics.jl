module SpinPLDM

using HDF5
using ..Utilities
using ..SolventsX, ..Systems, ..SpectralDensities
using LinearAlgebra: diagm

const references = """
- Mannouch, J. R.; Richardsion, J. O. A partially linearised spin-mapping approach for non-adiabatic dynamics. I. Derivation of the theory. J. Chem. Phys. 2020 153, 194109."""

struct SpinPLDMSysPhaseSpace <: SolventsX.PhaseSpace
    Xf::Vector{<:Real}
    Pf::Vector{<:Real}
    Xb::Vector{<:Real}
    Pb::Vector{<:Real}
end

struct SpinPLDMSys <: Systems.SpinMappedSystem
    transform::Type{<:Systems.SWTransform}
    h::AbstractMatrix
    ρ₀::Union{Nothing,AbstractMatrix{<:Complex}}
    R²::Real
    γₛ::Real
    d::Integer
    bath::SolventsX.Solvent
    nsamples::Integer
end
function SpinPLDMSys(; transform::Type{<:Systems.SWTransform},
                     Hamiltonian::AbstractMatrix,
                     ρ₀::Union{Nothing,AbstractMatrix{<:Complex}},
                     bath::SolventsX.Solvent, nsamples::Integer)
    @assert nsamples == bath.nsamples
    d = size(Hamiltonian, 1)
    SpinPLDMSys(transform, Hamiltonian, ρ₀,
                Systems.R²(transform, d), Systems.γ(transform, d),
                d, bath, nsamples)
end


function Base.iterate(sys::SpinPLDMSys, state=1)
    state > sys.nsamples && return nothing

    bathps, _ = iterate(sys.bath, state)

    R = sqrt(sys.R²)
    Xf = randn(sys.d)
    Pf = randn(sys.d)
    Xb = randn(sys.d)
    Pb = randn(sys.d)

    sqΣf = sqrt(sum(Xf.^2 + Pf.^2))
    sqΣb = sqrt(sum(Xb.^2 + Pb.^2))

    Xf *= R/sqΣf
    Pf *= R/sqΣf
    Xb *= R/sqΣb
    Pb *= R/sqΣb

    (SpinPLDMSysPhaseSpace(Xf, Pf, Xb, Pb), bathps), state+1
end
Base.eltype(::SpinPLDMSys) = SpinPLDMSysPhaseSpace
Base.length(s::SpinPLDMSys) = s.nsamples
Base.firstindex(::SpinPLDMSys) = 1
Base.getindex(s::SpinPLDMSys, n::Integer) = iterate(s, n)[1]



transform_op_fwd(sys::SpinPLDMSys, op::AbstractMatrix, ps::SpinPLDMSysPhaseSpace) =
    Systems.transform_op(sys, op, ps.Xf, ps.Pf)

transform_op_bwd(sys::SpinPLDMSys, op::AbstractMatrix, ps::SpinPLDMSysPhaseSpace) =
    Systems.transform_op(sys, op, ps.Xb, ps.Pb)

function Fbath(sys::SpinPLDMSys, sps::SpinPLDMSysPhaseSpace,
               q::AbstractVector{<:Real}, i::Integer)
    s = diagm(sys.bath.s[i])
    -sys.bath.ω[i].^2 .* q .+
        sys.bath.c[i] *
          (transform_op_fwd(sys, s, sps) +
           transform_op_bwd(sys, s, sps)) / 2
end

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

function propagate_trajectory(sys::SpinPLDMSys,
                              sps0::SpinPLDMSysPhaseSpace,
                              bps0::SolventsX.PhaseSpace,
                              dt::Real, ntimes::Integer)
    XPf = [ sps0.Xf; sps0.Pf ]
    XPb = [ sps0.Xb; sps0.Pb ]
    x = bps0.q
    p = bps0.p
    d = sys.d
    d² = d^2

    ρ = isnothing(sys.ρ₀) ? nothing : zeros(ComplexF64, ntimes+1,d,d)
    U = diagm(ones(ComplexF64, sys.d))

    build_ρ!(t) = if !isnothing(ρ)
        wf = transform_kernel(sys, sps0.Xf, sps0.Pf, XPf[1:d], XPf[d+1:2d], U)
        wb = transform_kernel(sys, sps0.Xb, sps0.Pb, XPb[1:d], XPb[d+1:2d], U)
        ρ[t,:,:] = d² * (wf * sys.ρ₀ * wb' + wb * sys.ρ₀ * wf') / 2
    end

    build_ρ!(1)

    δtₓ = dt / 100
    N½ = dt / 2 / δtₓ
    propagate_xp!() = let sps = SpinPLDMSysPhaseSpace(XPf[1:d], XPf[d+1:2d],
                                                      XPb[1:d], XPb[d+1:2d])
        for b in 1:sys.bath.nbaths
            for _ in 1:N½
                p[b] = p[b] .+ 0.5 * Fbath(sys, sps, x[b], b) * δtₓ
                x[b] = x[b] .+ p[b] * δtₓ
                p[b] = p[b] .+ 0.5 * Fbath(sys, sps, x[b], b) * δtₓ
            end
        end
    end

    bs = sys.bath
    svecs = map(diagm, bs.s)
    LXP = zeros(2d,2d)
    for t in 2:ntimes+1
        propagate_xp!()

        V = sys.h - mapreduce((b, x) -> sum(bs.c[b] .* x) * svecs[b], +,
                              1:bs.nbaths, x)
        LXP[1:d,d+1:2d] = V
        LXP[d+1:2d,1:d] = -V
        eLXP = exp(LXP * dt)
        XPf = eLXP * XPf
        XPb = eLXP * XPb
        U = exp(-im * V * dt) * U

        propagate_xp!()

        build_ρ!(t)
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
    Threads.@threads for (sps0, bps0) in sys
        ρᵢ = propagate_trajectory(sys, sps0, bps0, dt, ntimes)
        lock(mutlock) do
            ndone += 1
            isnothing(ρ) || (ρ += ρᵢ)
            verbose && ndone % nthreads == 0 &&
                @info "Trajectories complete: $(100ndone / length(sys))"
        end
    end

    if !isnothing(ρ)
        ρ /= length(sys)
        if !isnothing(outputρ)
            outputρ["rho"] = ρ
            flush(outputρ)
        end
    end

    ρ
end

function propagate(; Hamiltonian::Matrix{<:Complex}, Jw::Vector{T},
                   β::Real, num_osc::Vector{<:Integer}, svec::Matrix{<:Real},
                   ρ0::Union{Nothing,Matrix{<:Complex}}, dt::Real,
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
