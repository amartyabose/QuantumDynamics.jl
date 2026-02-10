module SpinLSC

using HDF5
using ..Utilities
using ..TTM
using ..SolventsX, ..Systems, ..SpectralDensities
using LinearAlgebra: diagm, isdiag, diag
using LinearAlgebra: I as Id

const references = """
- Runeson, J. E.; Richardson, J. O. Spin-mapping approach for non-adiabatic molecular dynamics. J. Chem. Phys. 2019, 151, 044119.
- Runeson, J. E.; Richardson, J. O. Generalized spin mapping for quantum-classical dynamics. J. Chem. Phys. 2020, 152, 084110."""

struct SpinLSCSysPhaseSpace <: SolventsX.PhaseSpace
    X::Vector{Float64}
    P::Vector{Float64}
end

struct SpinLSCSys <: Systems.SpinMappedSystem
    transform::Type{<:Systems.SWTransform}
    h::AbstractMatrix
    ρ₀::Union{Nothing,AbstractMatrix{<:Complex}}
    R²::Float64
    γₛ::Float64
    d::Integer
    bath::SolventsX.HarmonicBathX
    focused_n::Integer
    nsamples::Integer
end
function SpinLSCSys(; transform::Type{<:Systems.SWTransform},
                    Hamiltonian::AbstractMatrix,
                    ρ₀::Union{Nothing,AbstractMatrix{<:Complex}},
                    bath::SolventsX.HarmonicBathX, focused=false,
                    nsamples::Integer)
    @assert nsamples == bath.nsamples
    focused && @assert (isnothing(ρ₀) == false &&
        isdiag(ρ₀) && sum(diag(ρ₀) .== 1) == 1)

    focused_n = focused ? findfirst(diag(ρ₀) .== 1) : -1

    d = size(Hamiltonian, 1)
    SpinLSCSys(transform, Hamiltonian, ρ₀,
               Systems.R²(transform, d), Systems.γ(transform, d),
               d, bath, focused_n, nsamples)
end

# NOTE: For now, we only support ρ₀ = |n⟩⟨n|.
# To support other forms of ρ₀, we need to run multiple trajectories,
# and they are only allowed for W transform.
function focused_sample(sys::SpinLSCSys)
    r = fill(sqrt(sys.γₛ), sys.d)
    r[sys.focused_n] = sqrt(sys.γₛ + 2)

    θ = 2rand(sys.d)

    sinθ, cosθ = sincospi.(θ)

    r .* cosθ, r .* sinθ
end

function Base.iterate(sys::SpinLSCSys, state=1)
    state > sys.nsamples && return nothing

    bathps, _ = iterate(sys.bath, state)
    X, P = sys.focused_n > 0 ? focused_sample(sys) : Systems.sample_XP(sys)

    (SpinLSCSysPhaseSpace(X, P), bathps), state+1
end
Base.eltype(::SpinLSCSys) = SpinLSCSysPhaseSpace
Base.length(s::SpinLSCSys) = s.nsamples
Base.firstindex(::SpinLSCSys) = 1
Base.getindex(s::SpinLSCSys, n::Integer) = iterate(s, n)[1]



"Perform the relevant SW transform of the operator `op`."
Systems.transform_op(sys::SpinLSCSys, op::Union{AbstractMatrix,AbstractVector}, ps::SpinLSCSysPhaseSpace) =
    Systems.transform_op(sys, op, ps.X, ps.P)

const transform_op = Systems.transform_op

"""
    reconstruct_bare_ρ(sys::SpinLSCSys, sps::SpinLSCSysPhaseSpace)

Reconstruct the bare density matrix for the phase space point.

This does NOT multiply by the number of system's dof nor the
Stratonovich–Weyl transform of the initial density matrix.
"""
function reconstruct_bare_ρ(sys::SpinLSCSys, sps::SpinLSCSysPhaseSpace)
    if sys.focused_n < 0
        dual_trans = Systems.dual(sys.transform)
        rescale = Systems.rescale_factor(sys.transform, dual_trans, sys.d)
        γ = Systems.γ(dual_trans, sys.d)
    else
        rescale = 1
        γ = sys.γₛ
    end

    X̄ = rescale * sps.X
    P̄ = rescale * sps.P

    # \[ [\ket{k}\bra{j}]_s = \frac{(X_k - i P_k) (X_j + i P_j) - \delta_kj \gamma_s}{2} \]
    0.5 * ((X̄ + im * P̄) * (X̄ + im * P̄)' - Id(sys.d) * γ)
end

"""
    sampling_func_mult(sys::SpinLSCSys)

Return the sampling function's constant multipler for `sys`.

This returns 1 for focused sampling, and the number of dof of the
system for full-sphere sampling.
"""
sampling_func_mult(sys::SpinLSCSys) = sys.focused_n > 0 ? 1 : sys.d

"""
    sampling_weight(sys::SpinLSCSys, ρ₀::AbstractMatrix{<:Complex}, ps::SpinLSCSysPhaseSpace)

Return the sampling weight for `sys` with initial density matrix `ρ₀`.

This returns `d × [ρ₀]ₛ` for full-sphere sampling, 1 for focused
sampling.
"""
function sampling_weight(sys::SpinLSCSys, ρ₀::AbstractMatrix{<:Complex},
                         ps::SpinLSCSysPhaseSpace)
    (sys.focused_n > 0 ? 1.0 : transform_op(sys, ρ₀, ps)) * sampling_func_mult(sys)
end

function update_dynmap!(U0e::AbstractMatrix{<:Complex},
                        bareρ::AbstractMatrix{<:Complex},
                        sys::SpinLSCSys,
                        sps0::SpinLSCSysPhaseSpace)
    d = size(bareρ, 1)
    d² = size(U0e, 1)

    @inbounds for i in 1:d
        ρ = zeros(ComplexF64, d,d)
        ρ[i,i] = 1.0
        weight = sampling_weight(sys, ρ, sps0)
        U0e[:,i+(i-1)*d] = Utilities.density_matrix_to_vector(bareρ * weight)

        for j in i+1:d
            ρ .= 0.0
            ρ[i,j] = 1.0
            weight = sampling_weight(sys, ρ, sps0)
            n = (i-1)*d + j
            U0e[:,n] = Utilities.density_matrix_to_vector(bareρ * weight)

            # Column corresponding to ρ₀ = ρⱼᵢ.
            m = n + (j-i) * (d-1)

            # Elements corresponding to diagonal terms such as
            # ⟨a|U|i⟩⟨j|U†|a⟩ are simple complex conjugates of each
            # other.
            U0e[1:d+1:d²,m] = conj(U0e[1:d+1:d²,n])

            # Equality corresponding to off-diagonal terms is
            # (⟨a|U|i⟩⟨j|U†|b⟩)* = ⟨b|U|j⟩⟨i|U†|a⟩
            for a = 1:d, b = a+1:d
                # Index for ⟨a|U|i⟩⟨j|U†|b⟩
                l = (b-1)*d + a
                # Index for ⟨b|U|j⟩⟨i|U†|a⟩
                k = (a-1)*d + b
                U0e[k,m] = conj(U0e[l,n])
                U0e[l,m] = conj(U0e[k,n])
            end
        end
    end
end

"Calculate the system force on the baths and store it in `f`."
function Fbath!(sys::SpinLSCSys, ps::SpinLSCSysPhaseSpace, f::Vector{Vector{Float64}})
    @inbounds for b in eachindex(sys.bath.c)
        sₛ = transform_op(sys, sys.bath.s[b], ps)
        @. f[b] = sₛ * sys.bath.c[b]
    end
end

function propagate_trajectory(sys::SpinLSCSys, sps0::SpinLSCSysPhaseSpace,
                              bps0::SolventsX.PhaseSpace, dt::Real, ntimes::Integer,
                              build_dynamical_map::Bool=false)
    XP = [ sps0.X; sps0.P ]
    bps = bps0
    d = sys.d

    U0e = build_dynamical_map && sys.focused_n < 0 ? zeros(ComplexF64, ntimes,sys.d^2,sys.d^2) : nothing
    ρ = isnothing(sys.ρ₀) ? nothing : zeros(ComplexF64, ntimes+1,d,d)

    if !isnothing(ρ)
        w₀ = sampling_weight(sys, sys.ρ₀, sps0)
        ρ[1,:,:] = w₀ * reconstruct_bare_ρ(sys, sps0)

        isnothing(U0e) || (ρ₀ᵥ = Utilities.density_matrix_to_vector(sys.ρ₀))
    end

    dt2 = dt / 2
    bs = sys.bath
    svecs = map(diagm, bs.s)
    LXP = zeros(2d,2d)
    sₛc = similar.(bs.c)
    @inbounds for t in 2:ntimes+1
        sps = SpinLSCSysPhaseSpace(XP[1:d], XP[d+1:2d])
        Fbath!(sys, sps, sₛc)
        bps = SolventsX.propagate_forced_bath(bs, bps, sₛc, dt2, 1)

        LXP[1:d,d+1:2d] = @views sys.h - mapreduce((b, x) -> sum(bs.c[b] .* x) .* svecs[b], +, 1:bs.nbaths, bps.q)
        LXP[d+1:2d,1:d] = -LXP[1:d,d+1:2d]
        XP = exp(LXP * dt) * XP

        sps = SpinLSCSysPhaseSpace(XP[1:d], XP[d+1:2d])
        Fbath!(sys, sps, sₛc)
        bps = SolventsX.propagate_forced_bath(bs, bps, sₛc, dt2, 1)

        bareρ = reconstruct_bare_ρ(sys, sps)
        if !isnothing(U0e)
            update_dynmap!(view(U0e, t-1,:,:), bareρ, sys, sps0)
            if !isnothing(sys.ρ₀)
                ρ[t,:,:] = Utilities.density_matrix_vector_to_matrix(U0e[t-1,:,:] * ρ₀ᵥ)
            end
        else
            ρ[t,:,:] = w₀ * bareρ
        end
    end

    U0e, ρ
end

function propagate_trajectories(sys::SpinLSCSys, dt::Real, ntimes::Integer;
                                output::Union{Nothing,HDF5.Group}=nothing, verbose::Bool=false,
                                build_dynamical_map::Bool=false, kwargs...)
    U0e = build_dynamical_map && sys.focused_n < 0 ? zeros(ComplexF64, ntimes,sys.d^2,sys.d^2) : nothing
    ρ = isnothing(sys.ρ₀) ? nothing : zeros(ComplexF64, ntimes+1,sys.d,sys.d)

    outputρ = if !isnothing(output) && haskey(kwargs, :outgroup)
        Utilities.create_and_select_group(output, kwargs[:outgroup])
    else
        nothing
    end

    mutlock = ReentrantLock()
    ndone = 0
    nthreads = Threads.nthreads()
    stats = @timed Threads.@threads for (sps0, bps0) in sys
        U0eᵢ, ρᵢ = propagate_trajectory(Verlet, sys, sps0, bps0, dt, ntimes, build_dynamical_map)
        lock(mutlock) do
            isnothing(U0e) || (U0e += U0eᵢ)
            isnothing(ρ)   || (ρ += ρᵢ)
            ndone += 1
            verbose && ndone % nthreads == 0 &&
                @info "Trajectories complete: $(100ndone / length(sys))%"
        end
    end
    @info "All trajectories complete\n" *
        "Time taken = $(round(stats.time; digits=3)) sec; memory allocated = $(round(stats.bytes / 1e9; digits=3)) GB; gc time = $(round(stats.gctime; digits=3)) sec"

    if !isnothing(U0e)
        U0e /= length(sys)
        if !isnothing(output)
            output["U0e"] = U0e
            output["T0e"] = TTM.get_Ts(U0e)
            flush(output)
        end
    end

    if !isnothing(ρ)
        ρ /= length(sys)
        if !isnothing(outputρ)
            outputρ["rho"] = ρ
            flush(outputρ)
        end
    end

    U0e, ρ
end



"""
    propagate(; Hamiltonian::Matrix{<:Complex}, Jw::Vector{T},
             β::Real, num_osc::Vector{<:Integer}, svec::Matrix{<:Real},
             ρ0::Union{Nothing,Matrix{<:Complex}}, dt::Real,
             ntimes::Real, transform::Type{<:Systems.SWTransform},
             nmc::Integer, focused::Bool=false, verbose::Bool=false,
             build_dynamical_map::Bool=false, kwargs...) where {T<:SpectralDensities.SpectralDensity}

Propagate the system using the spin-mapped LSC method.

Arguments:
- `ρ0`: the initial density matrix.  If it is `nothing`, then build
  only the dynamical map
- `Hamiltonian`: the Hamiltonian of the sub-system
- `Jw`: list of spectral densities
- `β`: the inverse temperature of the baths
- `num_osc`: the number of oscillators for each bath
- `svec`: diagonal elements of system operators through which the
  corresponding baths interact
- `transform`: the Stratonovich–Weyl transform to use for the
  Hamiltonian
- `dt`: the time step for the propagation
- `focused`: should focused initial sampling be used
- `nmc`: the number of Monte-Carlo samples

Propagate the density matrix and build the dynamical map by doing
linearised semiclassical propagator using the given Stratonovich–Weyl
transform for the Hamiltonian of the systemnnnnnnnnnnnn.
"""
function propagate(; Hamiltonian::Matrix{<:Complex}, Jw::Vector{T},
                   β::Real, num_osc::Vector{<:Integer}, svec::Matrix{<:Real},
                   ρ0::Union{Nothing,Matrix{<:Complex}}, dt::Real,
                   ntimes::Real, transform::Type{<:Systems.SWTransform},
                   nmc::Integer, focused::Bool=false, verbose::Bool=false,
                   build_dynamical_map::Bool=false, kwargs...) where {T<:SpectralDensities.SpectralDensity}
    nbaths = length(Jw)
    c = Vector{Vector{Float64}}(undef, nbaths)
    ω = Vector{Vector{Float64}}(undef, nbaths)
    s = Vector{Vector{Float64}}(undef, nbaths)

    for n in 1:nbaths
        ω[n], c[n] = SpectralDensities.discretize(Jw[n], num_osc[n])
        s[n] = svec[n,:]
    end

    bath = SolventsX.HarmonicBathX(; β, ω, c, svecs=s, nsamples=nmc)
    sys = SpinLSCSys(; transform, Hamiltonian, ρ₀=ρ0, bath, focused, nsamples=nmc)

    propagate_trajectories(sys, dt, ntimes; verbose, build_dynamical_map, kwargs...)
end

end
