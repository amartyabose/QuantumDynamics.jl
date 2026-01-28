module SpinLSC

using HDF5
using ..Utilities
using ..TTM
using ..SolventsX, ..Systems, ..SpectralDensities
using LinearAlgebra: diagm
import OrdinaryDiffEq as ODE

const references = """
- Runeson, J. E.; Richardson, J. O. Spin-mapping approach for non-adiabatic molecular dynamics. J. Chem. Phys. 2019, 151, 044119.
- Runeson, J. E.; Richardson, J. O. Generalized spin mapping for quantum-classical dynamics. J. Chem. Phys. 2020, 152, 084110."""

struct SpinLSCSysPhaseSpace <: SolventsX.PhaseSpace
    X::Vector{<:Real}
    P::Vector{<:Real}
end

struct SpinLSCSys <: System.SpinMappedSystem
    data::System.SpinMappedSystem
end
function SpinLSCSys(; transform::Type{<:Systems.SWTransform},
                    Hamiltonian::AbstractMatrix,
                    ρ₀::Union{Nothing,AbstractMatrix{<:Complex}},
                    bath::Solvent, nsamples::Integer)
    SpinLSCSys(
        System.SpinMappedSystem(; transform, Hamiltonian, ρ₀,
                                bath, nsamples))
end
Base.getproperty(sys::SpinLSCSys, s::Symbol) = getproperty(sys.data, s)

function Base.iterate(sys::SpinLSCSys, state=1)
    state > sys.nsamples && return nothing

    bathps, _ = iterate(sys.bath)

    R = sqrt(sys.R²)
    X = randn(sys.d)
    P = randn(sys.d)

    sqΣ = sqrt(sum(X.^2 .+ P.^2))

    (SpinLSCSysPhaseSpace(X * R / sqΣ, P * R / sqΣ), bathps), state+1
end
Base.eltype(::SpinLSCSys) = SpinLSCSysPhaseSpace
Base.length(s::SpinLSCSys) = s.nsamples



"Perform the relevant SW transform of the operator `op`."
transform_op(sys::SpinLSCSys, op::AbstractMatrix, ps::SpinLSCSysPhaseSpace) =
    Systems.transform_op(sys, op, ps.X, ps.P)

"Calculate the force on the `i`th bath."
function Fbath(sys::SpinLSCSys, sps::SpinLSCSysPhaseSpace,
               q::AbstractVector{<:Real}, i::Integer)
    -sys.bath.ω[i].^2 .* q .+
        transform_op(sys, diagm(sys.bath.s[i]), sps) * sys.bath.c[i]
end

function H(sys::SpinLSCSys, sps::SpinLSCSysPhaseSpace, bps::SolventsX.PhaseSpace)
    bs = sys.bath
    transform_op(sys, sys.h, sps) +
        mapreduce((b, x, p) -> 0.5 * p.^2 + 0.5 * bs.ω[b].^2 * x.^2, +,
                  1:sys.bath.nbaths, bps.q, bps.p) +
       -mapreduce((b,x) -> transform_op(sys, sum(bs.c[b] .* x) * diagm(bs.s[b]), sps), +,
                  1:sys.bath.nbaths, bps.q)
end

"""Reconstruct the bare density matrix for the phase space point.

This does NOT multiply the Stratonovich–Weyl transform of the initial
density matrix.  Intended to be used when constructing the dynamical
maps.
"""
function reconstruct_bare_ρ(sys::SpinLSCSys, sps::SpinLSCSysPhaseSpace)
    ρ = zeros(ComplexF64, sys.d,sys.d)

    dual_trans = Systems.dual(sys.transform)
    rescale = Systems.rescale_factor(sys.transform, dual_trans, sys.d)
    γs̄ = Systems.γ(dual_trans, sys.d)

    X̄ = rescale * sps.X
    P̄ = rescale * sps.P

    for j = 1:sys.d
        # Some formulae:
        # \[ [\ket{k}\bra{j}]_s = \frac{(X_k - i P_k) (X_j + i P_j)}{2} \]
        # \[ [\ket{j}\bra{j}]_s = \frac{X_j^2 + P_j^2 - \gamma_s}{2} \]
        # with the appropriate rescaling of X_j's.

        ρ[j,j] = 0.5 * (X̄[j].^2 + P̄[j].^2 - γs̄)
        for k = j+1:sys.d
            ρ[j,k] = 0.5 * (X̄[k] - im * P̄[k]) * (X̄[j] + im * P̄[j])
            ρ[k,j] = 0.5 * (X̄[j] - im * P̄[j]) * (X̄[k] + im * P̄[k])
        end
    end

    ρ * sys.d
end



abstract type SpinLSCSolver end
abstract type RK4 <: SpinLSCSolver end

function HX(sys::SpinLSCSys, sps::SpinLSCSysPhaseSpace, bps::SolventsX.PhaseSpace)
    ddX = zeros(sys.d)

    bs = sys.bath
    for i in 1:sys.d
        # NOTE: The diagonal term of the Hamiltonian enters naturally
        # from the summation.
        ddX[i] = 0.5 * sum(sys.h[i,:] .* (sps.X .+ im * sps.P)) +
            0.5 * sum(sys.h[:,i] .* (sps.X .- im * sps.P)) +
            -mapreduce((b,x) -> sum(bs.c[b] .* x) * bs.s[b][i] * sps.X[i], +,
                       1:sys.bath.nbaths, bps.q)
    end

    ddX
end

function HP(sys::SpinLSCSys, sps::SpinLSCSysPhaseSpace, bps::SolventsX.PhaseSpace)
    ddP = zeros(sys.d)

    bs = sys.bath
    for i in 1:sys.d
        # NOTE: The diagonal term of the Hamiltonian enters naturally
        # from the summation.
        ddP[i] = 0.5 * sum(sys.h[:,i] .* (sps.P .+ im * sps.X)) +
            0.5 * sum(sys.h[i,:] .* (sps.P .- im * sps.X)) +
            -mapreduce((b, x) -> sum(bs.c[b] .* x) * bs.s[b][i] * sps.P[i], +,
                       1:sys.bath.nbaths, bps.q)
    end

    ddP
end

function propagate_xpXP!(du, u, p, t)
    sys, xis, pis = p
    d = sys.d

    sps = SpinLSCSysPhaseSpace(u[1:d], u[d+1:2d])
    bps = SolventsX.HarmonicPhaseSpaceX([ u[i] for i in xis ], [ u[i] for i in pis ])

    du[1:d] = HP(sys, sps, bps)     # Ẋ
    du[d+1:2d] = -HX(sys, sps, bps) # Ṗ

    for n in 1:sys.bath.nbaths
        du[xis[n]] = bps.p[n]   # ẋ
        du[pis[n]] = Fbath(sys, sps, bps.q[n], n)
    end
end

function bathinds(sys::SpinLSCSys)
    m = length.(sys.bath.ω)
    total = sum(m)
    xis = Vector{UnitRange{Int}}(undef, sys.bath.nbaths)
    pis = Vector{UnitRange{Int}}(undef, sys.bath.nbaths)

    prevstart = 2sys.d
    for (n,len) in enumerate(m)
        xis[n] = prevstart+1:prevstart+len
        pis[n] = total+prevstart+1:total+prevstart+len
        prevstart = prevstart+len
    end

    xis, pis
end

function update_dynmap!(U0e::AbstractMatrix{<:Complex},
                        bareρ::AbstractMatrix{<:Complex},
                        sys::SpinLSCSys,
                        sps0::SpinLSCSysPhaseSpace)
    d = size(bareρ, 1)
    d² = size(U0e, 1)

    for i in 1:d
        ρ = zeros(ComplexF64, d,d)
        ρ[i,i] = 1.0
        ρₛ = transform_op(sys, ρ, sps0)
        U0e[:,i+(i-1)*d] = Utilities.density_matrix_to_vector(bareρ * ρₛ)

        for j in i+1:d
            ρ .= 0.0
            ρ[i,j] = 1.0
            ρₛ = transform_op(sys, ρ, sps0)
            n = (i-1)*d + j
            U0e[:,n] = Utilities.density_matrix_to_vector(bareρ * ρₛ)

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

function build_dynmap_ρ(sol::ODE.ODESolution)
    sys, _ = sol.prob.p
    Nₜ = length(sol.u)
    d = sys.d
    sps0 = SpinLSCSysPhaseSpace(sol.u[1][1:d], sol.u[1][d+1:2d])

    U0e = zeros(ComplexF64, Nₜ-1,sys.d^2,sys.d^2)
    if !isnothing(sys.ρ₀)
        ρ = zeros(ComplexF64, Nₜ,d,d)
        ρ[1,:,:] = transform_op(sys, sys.ρ₀, sps0) * reconstruct_bare_ρ(sys, sps0)
        ρ₀ᵥ = Utilities.density_matrix_to_vector(sys.ρ₀)
    end

    for t in 2:Nₜ
        sps = SpinLSCSysPhaseSpace(sol.u[t][1:d], sol.u[t][d+1:2d])
        bareρ = reconstruct_bare_ρ(sys, sps)
        update_dynmap!(view(U0e, t-1,:,:), bareρ, sys, sps0)

        if !isnothing(sys.ρ₀)
            ρ[t,:,:] = Utilities.density_matrix_vector_to_matrix(U0e[t-1,:,:] * ρ₀ᵥ)
        end
    end

    (U0e, isnothing(sys.ρ₀) ? nothing : ρ)
end

function propagate_trajectories(::Type{RK4}, sys::SpinLSCSys, dt::Real, ntimes::Integer;
                                output::Union{Nothing,HDF5.Group}=nothing, verbose::Bool=false, kwargs...)
    xis, pis = bathinds(sys)

    outputρ = if !isnothing(output) && haskey(kwargs, :outgroup)
        Utilities.create_and_select_group(output, kwargs[:outgroup])
    else
        nothing
    end

    probfn(p, i, _) = begin
        ps, _ = iterate(sys, i)
        sps, bps = ps
        ODE.remake(p, u0=vcat(sps.X, sps.P, bps.q..., bps.p...))
    end
    outputfn(sol, _) = (build_dynmap_ρ(sol), false)

    done = 0
    reducefn(data, us, I) = begin
        done += length(I)
        verbose && @info "Trajectories completed: $(done * 100 / length(sys))%"

        U0e = data[1] + sum(getindex.(us,1))
        ρ = isnothing(data[2]) ? nothing : data[2] + sum(getindex.(us, 2))

        ((U0e, ρ), false)
    end

    ensemble = ODE.EnsembleProblem(
        ODE.ODEProblem(propagate_xpXP!,
                       zeros(2sys.d+2sum(length.(sys.bath.ω))),
                       (0.0, ntimes*dt),
                       (sys, xis, pis));
        output_func=outputfn,
        prob_func=probfn,
        reduction=reducefn,
        u_init=(zeros(ComplexF64, ntimes,sys.d^2,sys.d^2),
                isnothing(sys.ρ₀) ? nothing : zeros(ComplexF64, ntimes+1,sys.d,sys.d)))

    sol = ODE.solve(ensemble, ODE.RK4(), ODE.EnsembleThreads();
                    dt, saveat=dt, trajectories=length(sys),
                    batch_size=Threads.nthreads())

    U0e = sol.u[1] / length(sys)
    ρ = isnothing(sys.ρ₀) ? nothing : sol.u[2] / length(sys)

    if !isnothing(output)
        output["U0e"] = U0e
        output["T0e"] = TTM.get_Ts(U0e)
        flush(output)
    end

    if !isnothing(sys.ρ₀) && !isnothing(outputρ)
            outputρ["rho"] = ρ
            flush(outputρ)
    end

    U0e, ρ
end



abstract type Verlet <: SpinLSCSolver end

function propagate_trajectory(::Type{Verlet}, sys::SpinLSCSys,
                              sps0::SpinLSCSysPhaseSpace,
                              bps0::SolventsX.PhaseSpace,
                              dt::Real, ntimes::Integer)
    XP = [ sps0.X; sps0.P ]
    x = bps0.q
    p = bps0.p
    d = sys.d

    U0e = zeros(ComplexF64, ntimes,sys.d^2,sys.d^2)
    isnothing(sys.ρ₀) || (ρ₀ᵥ = Utilities.density_matrix_to_vector(sys.ρ₀))
    if !isnothing(sys.ρ₀)
        ρ = zeros(ComplexF64, ntimes+1,d,d)
        ρ[1,:,:] = transform_op(sys, sys.ρ₀, sps0) * reconstruct_bare_ρ(sys, sps0)
    end

    build_dynmap_ρ!(t) = begin
        sps = SpinLSCSysPhaseSpace(XP[1:d], XP[d+1:2d])
        bareρ = reconstruct_bare_ρ(sys, sps)
        update_dynmap!(view(U0e, t-1,:,:), bareρ, sys, sps0)

        if !isnothing(sys.ρ₀)
            ρ[t,:,:] = Utilities.density_matrix_vector_to_matrix(U0e[t-1,:,:] * ρ₀ᵥ)
        end
    end

    δtₓ = dt / 100
    N½ = dt / 2 / δtₓ
    propagate_xp!(t) = for b in 1:sys.bath.nbaths
        for _ in 1:N½
            p[b] = p[b] .+ 0.5 * Fbath(sys, SpinLSCSysPhaseSpace(XP[1:d], XP[d+1:2d]), x[b], b) * δtₓ
            x[b] = x[b] .+ p[b] * δtₓ
            p[b] = p[b] .+ 0.5 * Fbath(sys, SpinLSCSysPhaseSpace(XP[1:d], XP[d+1:2d]), x[b], b) * δtₓ
        end
    end

    bs = sys.bath
    LXP = zeros(2d,2d)
    for t in 2:ntimes+1
        propagate_xp!(t)

        V = sys.h - mapreduce((b, x) -> sum(bs.c[b] .* x) * diagm(bs.s[b]), +,
                              1:bs.nbaths, x)
        LXP[1:d,d+1:2d] = V
        LXP[d+1:2d,1:d] = -V
        XP = exp(LXP * dt) * XP

        propagate_xp!(t)

        build_dynmap_ρ!(t)
    end

    U0e, isnothing(sys.ρ₀) ? nothing : ρ
end

function propagate_trajectories(::Type{Verlet}, sys::SpinLSCSys, dt::Real, ntimes::Integer;
                                output::Union{Nothing,HDF5.Group}=nothing, verbose::Bool=false, kwargs...)
    U0e = zeros(ComplexF64, ntimes,sys.d^2,sys.d^2)
    isnothing(sys.ρ₀) || (ρ = zeros(ComplexF64, ntimes+1,sys.d,sys.d))

    outputρ = if !isnothing(output) && haskey(kwargs, :outgroup)
        Utilities.create_and_select_group(output, kwargs[:outgroup])
    else
        nothing
    end

    batches = Iterators.partition(1:length(sys), Threads.nthreads())
    for samples in batches
        tasks = map(samples) do state
            ps, _ = iterate(sys, state)
            sps0, bps0 = ps
            Threads.@spawn(propagate_trajectory(Verlet, sys, sps0, bps0, dt, ntimes))
        end
        solns = fetch.(tasks)

        U0e += sum(getindex.(solns, 1))
        if !isnothing(sys.ρ₀)
            ρ += sum(getindex.(solns, 2))
        end
        verbose && @info "Trajectories complete: $(samples[end] * 100 / length(sys))%"
    end

    U0e /= length(sys)

    if !isnothing(output)
        output["U0e"] = U0e
        output["T0e"] = TTM.get_Ts(U0e)
        flush(output)
    end

    if !isnothing(sys.ρ₀)
        ρ /= length(sys)
        if !isnothing(outputρ)
            outputρ["rho"] = ρ
            flush(outputρ)
        end
    end

    U0e, isnothing(sys.ρ₀) ? nothing : ρ
end



"""
    propagate(; Hamiltonian::Matrix{<:Complex}, Jw::Vector{T},
             β::Real, num_bath_modes::Vector{<:Integer}, svec::Matrix{<:Real},
             ρ0::Union{Nothing,Matrix{<:Complex}}, dt::Real,
             ntimes::Real, transform::Type{<:Systems.SWTransform},
             nmc::Integer, solver::Type{<:SpinLSCSolver}, verbose::Bool=false,
             kwargs...) where {T<:SpectralDensities.SpectralDensity}

Propagate the system using the spin-mapped LSC method.

Arguments:
- `ρ0`: initial reduced density matrix
- `Hamiltonian`: the Hamiltonian of the sub-system
- `Jw`: list of spectral densities
- `β`: the inverse temperature of the baths
- `num_bath_modes`: a list of discretisation points for each bath
- `svec`: diagonal elements of system operators through which the
          corresponding baths interact
- `ρ0`: the initial density matrix.  If it is `nothing`, then build
        only the dynamical map.
- `dt`: the time step for the propagation
- `transform`: the Stratonovich–Weyl transform to use to calculate the
  correlation functions
- `solver`: the algorithm to use to solve the EOMs
- `nmc`: the number of Monte-Carlo samples

Propagate the density matrix and build the dynamical map by doing
linearised semiclassical propagator using the given Stratonovich–Weyl
transform.
"""
function propagate(; Hamiltonian::Matrix{<:Complex}, Jw::Vector{T},
                   β::Real, num_bath_modes::Vector{<:Integer}, svec::Matrix{<:Real},
                   ρ0::Union{Nothing,Matrix{<:Complex}}, dt::Real,
                   ntimes::Real, transform::Type{<:Systems.SWTransform},
                   nmc::Integer, solver::Type{<:SpinLSCSolver}, verbose::Bool=false,
                   kwargs...) where {T<:SpectralDensities.SpectralDensity}
    nbaths = length(Jw)
    c = Vector{Vector{Float64}}(undef, nbaths)
    ω = Vector{Vector{Float64}}(undef, nbaths)
    s = Vector{Vector{Float64}}(undef, nbaths)

    for n in 1:nbaths
        ω[n], c[n] = SpectralDensities.discretize(Jw[n], num_bath_modes[n])
        s[n] = svec[n,:]
    end

    bath = SolventsX.HarmonicBathX(; β, ω, c, svecs=s, nsamples=nmc)
    sys = SpinLSCSys(; transform, Hamiltonian, ρ₀=ρ0, bath, nsamples=nmc)

    propagate_trajectories(solver, sys, dt, ntimes; verbose, kwargs...)
end

end
