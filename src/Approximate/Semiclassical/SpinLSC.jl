module SpinLSC

using HDF5
using ..Utilities
using ..SolventsX
using ..System
using ..SpectralDensities
using LinearAlgebra: diag, diagm
import OrdinaryDiffEq as ODE


const references = """
"""

"Abstract type for all Stratonovich–Weyl transforms."
abstract type SWTransform end

abstract type WTransform <: SWTransform end
abstract type QTransform <: SWTransform end
abstract type PTransform <: SWTransform end

R²(::Type{WTransform}, d::Integer) = 2sqrt(d+1)
R²(::Type{QTransform}, d::Integer) = 2.0
R²(::Type{PTransform}, d::Integer) = 2.0 * (d+1)

γ(t::Type{<:SWTransform}, d::Integer) = (R²(t, d) - 2.0) / d
γ(::Type{QTransform}, d::Integer) = 0.0

dual(W::Type{WTransform}) = W
dual(::Type{PTransform}) = QTransform
dual(::Type{QTransform}) = PTransform

rescale_factor(s::Type{<:SWTransform}, s̄::Type{<:SWTransform}, d::Integer) = sqrt(R²(s̄, d) / R²(s, d))
rescale_factor(::Type{WTransform}, ::Type{WTransform}, ::Integer) = 1



struct SpinMappedSysPhaseSpace <: PhaseSpace
    X::Vector{Float64}
    P::Vector{Float64}
end

struct SpinMappedSystem <: CompositeSystem
    transform::Type{<:SWTransform}
    h::AbstractMatrix
    ρ₀::Union{Nothing,AbstractMatrix{<:Complex}}
    R²::Float64
    γₛ::Float64
    d::Integer
    bath::Solvent
    nsamples::Integer
end
function SpinMappedSystem(; transform::Type{<:SWTransform}, hamiltonian::AbstractMatrix,
                          ρ₀::AbstractMatrix, bath::Solvent, nsamples::Integer)
    @assert nsamples == bath.nsamples
    d = size(hamiltonian, 1)
    SpinMappedSystem(transform, hamiltonian, ρ₀, R²(transform, d), γ(transform, d), d, bath, nsamples)
end

function Base.iterate(sys::SpinMappedSystem, state=1)
    state > sys.nsamples && return nothing

    bathps = next(sys.bath)

    R = sqrt(sys.R²)
    X = randn(d)
    P = randn(d)

    sqΣ = sqrt(sum(X.^2 .+ P.^2))

    SpinMappedSysPhaseSpace(X * R / sqΣ, P * R / sqΣ), bathps
end
Base.eltype(::SpinMappedSystem) = SpinMappedSysPhaseSpace
Base.length(s::SpinMappedSystem) = s.nsamples



"Perform the relevant SW transform of the operator `op`."
function transform_op(sys::SpinMappedSystem, op::AbstractMatrix, ps::SpinMappedSysPhaseSpace)
    opt = 0.5 * sum(diag(op) .* (ps.X.^2 .+ ps.P.^2 .- sys.γₛ))

    for n = 1:sys.d, m = n+1:sys.d
        opt += 0.5 * op[n,m] * (ps.X[n] - im * ps.P[n]) * (ps.X[m] + im * ps.P[m])
        opt += 0.5 * op[m,n] * (ps.X[m] - im * ps.P[m]) * (ps.X[n] + im * ps.P[n])
    end

    opt
end

"Calculate the force on the `i`th bath."
function Fbath(sys::SpinMappedSystem, sps::SpinMappedSysPhaseSpace,
               q::AbstractVector{Float64}, i::Integer)
    -sys.bath.ω[i].^2 * q .+
        transform_op(sys, diagm(sys.bath.s[i]), sps) * sys.bath.c[i]
end

function H(sys::SpinMappedSystem, sps::SpinMappedSysPhaseSpace, bps::PhaseSpace)
    bs = sys.bath
    transform_op(sys, sys.h, sps) +
        mapreduce((b, x, p) -> 0.5 * p.^2 + 0.5 * bs.ω[b].^2 * x.^2, +, 1:sys.bath.nbaths, bps.q, bps.p) +
        -mapreduce((b,x) -> transform_op(sys, sum(bs.c[b] .* x) * diagm(bs.s[b]), sps), +, 1:sys.bath.nbaths, bps.q)
end

"""Reconstruct the bare density matrix for the phase space point.

This does NOT multiply the Stratonovich–Weyl transform of the initial
density matrix.  Intended to be used when constructing the dynamical
maps.
"""
function reconstruct_bare_ρ(sys::SpinMappedSystem, sps::SpinMappedSysPhaseSpace)
    ρ = zeros(ComplexF64, sys.d,sys.d)

    dual_trans = dual(sys.transform)
    rescale = rescale_factor(sys.transform, dual_trans, sys.d)
    γs̄ = γ(dual_trans, sys.d)

    X̄ = rescale * sps.X
    P̄ = rescale * sps.P

    for j = 1:sys.d
        # Some formulae:
        # \[ [\ket{k}\bra{j}]_s = \frac{(X_k - i P_k) (X_j + i P_j)}{2} \]
        # \[ [\ket{j}\bra{j}]_s = \frac{X_j^2 + P_j^2 - \gamma_s}{2} \]
        # with the approriate rescaling of X_j's.

        ρ[j,j] = 0.5 * (X̄[j].^2 + P̄[j].^2 - γ)
        for k = j+1:sys.d
            ρ[j,k] = 0.5 * (X̄[k] - im * P̄[k]) * (X̄[j] + im * P̄[j])
            ρ[k,j] = 0.5 * (X̄[j] - im * P̄[j]) * (X̄[k] + im * P̄[k])
        end
    end

    ρ * d
end



abstract type SpinLSCSolver end
abstract type RK4 <: SpinLSCSolver end

function HX(sys::SpinMappedSystem, sps::SpinMappedSysPhaseSpace, bps::PhaseSpace)
    ddX = zeros(sys.d)

    bs = sys.bath
    for i in 1:sys.d
        # NOTE: The diagonal term of the Hamiltonian enters naturally
        # from the summation.
        ddX[i] = 0.5 * sum(sys.h[i,:] .* (sps.X .+ im * sys.P)) +
            0.5 * sum(sys.h[:,i] .* (sys.X .- im * sys.P)) +
            -mapreduce((b,x) -> sum(bs.c[b] .* x) * bs.s[b][i] * X[i], +, 1:sys.bath.nbaths, bps.q)
    end

    ddX
end

function HP(sys::SpinMappedSystem, sps::SpinMappedSysPhaseSpace, bps::PhaseSpace)
    ddP = zeros(sys.d)

    bs = sys.bath
    for i in 1:sys.d
        # NOTE: The diagonal term of the Hamiltonian enters naturally
        # from the summation.
        ddP[i] = 0.5 * sum(sys.h[i,:] .* (sps.P .+ im * sps.X)) +
            0.5 * (sum.h[:,i] .* (sps.P .- im * sps.X)) +
            -mapreduce((b, x) -> sum(bs.c[b] .* x) * bs.s[b][i] * X[i], +,
                       1:sys.bath.nbaths, bps.q)
    end

    ddP
end

function propagate_xpXP!(du, u, p, t)
    sys, xis, pis = p
    d = sys.d

    sps = SpinMappedSysPhaseSpace(u[1:d], u[d+1:2d])
    bps = HarmonicPhaseSpaceX([ u[i] for i in xis ], [ u[i] for i in pis ])

    du[1:d] = HP(sys, sps, bps)     # Ẋ
    du[d+1:2d] = -HX(sys, sps, bps) # Ṗ

    for n in 1:sys.bath.nbaths
        du[xis[n]] = bps.p[n]   # ẋ
        du[pis[n]] = Fbath(sys, sps, bps.x[n], n)
    end
end

function bathinds(sys::SpinMappedSystem)
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

function build_dynmap(sol::ODE.ODESolution)
    sys, _ = sol.prob.p
    Nₜ = length(sol.u)
    d = sys.d

    U0e = zeros(ComplexF64, Nₜ,2d,2d)

    sps0 = SpinMappedSysPhaseSpace(sol.u[1][1:d], sol.u[1][d+1:2d])
    for t in 1:N
        sps = SpinMappedSysPhaseSpace(sol.u[t][1:d], sol.u[t][d+1:2d])
        bareρ = reconstruct_bare_ρ(sys, sps)

        for i in 1:2d
            ρᵥ = zeros(ComplexF64, 2d)
            ρᵥ[i] = 1.0
            ρₛ = transform_op(sys,
                              Utilities.density_matrix_vector_to_matrix(ρᵥ),
                              sps0)
            U0e[t,:,i] = Utilities.density_matrix_to_vector(ρₛ * bareρ)
        end
    end

    U0e
end

function reconstruct_ρ(sol::ODE.ODESolution)
    sys, _ = sol.prob.p
    Nₜ = length(sol.u)
    d = sys.d

    ρ = zeros(ComplexF64, Nₜ,sys.d,sys.d)
    sps0 = SpinMappedSysPhaseSpace(sol.u[1][1:d], sol.u[1][d+1:2d])
    ρₛ = transform_op(sys, sys.ρ₀, sps0)
    for t in 1:N
        ρ[t,:,:] = ρₛ * reconstruct_bare_ρ(sys,
                                           SpinMappedSysPhaseSpace(
                                               sol.u[t][1:d],
                                               sol.u[t][d+1:2d]))
    end

    ρ
end

function propagate_trajectories(::Type{RK4}, sys::SpinMappedSystem, dt::Float64, ntimes::Integer;
                                output::Union{Nothing,HDF5.Group}=nothing, verbose::Bool=false, kwargs...)
    xis, pis = bathinds(sys)

    outputρ = if !isnothing(output) && haskey(kwargs, :outgroup)
        Utilites.create_and_select_group(output, kwargs[:outgroup])
    else
        nothing
    end

    if !isnothing(output)
        Utilities.check_or_insert_value(output, "U0e", zeros(ComplexF64, ntimes+1,2sys.d,2sys.d))
        Utilities.check_or_insert_value(output, "samples_done", 0)
        !isnothing(sys.ρ₀) && !isnothing(outputρ) &&
            Utilities.check_or_insert_value(outputρ, "rho", zeros(ComplexF64, ntimes+1,sys.d,sys.d))
    end

    probfn(p, _, _) = begin
        sps, bps = next(sys)
        ODE.remake(p, u0=vcat(sps.X, sps.P, bps.q..., bps.p...))
    end
    outputfn(sol, _) = begin
        (build_dynmap(sol),
         isnothing(sys.ρ₀) ? nothing, reconstruct_ρ(sol))
    end

    done = 0
    reducefn(data, us, I) = begin
        done += length(I)
        verbose && @info "Trajectories completed: $(done * 100 / N)%"

        if !isnothing(output)
            output["U0e"] = (data[1] + sum(getindex.(us, 1))) / done
            output["samples_done"] = done
            flush(output)
        end

        if !isnothing(data[2]) && !isnothing(outputρ)
            outputρ["rho"] = (data[2] + sum(getindex.(us, 2))) / done
            flush(outputρ)
        end

        if isnothing(us[1][end])
            (data[1] + sum(getindex.(us, 1)), nothing)
        else
            (data[1] + sum(getindex.(us, 1)), data[2] + sum(getindex.(us, 2)))
        end
    end

    ensemble = ODE.EnsembleProblem(
        ODE.ODEProblem(propagate_xpXP!,
                       zeros(2d+2sum(length.(sys.bath.ω))),
                       (0.0, ntimes*dt),
                       (sys, xis, pis));
        output_func=outputfn,
        prob_func=probfn,
        reduction=reducefn,
        u_init=(zeros(ComplexF64, ntimes+1,2sys.d,2sys.d),
                isnothing(sys.ρ₀) ? nothing : zeros(ComplexF64, ntimes+1,sys.d,sys.d)))

    sol = ODE.solve(ensemble, ODE.RK4(), ODE.EnsembleThreads();
                    dt=dt, saveat=dt, trajectories=length(sys),
                    batch_size=Threads.nthreads()*2)

    (sol.u[1] / length(sys),
     isnothing(sys.ρ₀) ? nothing : sol.u[2] / length(sys))
end



abstract type Verlet <: SpinLSCSolver end

function propagate_trajectory(::Type{Verlet}, sys::SpinMappedSystem,
                              sps0::SpinMappedSysPhaseSpace,
                              bps0::PhaseSpace,
                              dt::Float64, ntimes::Integer)
    X = sps0.X
    P = sps0.P
    x = bps0.q
    p = bps0.p
    d = sys.d

    U0e = zeros(ComplexF64, ntimes+1,2d,2d)
    isnothing(sys.ρ₀) || (ρ = zeros(ComplexF64, ntimes+1,d,d))

    build_dynmap_ρ!(t) = begin
        sps = SpinMappedSysPhaseSpace(X, P)
        bareρ = reconstruct_bare_ρ(sys, sps)

        for i in 1:2d
            ρᵥ = zeros(ComplexF64, 2d)
            ρᵥ[i] = 1.0
            ρₛ = transform_op(sys,
                              Utilities.density_matrix_vector_to_matrix(ρᵥ),
                              sps0)
            U0e[t,:,i] = Utilities.density_matrix_to_vector(ρₛ * bareρ)
        end

        if !isnothing(sys.ρ₀)
            ρ[t,:,:] = transform_op(sys, sys.ρ₀, sps0) * bareρ
        end
    end

    δtₓ = dt / 100
    N½ = dt / 2 / δtₓ
    propagate_xp!(t) = for b in 1:sys.bath.nbaths
        for _ in 1:N½
            p[b] = p[b] .+ 0.5 * Fbath(sys, SpinMappedSysPhaseSpace(X[:,t], P[:,t]), x[b], b) * δtₓ
            x[b] = x[b] .+ p[b] * δtₓ
            p[b] = p[b] .+ 0.5 * Fbath(sys, SpinMappedSysPhaseSpace(X[:,t], P[:,t]), x[b], b) * δtₓ
        end
    end

    bs = sys.bath
    LXP = zeros(2d,2d)
    for t in 2:ntimes+1
        propagate_xp!(t)

        V = sys.h - mapreduce((b, x) -> sum(bs.c[b] .* x) * diagm(bs.s[b]), +, bs, x)
        LXP[1:d,d+1:2d] = V
        LXP[d+1:2d,1:d] = -V
        XPt = exp(LXP * dt) * [ X, P ]
        X, P = Xpt[1:d], Xpt[d+1:2d]

        propagate_xp!(t)

        build_dynmap_ρ!(t)
    end

    U0e, isnothing(sys.ρ₀) ? nothing : ρ
end

function propagate_trajectories(::Type{Verlet}, sys::SpinMappedSystem, dt::Float64, ntimes::Integer;
                                output=Union{Nothing,HDF5.Group}=nothing, verbose::Bool=false, kwargs...)
    U0e = zeros(ComplexF64, ntimes+1,2sys.d,2sys.d)
    isnothing(sys.ρ₀) || (ρ = zeros(ComplexF64, ntimes+1,sys.d,sys.d))

    nthreads = Threads.nthreads()

    outputρ = if !isnothing(output) && haskey(kwargs, :outgroup)
        Utilities.create_and_select_group(output, kwargs[:outgroup])
    else
        nothing
    end

    if !isnothing(output)
        Utilities.check_or_insert_value(output, "U0e", U0e)
        Utilities.check_or_insert_value(output, "samples_done", 0)
        !isnothing(sys.ρ₀) && !isnothing(outputρ) &&
            Utilities.check_or_insert_value(outputρ, "rho", ρ)
    end

    done = 0
    for i in 1:cld(length(sys), threads)
        tasks = map(1:nthreads) do j
            sps0, bps0 = next(sys)
            Threads.@spawn(propagate_trajectory(Verlet, sys, sps0, bps0))
        end
        solns = fetch.(tasks)
        done += nthreads

        Ue0 += sum(solns[1])
        if !isnothing(output)
            output["samples_done"] = done
            output["U0e"] = U0e / done
            flush(output)
        end

        if !isnothing(sys.ρ₀)
            ρ += sum(solns[2])
            if !isnothing(outputρ)
                outputρ["rho"] = ρ / done
                flush(outputρ)
            end
        end
        verbose && @info "Trajectories complete: $(done * 100 / length(sys))%"
    end

    Ue0 / length(sys), isnothing(sys.ρ₀) ? nothing : ρ / length(sys)
end



"""
    propagate(; Hamiltonian::Matrix{<:Complex}, Jw::Vector{T},
             β::Real, num_bath_modes::Vector{<:Int}, svec::Matrix{Float64},
             ρ0::Union{Nothing,Matrix{<:Complex}}, dt::Real,
             ntimes::Real, transform::Type{<:SWTransform},
             nmc::Int, solver::Type{<:SpinLSCSolver}, verbose::Bool=false,
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
                   β::Real, num_bath_modes::Vector{<:Int}, svec::Matrix{Float64},
                   ρ0::Union{Nothing,Matrix{<:Complex}}, dt::Real,
                   ntimes::Real, transform::Type{<:SWTransform},
                   nmc::Int, solver::Type{<:SpinLSCSolver}, verbose::Bool=false,
                   kwargs...) where {T<:SpectralDensities.SpectralDensity}
    nbaths = length(Jw)
    c = Vector{Vector{Float64}}(undef, nbaths)
    ω = Vector{Vector{Float64}}(undef, nbaths)
    s = Vector{Vector{Float64}}(undef, nbaths)

    for n in 1:nbaths
        ω[n], c[n] = SpectralDensities.discretize(Jw[n], num_bath_modes[n])
        s[n] = svec[n,:]
    end

    bath = HarmonicBathX(β, ω, c, svecs=s, nsamples=nmc)
    sys = SpinMappedSystem(transform, hamiltonian=Hamiltonian, ρ₀=ρ0, bath, nsamples=nmc)

    propagate_trajectories(solver, sys, dt, ntimes; verbose, kwargs...)
end

end
