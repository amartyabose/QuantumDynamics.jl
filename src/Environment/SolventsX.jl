"Another take at the general description of different types of solvents."
module SolventsX

using ..Solvents: PhaseSpace, Solvent
using Distributions: MvNormal
using LinearAlgebra: diagm

struct HarmonicPhaseSpaceX <: PhaseSpace
    q::Vector{Vector{Float64}}
    p::Vector{Vector{Float64}}
end

struct HarmonicBathX <: Solvent
    β::Float64
    ω::Vector{Vector{Float64}}
    c::Vector{Vector{Float64}}
    s::Vector{Vector{Float64}}
    distq::Vector{MvNormal}
    distp::Vector{MvNormal}
    nsamples::Integer
    nbaths::Integer
end

function HarmonicBathX(; β::Float64, ω::Vector{Vector{Float64}},
                         c::Vector{Vector{Float64}},
                         svecs::Vector{Vector{Float64}},
                         nsamples::Integer)
    @assert length(c) == length(ω)
    @assert length(c) == length(svecs)

    distq = Vector{MvNormal}(undef, length(ω))
    distp = Vector{MvNormal}(undef, length(ω))

    for b in 1:length(ω)
        cth = coth.(0.5 * ω[b] * β)
        distq[b] = MvNormal(diagm(cth ./ ω[b] / 2))
        distp[b] = MvNormal(diagm(cth .* ω[b] / 2))
    end

    HarmonicBathX(β, ω, c, svecs, distq, distp, nsamples, length(ω))
end

function Base.iterate(bath::HarmonicBathX, state=1)
    state > bath.nsamples && return nothing

    (HarmonicPhaseSpaceX(
        map(dist -> rand(dist), bath.distq),
        map(dist -> rand(dist), bath.distp)),
     state+1)
end
Base.eltype(::HarmonicBathX) = HarmonicPhaseSpaceX
Base.length(b::HarmonicBathX) = b.nsamples

raw"""
    propagate_forced_bath(bath::HarmonicBathX, bps::HarmonicPhaseSpaceX,
                          f::Vector{Vector{Float64}}, dt::Real, _::Integer)

Propagate the `bath` subject to constant force `f` from the system.

This assumes the Hamiltonian of the oscillator is of the form
``H = \frac{p^2}{2} + \frac{\omega^2 q^2}{2} - f q.``
"""
function propagate_forced_bath(bath::HarmonicBathX, bps::HarmonicPhaseSpaceX,
                               f::Vector{Vector{Float64}}, dt::Real, _::Integer)
    q = similar.(bps.q)
    p = similar.(bps.p)

    @inbounds for b in eachindex(bps.q)
        sinωt = @. sin(bath.ω[b] * dt)
        cosωt = @. cos(bath.ω[b] * dt)
        disp = @. f[b] / bath.ω[b]^2
        @. q[b] =  (bps.q[b] - disp) * cosωt + bps.p[b] * sinωt / bath.ω[b] + disp
        @. p[b] = -(bps.q[b] - disp) * bath.ω[b] * sinωt + bps.p[b] * cosωt
    end

    HarmonicPhaseSpaceX(q, p)
end

"""
    bath_force(bath::Solvent, state::PhaseSpace, n::Integer)

Return the `n`th bath's (self-)force at the phase space point given by `state`.
"""
function bath_force(bath::Solvent, state::PhaseSpace, n::Integer) end

bath_force(bath::HarmonicBathX, state::HarmonicPhaseSpaceX, n::Integer) = -bath.ω[n]^2 * state.q[n]

raw"""
    propagate_forced_bath(bath::Solvent, state::PhaseSpace,
                          f::Vector{Vector{Float64}}, dt::Real,
                          ntimes::Integer)

Propagate `bath` subject to constant force `f` from the system.

The Hamiltonian for a single bath is assumed to be of the form
``H = \sum_i \frac{p_i^2}{2} + V(\{q_i\}) - f_i q_i``
where `V(\{q_i\})` is the state-independent potential of the bath.

The time step of the propagation is given by `dt`, and the number of
Verlet steps by `ntimes`.
"""
function propagate_forced_bath(bath::Solvent, state::PhaseSpace,
                               f::Vector{Vector{Float64}}, dt::Real,
                               ntimes::Integer)
    q = similar.(state.q)
    p = similar.(state.p)

    @inbounds for _ in 1:ntimes
        for b in eachindex(state.q)
            @. p[b] += 0.5 * (bath_force(bath, state, b) + f[b]) * dt
            @. q[b] += p[b] * dt
        end
        state = eltype(bath)(q, p)
        for b in eachindex(state.p)
            @. p[b] += 0.5 * (bath_force(bath, state, b) + f[b]) * dt
        end
    end

    eltype(bath)(q, p)
end

end
