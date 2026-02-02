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

end
