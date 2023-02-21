"General description of different types of solvents."
module Solvents

using Distributions, LinearAlgebra

"Abstract type for every solvent. Every solvent needs to implement `Base.iterate`, which returns the next sample of the phase-space, and a `propagate_trajectory` method which propagates a given phasespace point using classical mechanics."
abstract type Solvent end

"Abstract type for all phase spaces. Each `Solvent` has an associated phase-space."
abstract type PhaseSpace end

struct HarmonicPhaseSpace <: PhaseSpace
    q::AbstractMatrix{Float64}
    p::AbstractMatrix{Float64}
end

struct HarmonicBath <: Solvent
    β::Float64
    ω::Vector{Float64}
    c::Vector{Float64}
    sys_op::Vector{Float64}
    distq::FullNormal
    distp::FullNormal
    num_samples::Int
    eqm_center::Vector{Float64}
end
function HarmonicBath(β::Float64, ω::Vector{Float64}, c::Vector{Float64}, sys_op::Vector{Float64}, num_samples::Int)
    nmodes = length(ω)
    σp2 = ω ./ 2 .* coth.(ω .* β ./ 2)
    distp = MvNormal(repeat([0.0], nmodes), diagm(σp2))
    σq2 = coth.(ω .* β ./ 2) ./ (2 .* ω)
    distq = MvNormal(repeat([0.0], nmodes), diagm(σq2))
    eqm_center = c ./ ω .^ 2
    HarmonicBath(β, ω, c, sys_op, distq, distp, num_samples, eqm_center)
end
function Base.iterate(hb::HarmonicBath, state=1)
    state <= hb.num_samples ? (HarmonicPhaseSpace(rand(hb.distq, 1), rand(hb.distp, 1)), state + 1) : nothing
end
Base.eltype(hb::HarmonicBath) = HarmonicPhaseSpace
Base.length(hb::HarmonicBath) = hb.num_samples
function propagate_trajectory(hb::HarmonicBath, state::HarmonicPhaseSpace, dt::Float64, ntimes::Int64, ref_pos_mod)
    q0, p0 = state.q, state.p
    q = copy(q0)
    p = copy(p0)
    nsys = length(hb.sys_op)
    bath_contribution = 0.5 * sum(hb.c .^ 2 ./ hb.ω .^ 2) .* hb.sys_op .^ 2
    energy = zeros(Float64, ntimes + 1, nsys)
    energy[1, :] .= -sum(hb.c .* q0) .* hb.sys_op .- bath_contribution
    for i = 1:ntimes
        q .= (q0 .- ref_pos_mod[i] .* hb.eqm_center) .* cos.(hb.ω .* dt) .+ p0 ./ hb.ω .* sin.(hb.ω .* dt) .+ ref_pos_mod[i] .* hb.eqm_center
        p .= p0 .* cos.(hb.ω .* dt) .- (q0 .- ref_pos_mod[i] .* hb.eqm_center) .* hb.ω .* sin.(hb.ω .* dt)
        energy[i+1, :] .= -sum(hb.c .* q) * hb.sys_op .- bath_contribution
        q0 .= q
        p0 .= p
    end
    energy, HarmonicPhaseSpace(q0, p0)
end

end