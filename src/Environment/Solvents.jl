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
    num_baths::Int
    ω::Vector{Float64}
    c::Vector{Float64}
    sys_op::Matrix{Float64}
    distq::FullNormal
    distp::FullNormal
    num_samples::Int
    eqm_center::Vector{Float64}
    mappings::Vector{Int64}
end
function HarmonicBath(β::Float64, ωs::Vector{Vector{Float64}}, cs::Vector{Vector{Float64}}, sys_ops::AbstractMatrix{Float64}, num_samples::Int)
    @assert length(ωs) == length(cs)
    @assert length(ωs) == size(sys_ops, 1)
    ω = vcat(ωs...)
    c = vcat(cs...)
    mappings = reduce(vcat, [repeat([j], length(ωs[j])) for j in eachindex(ωs)])
    σp2 = ω ./ 2 .* coth.(ω .* β ./ 2)
    distp = MvNormal(repeat([0.0], length(ω)), diagm(σp2))
    σq2 = coth.(ω .* β ./ 2) ./ (2 .* ω)
    distq = MvNormal(repeat([0.0], length(ω)), diagm(σq2))
    eqm_center = c ./ ω .^ 2
    HarmonicBath(β, length(ωs), ω, c, sys_ops, distq, distp, num_samples, eqm_center, mappings)
end
function Base.iterate(hb::HarmonicBath, state=1)
    state <= hb.num_samples ? (HarmonicPhaseSpace(rand(hb.distq, 1), rand(hb.distp, 1)), state + 1) : nothing
end
Base.eltype(::HarmonicBath) = HarmonicPhaseSpace
Base.length(hb::HarmonicBath) = hb.num_samples
function propagate_trajectory(hb::HarmonicBath, state::HarmonicPhaseSpace, dt::Float64, ntimes::Int64, ref_pos_mod::Vector{Float64})
    q0, p0 = state.q, state.p
    q = copy(q0)
    p = copy(p0)
    nsys = size(hb.sys_op, 2)
    energy = zeros(Float64, ntimes + 1, nsys)
    bath_contribution = zeros(nsys)
    @inbounds begin
        for (j, ω, c) in zip(hb.mappings, hb.ω, hb.c)
            bath_contribution .+= 0.5 * c^2 / ω^2 * hb.sys_op[j, :].^2
        end

        for j = 1:ntimes+1
            energy[j, :] .= bath_contribution
        end
        for (j, nb) in enumerate(hb.mappings)
            energy[1, :] += - hb.sys_op[nb, :]  * hb.c[j] * q[j]
        end
        for i = 1:ntimes
            for (j, nb) in enumerate(hb.mappings)
                q[j] = (q0[j] - ref_pos_mod[nb] * hb.eqm_center[j]) * cos(hb.ω[j] * i * dt) + p0[j] / hb.ω[j] * sin(hb.ω[j] * i * dt) + ref_pos_mod[nb] * hb.eqm_center[j]
                p[j] = p0[j] * cos(hb.ω[j] * i * dt) - (q0[j] - ref_pos_mod[nb] * hb.eqm_center[j]) * hb.ω[j] * sin(hb.ω[j] * i * dt)
                energy[i+1, :] += - hb.sys_op[nb, :]  * hb.c[j] * q[j]
            end
        end
    end
    energy, HarmonicPhaseSpace(q, p)
end
# dedt(hb::HarmonicBath, state::HarmonicPhaseSpace) = -sum(hb.c .* state.p) .* hb.sys_op

end