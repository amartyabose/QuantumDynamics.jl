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
    num_modes::Int
    ω::Vector{Float64}
    c::Vector{Float64}
    sys_op::Matrix{Float64}
    distq::FullNormal
    distp::FullNormal
    num_samples::Int
    eqm_center::Vector{Float64}
end
function HarmonicBath(β::Float64, ωs::Vector{Vector{Float64}}, cs::Vector{Vector{Float64}}, sys_ops::AbstractMatrix{Float64}, num_samples::Int)
    @assert length(ωs) == length(cs)
    @assert length(ωs) == size(sys_ops, 1)
    nmodes = 0
    for ω in ωs
        nmodes += length(ω)
    end
    ω = vcat(ωs...)
    c = vcat(cs...)
    sys_op = zeros(nmodes, size(sys_ops, 2))
    ind = 1
    for (j, ω) in enumerate(ωs)
        for _ = 1:length(ω)
            sys_op[ind, :] = sys_ops[j, :]
            ind += 1
        end
    end
    σp2 = ω ./ 2 .* coth.(ω .* β ./ 2)
    distp = MvNormal(repeat([0.0], nmodes), diagm(σp2))
    σq2 = coth.(ω .* β ./ 2) ./ (2 .* ω)
    distq = MvNormal(repeat([0.0], nmodes), diagm(σq2))
    eqm_center = c ./ ω .^ 2
    HarmonicBath(β, nmodes, ω, c, sys_op, distq, distp, num_samples, eqm_center)
end
function Base.iterate(hb::HarmonicBath, state=1)
    state <= hb.num_samples ? (HarmonicPhaseSpace(rand(hb.distq, 1), rand(hb.distp, 1)), state + 1) : nothing
end
Base.eltype(hb::HarmonicBath) = HarmonicPhaseSpace
Base.length(hb::HarmonicBath) = hb.num_samples
function propagate_trajectory(hb::HarmonicBath, state::HarmonicPhaseSpace, dt::Float64, ntimes::Int64, ref_pos_mod::Vector{Float64})
    q0, p0 = state.q, state.p
    q = copy(q0)
    p = copy(p0)
    nsys = size(hb.sys_op, 2)
    bath_contribution = 0.5 * transpose(hb.sys_op) .^ 2 * (hb.c .^ 2 ./ hb.ω .^ 2) |> collect
    # bath_contribution = zeros(nsys)
    # for (j, (ω, c)) in enumerate(zip(hb.ω, hb.c))
    #     bath_contribution .+= 0.5 * c^2 / ω^2 * hb.sys_op[j, :].^2
    # end
    # @show bath_contribution
    energy = zeros(Float64, ntimes + 1, nsys)
    dedt = zeros(Float64, ntimes + 1, nsys)
    for j = 1:ntimes+1
        energy[j, :] .= bath_contribution
    end
    energy[1, :] .+= -transpose(hb.sys_op) * (hb.c .* q0) |> collect
    dedt[1, :] .+= -transpose(hb.sys_op) * (hb.c .* p0) |> collect
    for i = 1:ntimes
        q .= (q0 .- ref_pos_mod .* hb.eqm_center) .* cos.(hb.ω .* i .* dt) .+ p0 ./ hb.ω .* sin.(hb.ω .* i .* dt) .+ ref_pos_mod .* hb.eqm_center
        p .= p0 .* cos.(hb.ω .* i .* dt) .- (q0 .- ref_pos_mod .* hb.eqm_center) .* hb.ω .* sin.(hb.ω .* i .* dt)
        energy[i+1, :] .+= -transpose(hb.sys_op) * (hb.c .* q) |> collect
        dedt[i+1, :] .+= -transpose(hb.sys_op) * (hb.c .* p) |> collect
    end
    energy, dedt, HarmonicPhaseSpace(q, p)
end
dedt(hb::HarmonicBath, state::HarmonicPhaseSpace) = -sum(hb.c .* state.p) .* hb.sys_op

end