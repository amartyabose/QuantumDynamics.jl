"Collection of spectral densities commonly used to describe solvents."
module SpectralDensities

using DelimitedFiles
using ..Utilities

abstract type SpectralDensity end

abstract type ContinuousSpectralDensity <: SpectralDensity end
abstract type DiscreteOscillators <: SpectralDensity end

abstract type AnalyticalSpectralDensity <: ContinuousSpectralDensity end
(sd::AnalyticalSpectralDensity)(ω::Real) = evaluate(sd, ω)
eval_spectrum(sd::AnalyticalSpectralDensity, ω::Real, β::Real) = ω == 0.0 ? eval_spectrum_at_zero(sd) : 2.0 * sd(ω) / (1 - exp(-β * ω))

struct ExponentialCutoff <: AnalyticalSpectralDensity
    ξ::Real
    ωc::Real
    Δs::Real
    n::Real
    ωmax::Real
end
"""
    ExponentialCutoff(; ξ, ωc, n=1.0, Δs=2.0)
Construct a model spectral density with an exponential cutoff.

``J(ω) = \\frac{2π}{Δs^2} ξ \\frac{ω^n}{ω_c^{n-1}} \\exp\\left(-\\frac{ω}{ωc}\\right)``

where `Δs` is the distance between the two system states. The model is Ohmic if `n = 1`, sub-Ohmic if `n < 1`, and super-Ohmic if `n > 1`.
"""
ExponentialCutoff(; ξ::Float64, ωc::Float64, n=1.0, Δs=2.0) = ExponentialCutoff(ξ, ωc, Δs, n, 30 * ωc)
evaluate(sd::ExponentialCutoff, ω::Real) = 2π / sd.Δs^2 * sd.ξ * sign(ω) * abs(ω)^sd.n * sd.ωc^(1 - sd.n) * exp(-abs(ω) / sd.ωc)
eval_spectrum_at_zero(sd::ExponentialCutoff) = sd.n == 1 ? 2.0 * 2π / sd.Δs^2 * sd.ξ : 0

function discretize(sd::ExponentialCutoff, num_osc::Int)
    @assert sd.n == 1 "Discretization works only for sd.n=1"
    ω = zeros(num_osc)
    c = zeros(num_osc)
    if sd.n != 1
        return ω, c
    end

    ωmax = 5 * sd.ωc
    for i = 1:num_osc
        ω[i] = -sd.ωc * log(1 - i * (1 - exp(-ωmax / sd.ωc)) / num_osc)
        c[i] = sqrt(sd.ξ * sd.ωc * (1 - exp(-ωmax / sd.ωc)) / num_osc) * ω[i]
    end
    ω, c
end

struct DrudeLorentz <: AnalyticalSpectralDensity
    λ::Real
    γ::Real
    Δs::Real
    ωmax::Real
end
"""
    DrudeLorentz(; λ, γ, Δs=2.0)
Construct a model spectral density with a Drude-Lorentz cutoff.

``J(ω) = \\frac{2λ}{Δs^2} \\frac{ω γ}{ω^2 + γ^2}``

where `Δs` is the distance between the two system states.
"""
DrudeLorentz(; λ::Float64, γ::Float64, Δs=2.0, ωmax=1000 * γ) = DrudeLorentz(λ, γ, Δs, ωmax)
evaluate(sd::DrudeLorentz, ω::Real) = 2 * sd.λ / sd.Δs^2 * sign(ω) * abs(ω) * sd.γ / (abs(ω)^2 + sd.γ^2)
eval_spectrum_at_zero(sd::DrudeLorentz) = 2.0 * 2 * sd.λ / sd.Δs^2 * sd.γ

"""
    matsubara_decomposition(sd::DrudeLorentz, num_modes::Int, β::Float64)

Returns the Matsubara frequencies, `γ`, and the expansion coefficients, `c`.
"""
function matsubara_decomposition(sd::DrudeLorentz, num_modes::Int, β::Float64)
    γ = zeros(num_modes + 1)
    c = zeros(ComplexF64, num_modes + 1)
    γ[1] = sd.γ
    c[1] = sd.λ * sd.γ / sd.Δs^2 * (cot(β * sd.γ / 2.0) - 1im)
    for k = 2:num_modes+1
        γ[k] = 2 * (k - 1) * π / β
        c[k] = 4 * sd.λ / sd.Δs^2 * sd.γ / β * γ[k] / (γ[k]^2 - sd.γ^2)
    end

    γ, c
end

function discretize(sd::DrudeLorentz, num_osc::Int)
    js = 1:num_osc
    ω = sd.γ .* tan.(π / 2.0 .* (js .- 0.5) ./ num_osc)
    c = sqrt(sd.λ / (2 * num_osc)) .* ω
    ω, c
end

function tabulate(sd::T, full_real::Bool=true, npoints::Int=100001) where {T<:AnalyticalSpectralDensity}
    ω = full_real ? range(-sd.ωmax, stop=sd.ωmax, step=2 * sd.ωmax / npoints) : range(sd.ωmax / npoints, stop=sd.ωmax, step=sd.ωmax / npoints)
    ω, sd.(ω)
end

function reorganization_energy(sd::AnalyticalSpectralDensity)
    ω, jw = tabulate(sd)
    jw ./= ω
    Utilities.trapezoid(ω, jw) / 2π * sd.Δs^2
end


"""
    SpectralDensityTable <: ContinuousSpectralDensity

Spectral density provided in tabular form. Contains a vector of `ω`s and a vector corresponding to `jw`s.
"""
struct SpectralDensityTable <: ContinuousSpectralDensity
    ω::Vector{Float64}
    jw::Vector{Float64}
end
function read_jw(filename, delim; skipstart=0)
    w_jw = readdlm(filename, delim; skipstart)
    SpectralDensityTable(w_jw[:, 1], w_jw[:, 2])
end
function read_jw_over_w(filename, delim; skipstart=0)
    w_jw = readdlm(filename, delim; skipstart)
    SpectralDensityTable(w_jw[:, 1], w_jw[:, 2] .* w_jw[:, 1])
end

function tabulate(sd::SpectralDensityTable, full_real::Bool=true, npoints::Int=100001)
    sd.ω, sd.jw
end

function reorganization_energy(sd::SpectralDensityTable)
    jw = sd.jw ./ sd.ω
    Utilities.trapezoid(sd.ω, jw) / π
end

end