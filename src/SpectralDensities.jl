"Collection of spectral densities commonly used to describe solvents."
module SpectralDensities

abstract type SpectralDensity end

abstract type AnalyticalSpectralDensity <: SpectralDensity end
(sd::AnalyticalSpectralDensity)(ω::Real) = evaluate(sd, ω)
eval_spectrum(sd::AnalyticalSpectralDensity, ω::Real, β::Real) = ω==0.0 ? eval_spectrum_at_zero(sd) : 2.0 * sd(ω) / (1 - exp(-β * ω))

abstract type TabularSpectralDensity <: SpectralDensity end

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
eval_spectrum_at_zero(sd::ExponentialCutoff) = sd.n==1 ? 2.0 * 2π / sd.Δs^2 * sd.ξ : 0

struct DrudeLorentzCutoff <: AnalyticalSpectralDensity
    λ::Real
    γ::Real
    Δs::Real
    n::Real
    ωmax::Real
end
"""
    DrudeLorentzCutoff(; λ, γ, n=1.0, Δs=2.0)
Construct a model spectral density with a Drude-Lorentz cutoff.

``J(ω) = \\frac{2λ}{Δs^2} \\frac{ω^n γ^{2-n}}{ω^2 + γ^2}``

where `Δs` is the distance between the two system states. The model is Ohmic if `n = 1`, sub-Ohmic if `n < 1`, and super-Ohmic if `n > 1`.
"""
DrudeLorentzCutoff(; λ::Float64, γ::Float64, n=1.0, Δs=2.0) = DrudeLorentzCutoff(λ, γ, Δs, n, 100 * ωc)
evaluate(sd::DrudeLorentzCutoff, ω::Real) = 2 * sd.λ / sd.Δs^2 * sign(ω) * abs(ω)^sd.n * sd.γ^(2 - sd.n) / (abs(ω)^2 + sd.γ^2)
eval_spectrum_at_zero(sd::DrudeLorentzCutoff) = sd.n==1 ? 2.0 * 2 * sd.λ / sd.Δs^2 * sd.γ : 0

function tabulate(sd::T, full_real::Bool=true) where {T<:AnalyticalSpectralDensity}
    ω = full_real ? range(-sd.ωmax, stop=sd.ωmax, step=2 * sd.ωmax / 100001) : range(sd.ωmax / 100001, stop=sd.ωmax, step=sd.ωmax / 100001)
    ω, sd.(ω)
end

struct SpectralDensityTable <: TabularSpectralDensity
    ω::Vector{Float64}
    j::Vector{Float64}
end

end