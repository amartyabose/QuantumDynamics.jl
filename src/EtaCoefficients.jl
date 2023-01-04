module EtaCoefficients

using ..SpectralDensities, ..Utilities

function common_part(ω, sd, β, classical)
    if isnan(β) && !classical
        ans = zero(ω)
        for (i, freq) in enumerate(ω)
            if freq > 0
                ans[i] = sd ./ ω .^ 2.0 .* 2.0
            else
                ans[i] = 0
            end
        end
        return ans
    elseif !classical
        return sd ./ ω .^ 2 .* (2.0 ./ (1.0 .- exp.(-ω .* β)))
    elseif classical
        return sd ./ ω .^2 .* (2.0 / (ω .* β) + 1.0)
    end
end

"""
EtaCoefficients holds the various discretized η-coefficients required for a QuAPI-based simulation. These are the minimum number of coefficients required, stored using time-translational symmetry wherever possible.

The values are stored as follows:
`η00`: The self-interaction of the two terminal time points.
`ηmm`: The self-interaction of all intermediate points.
`η0m`: The interaction between a terminal and an intermediate point at different time separations.
`ηmn`: The interaction between two intermediate points at different time separations.
`η0e`: The interaction between the two terminal points at different time separations.
"""
struct EtaCoeffs
    η00::ComplexF64
    ηmm::ComplexF64
    η0m::Vector{ComplexF64}
    ηmn::Vector{ComplexF64}
    η0e::Vector{ComplexF64}
end

function calculate_η(ω, sd, β::Real, dt::Real, kmax::Int, classical::Bool=false, imaginary_only=false, discrete::Bool=false)
    common = common_part(ω, sd, β, classical)
    η00 = 1.0 / (2π) * Utilities.trapezoid(ω, common .* (1.0 .- exp.(-1im .* ω .* dt / 2.0)); discrete)
    ηmm = 1.0 / (2π) * Utilities.trapezoid(ω, common .* (1.0 .- exp.(-1im .* ω .* dt)); discrete)

    η0m = zeros(ComplexF64, kmax)
    ηmn = zeros(ComplexF64, kmax)
    η0e = zeros(ComplexF64, kmax)
    sin_quarter = sin.(ω * dt / 4.0)
    sin_half = sin.(ω * dt / 2.0)

    for k = 1:kmax
        η0m[k] = 2.0/π * Utilities.trapezoid(ω, common .* sin_quarter .* sin_half .* exp.(-1im .* ω .* (k-0.25) .* dt); discrete)
        η0e[k] = 2.0/π * Utilities.trapezoid(ω, common .* sin_quarter.^2 .* exp.(-1im .* ω .* (k-0.5) .* dt); discrete)
        ηmn[k] = 2.0/π * Utilities.trapezoid(ω, common .* sin_half.^2 .* exp.(-1im .* ω .* k .* dt); discrete)
    end

    imaginary_only ? EtaCoeffs(1im*imag(η00), 1im*imag(ηmm), 1im.*imag.(η0m), 1im.*imag.(ηmn), 1im.*imag.(η0e)) : EtaCoeffs(η00, ηmm, η0m, ηmn, η0e)
end

"""
    calculate_η(specdens<:SpectralDensities.AnalyticalSpectralDensity; β::Real, dt::Real, kmax::Int, classical::Bool=false, discrete::Bool=false)
Calculates the η-coefficients from an analytic spectral density and returns them as an object of the structure `EtaCoeffs`. The integrations involved are done using trapezoidal integration
"""
function calculate_η(specdens::T; β::Real, dt::Real, kmax::Int, classical::Bool=false, imaginary_only=false) where {T<:SpectralDensities.AnalyticalSpectralDensity}
    ω, sd = SpectralDensities.tabulate(specdens)
    calculate_η(ω, sd, β, dt, kmax, classical, imaginary_only, false)
end

"""
    calculate_η(specdens::T; β::Real, dt::Real, kmax::Int, classical::Bool=false, imaginary_only=false) where {T<:SpectralDensities.ContinuousSpectralDensity}
Calculates the η-coefficients from a a tabulated spectral density and returns them as an object of the structure `EtaCoeffs`.
"""
function calculate_η(specdens::T; β::Real, dt::Real, kmax::Int, classical::Bool=false, imaginary_only=false) where {T<:SpectralDensities.ContinuousSpectralDensity}
    calculate_η(specdens.ω, specdens.jw, β, dt, kmax, classical, imaginary_only, false)
end

"""
    calculate_η(specdens::T; β::Real, dt::Real, kmax::Int, classical::Bool=false, imaginary_only=false) where {T<:SpectralDensities.DiscreteOscillators}
Calculates the η-coefficients from a discretized set of harmonic modes and returns them as an object of the structure `EtaCoeffs`. The integrations involved are converted to sums over frequency modes.
"""
function calculate_η(specdens::T; β::Real, dt::Real, kmax::Int, classical::Bool=false, imaginary_only=false) where {T<:SpectralDensities.DiscreteOscillators}
    calculate_η(specdens.ω, specdens.jw, β, dt, kmax, classical, imaginary_only, true)
end

end