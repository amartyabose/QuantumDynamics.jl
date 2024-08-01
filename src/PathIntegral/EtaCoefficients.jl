module EtaCoefficients

using FLoops
using ..SpectralDensities, ..Utilities

function common_part(ω, sd, β, classical)
    elem_type = eltype(sd)
    if isnan(β) && !classical
        ans = zero(ω)
        for (i, freq) in enumerate(ω)
            if freq > 0
                ans[i] = sd ./ ω .^ 2 .* 2 * one(elem_type)
            else
                ans[i] = 0
            end
        end
        return ans
    elseif !classical
        return sd ./ ω .^ 2 .* (2 * one(elem_type) ./ (1 .- exp.(-ω .* β)))
    elseif classical
        return sd ./ ω .^ 2 .* (2 * one(elem_type) ./ (ω .* β) .+ 1)
    end
end

abstract type IFCoeffs end

"""
EtaCoefficients holds the various discretized η-coefficients required for a QuAPI-based simulation. These are the minimum number of coefficients required, stored using time-translational symmetry wherever possible.

The values are stored as follows:
- `η00`: The self-interaction of the two terminal time points.
- `ηmm`: The self-interaction of all intermediate points.
- `η0m`: The interaction between a terminal and an intermediate point at different time separations.
- `ηmn`: The interaction between two intermediate points at different time separations.
- `η0e`: The interaction between the two terminal points at different time separations.
"""
struct EtaCoeffs <: IFCoeffs
    η00::ComplexF64
    ηmm::ComplexF64
    η0m::Vector{ComplexF64}
    ηmn::Vector{ComplexF64}
    η0e::Vector{ComplexF64}
end

function calculate_η(ω, sd, β::Real, dt::Real, kmax::Int, classical::Bool=false, imaginary_only=false, discrete::Bool=false)
    common = common_part(ω, sd, β, classical)
    elem_type = eltype(sd)
    η00 = 1 / elem_type(2π) * Utilities.trapezoid(ω, common .* (1 .- exp.(-1im .* ω .* dt / 2)); discrete, exec=FLoops.ThreadedEx())
    ηmm = 1 / elem_type(2π) * Utilities.trapezoid(ω, common .* (1 .- exp.(-1im .* ω .* dt)); discrete, exec=FLoops.ThreadedEx())

    η0m = zeros(Complex{elem_type}, kmax)
    ηmn = zeros(Complex{elem_type}, kmax)
    η0e = zeros(Complex{elem_type}, kmax)
    sin_quarter = sin.(ω * dt / 4)
    sin_half = sin.(ω * dt / 2)

    for k = 1:kmax
        η0m[k] = 2 / elem_type(π) * Utilities.trapezoid(ω, common .* sin_quarter .* sin_half .* exp.(-1im .* ω .* (k - 1 / 4) .* dt); discrete, exec=FLoops.ThreadedEx())
        η0e[k] = 2 / elem_type(π) * Utilities.trapezoid(ω, common .* sin_quarter .^ 2 .* exp.(-1im .* ω .* (k - 1 / 2) .* dt); discrete, exec=FLoops.ThreadedEx())
        ηmn[k] = 2 / elem_type(π) * Utilities.trapezoid(ω, common .* sin_half .^ 2 .* exp.(-1im .* ω .* k .* dt); discrete, exec=FLoops.ThreadedEx())
    end

    imaginary_only ? EtaCoeffs(1im * imag(η00), 1im * imag(ηmm), 1im .* imag.(η0m), 1im .* imag.(ηmn), 1im .* imag.(η0e)) : EtaCoeffs(η00, ηmm, η0m, ηmn, η0e)
end

"""
    calculate_η(specdens::SpectralDensities.ContinuousSpectralDensity; β::Real, dt::Real, kmax::Int, imaginary_only=false)
Calculates the η-coefficients from an analytic spectral density and returns them as an object of the structure `EtaCoeffs`. The integrations involved are done using trapezoidal integration
"""
function calculate_η(specdens::SpectralDensities.ContinuousSpectralDensity; β::Real, dt::Real, kmax::Int, imaginary_only=false)
    ω, sd = SpectralDensities.tabulate(specdens)
    calculate_η(ω, sd, β, dt, kmax, specdens.classical, imaginary_only, false)
end

"""
    calculate_η(specdens::SpectralDensity.DiscreteOscillators; β::Real, dt::Real, kmax::Int, classical::Bool=false, imaginary_only=false)
Calculates the η-coefficients from a discretized set of harmonic modes and returns them as an object of the structure `EtaCoeffs`. The integrations involved are converted to sums over frequency modes.
"""
function calculate_η(specdens::SpectralDensities.DiscreteOscillators; β::Real, dt::Real, kmax::Int, imaginary_only=false)
    ω, sd = SpectralDensities.tabulate(specdens)
    calculate_η(ω, sd, β, dt, kmax, specdens.classical, imaginary_only, true)
end

"""
ZetaCoefficients holds the various discretized ζ-coefficients required for a QCPI-HBR-based simulation. These are the minimum number of coefficients required, stored using time-translational symmetry wherever possible.

The values are stored as follows:
- `ζ00`: The self-interaction of the two terminal time points.
- `ζmm`: The self-interaction of all intermediate points.
- `ζ0m`: The interaction between a terminal and an intermediate point at different time separations.
- `ζmn`: The interaction between two intermediate points at different time separations.
- `ζ0e`: The interaction between the two terminal points at different time separations.
"""
struct ZetaCoeffs <: IFCoeffs
    ζ00::Float64
    ζmm::Vector{Float64}
    ζ0m::Vector{Float64}
    ζme::Matrix{Float64}
    ζmn::Matrix{Float64}
    ζ0e::Vector{Float64}
end

function calculate_ζ(ω, sd, dt::Real, kmax::Int, discrete::Bool=false)
    sin_dt_4 = sin.(ω * dt / 4)
    sin_dt_2 = sin.(ω * dt / 2)
    sin_dt = sin.(ω * dt)
    one_minus_cos_dt_2 = 1 .- cos.(ω * dt / 2)
    sin_dt_4_sin_dt_2 = sin_dt_4 .* sin_dt_2

    common = 1.0 / π * sd ./ ω .^ 2
    elem_type = eltype(sd)
    ζ00 = Utilities.trapezoid(ω, common .* (0.5 * ω * dt .- sin_dt_2); discrete, exec=FLoops.ThreadedEx())

    ζmm = zeros(elem_type, 2)
    ζmm[1] = Utilities.trapezoid(ω, common .* (0.5 * ω * dt .+ sin_dt_2 .- sin_dt); discrete, exec=FLoops.ThreadedEx())
    ζmm[2] = ζ00

    ζ0e = zeros(elem_type, kmax)
    ζ0m = zeros(elem_type, kmax)
    ζme = zeros(elem_type, kmax, 2)
    ζmn = zeros(elem_type, kmax, 2)
    for k = 1:kmax
        ζ0e[k] = 2 * Utilities.trapezoid(ω, common .* one_minus_cos_dt_2 .* sin.((k - 0.5) * dt * ω))
        ζ0m[k] = 4 * Utilities.trapezoid(ω, common .* sin_dt_4_sin_dt_2 .* sin.((k - 0.25) * dt * ω))
        ζme[k, 1] = 2 * Utilities.trapezoid(ω, common .* one_minus_cos_dt_2 .* sin.(k * dt * ω))
        ζme[k, 2] = 2 * Utilities.trapezoid(ω, common .* one_minus_cos_dt_2 .* sin.((k - 0.5) * dt * ω))
        ζmn[k, 1] = 4 * Utilities.trapezoid(ω, common .* sin_dt_4_sin_dt_2 .* sin.((k + 0.25) * dt * ω))
        ζmn[k, 2] = ζ0m[k]
    end

    ZetaCoeffs(ζ00, ζmm, ζ0m, ζme, ζmn, ζ0e)
end

"""
    calculate_ζ(specdens::SpectralDensities.ContinuousSpectralDensity; dt::Real, kmax::Int)
Calculates the ζ-coefficients from an analytic spectral density and returns them as an object of the structure `EtaCoeffs`. The integrations involved are done using trapezoidal integration
"""
function calculate_ζ(specdens::SpectralDensities.ContinuousSpectralDensity; dt::Real, kmax::Int)
    ω, sd = SpectralDensities.tabulate(specdens)
    calculate_ζ(ω, sd, dt, kmax, false)
end

"""
    calculate_ζ(specdens::SpectralDensity.DiscreteOscillators; dt::Real, kmax::Int)
Calculates the ζ-coefficients from a discretized set of harmonic modes and returns them as an object of the structure `EtaCoeffs`. The integrations involved are converted to sums over frequency modes.
"""
function calculate_ζ(specdens::SpectralDensities.DiscreteOscillators; dt::Real, kmax::Int)
    ω, sd = SpectralDensities.tabulate(specdens)
    calculate_ζ(ω, sd, dt, kmax, true)
end

end