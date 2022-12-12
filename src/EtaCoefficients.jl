module EtaCoefficients

using ..SpectralDensities

function trapezoid(x, y, discrete::Bool)
    if discrete
        return sum(y)
    end
    sum = zero(y[1])
    for (a, b) in zip(y[2:end], y)
        sum += a + b
    end
    sum / 2 * (x[2] - x[1])
end

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

struct EtaCoeffs
    η00::ComplexF64
    ηmm::ComplexF64
    η0m::Vector{ComplexF64}
    ηmn::Vector{ComplexF64}
    η0e::Vector{ComplexF64}
end

function calculate_η(specdens::SpectralDensities.SpectralDensity; β::Real, dt::Real, kmax::Int, classical::Bool=false, discrete::Bool=false)
    ω, sd = SpectralDensities.tabulate(specdens)
    common = common_part(ω, sd, β, classical)
    η00 = 1.0 / (2π) * trapezoid(ω, common .* (1.0 .- exp.(-1im .* ω .* dt / 2.0)), discrete)
    ηmm = 1.0 / (2π) * trapezoid(ω, common .* (1.0 .- exp.(-1im .* ω .* dt)), discrete)

    η0m = zeros(ComplexF64, kmax)
    ηmn = zeros(ComplexF64, kmax)
    η0e = zeros(ComplexF64, kmax)
    sin_quarter = sin.(ω * dt / 4.0)
    sin_half = sin.(ω * dt / 2.0)

    for k = 1:kmax
        η0m[k] = 2.0/π * trapezoid(ω, common .* sin_quarter .* sin_half .* exp.(-1im .* ω .* (k-0.25) .* dt), discrete)
        η0e[k] = 2.0/π * trapezoid(ω, common .* sin_quarter.^2 .* exp.(-1im .* ω .* (k-0.5) .* dt), discrete)
        ηmn[k] = 2.0/π * trapezoid(ω, common .* sin_half.^2 .* exp.(-1im .* ω .* k .* dt), discrete)
    end

    EtaCoeffs(η00, ηmm, η0m, ηmn, η0e)
end

end