module BMatrix

using ....SpectralDensities, ..ComplexPISetup, ....Utilities

function compute_B!(B::AbstractMatrix{<:Complex}, ω::AbstractVector{<:AbstractFloat}, jw::AbstractVector{<:AbstractFloat}, tarr::AbstractVector{<:Complex}, npoints::Int, β::AbstractFloat)
    @inbounds begin
        common_part = jw ./ ω.^2
        coth_term = 1 .+ 2 ./ expm1.(β .* ω)
        tdiff = tarr[2:end] .- tarr[1:end-1]
        for k = 1:npoints
            for kp = 1:k-1
                B[k, kp] = 4 / π * Utilities.trapezoid(ω, common_part .* (cos.(ω .* (tarr[k+1] + tarr[k] - tarr[kp+1] - tarr[kp]) ./ 2) .* coth_term - 1im * sin.(ω .* (tarr[k+1] + tarr[k] - tarr[kp+1] - tarr[kp]) ./ 2)) .* sin.(ω .* tdiff[k] ./ 2) .* sin.(ω .* tdiff[kp] ./ 2))
                B[kp, k] = B[k, kp]
            end
            B[k, k] = 2 / π * Utilities.trapezoid(ω, common_part .* (sin.(ω .* tdiff[k] ./ 2) .* coth_term + 1im * cos.(ω .* tdiff[k] ./ 2)) .* sin.(ω .* tdiff[k] ./ 2))
        end
    end
end

function get_B_matrix(ω::AbstractVector{<:AbstractFloat}, j::AbstractVector{<:AbstractFloat}, β::Real, N, tarr)
    npoints = 2N + 2
    B = zeros(ComplexF64, npoints, npoints)
    compute_B!(B, ω, j, tarr, npoints, β)
    B
end

function get_B_matrix(J::SpectralDensities.SpectralDensity, β::Real, N, tarr)
    ω, j = SpectralDensities.tabulate(J, false)
    get_B_matrix(ω, j, β, N, tarr)
end

end
