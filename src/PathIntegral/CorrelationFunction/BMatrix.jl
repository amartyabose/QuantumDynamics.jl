module BMatrix

using ....SpectralDensities, ..ComplexPISetup, ....Utilities

function compute_B!(B::AbstractMatrix{<:Complex}, ω::AbstractVector{<:AbstractFloat}, common_part::AbstractVector{<:AbstractFloat}, tarr::AbstractVector{<:Complex}, npoints::Int, β::AbstractFloat)
    @inbounds begin
        for k = 1:npoints
            for kp = 1:k-1
                B[k, kp] = 4 / π * Utilities.trapezoid(ω, common_part .* cos.(ω .* (tarr[k+1] + tarr[k] - tarr[kp+1] - tarr[kp] + 1im * β) ./ 2) .* sin.(ω .* (tarr[k+1] - tarr[k]) ./ 2) .* sin.(ω .* (tarr[kp+1] - tarr[kp]) ./ 2))
                B[kp, k] = B[k, kp]
            end
            B[k, k] = 2 / π * Utilities.trapezoid(ω, common_part .* sin.(ω .* (tarr[k+1] - tarr[k] + 1im * β) ./ 2) .* sin.(ω .* (tarr[k+1] - tarr[k]) ./ 2))
        end
    end
end

common_part(ω, j, β) = j ./ (ω.^2 .* sinh.(ω .* β ./ 2))

function get_B_matrix(ω::AbstractVector{<:AbstractFloat}, j::AbstractVector{<:AbstractFloat}, β::Real, N, tarr)
    # common_part = j ./ (ω .^ 2 .* sinh.(ω .* β ./ 2))
    comm = common_part(ω, j, β)
    npoints = 2N + 2
    B = zeros(ComplexF64, npoints, npoints)
    if Utilities.trapezoid(ω, comm) ≈ 0.0
        return B
    end
    compute_B!(B, ω, comm, tarr, npoints, β)
    B
end

function get_B_matrix(J::SpectralDensities.SpectralDensity, β::Real, N, tarr)
    ω, j = SpectralDensities.tabulate(J, false)
    get_B_matrix(ω, j, β, N, tarr)
end

end