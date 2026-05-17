module Polaron

using ..SpectralDensities
using LinearAlgebra

function full_polaron_transform(; Hamiltonian::AbstractMatrix{<:Number}, Jw::AbstractVector{SpectralDensities.SpectralDensity}, svec=[1.0 -1.0], β::Real)
    H = copy(Hamiltonian)
    N = size(Hamiltonian, 1)
    Γ = zeros(N, N)
    for (i, J) in enumerate(Jw)
        λ = SpectralDensities.reorganization_energy(J)
        H .-= diagm(svec[i, :].^2 .* λ) 
	polaron_factor = SpectralDensities.polaron_shielding(J, β)
        for j=1:N, k=j+1:N
	    Δs = svec[i,j] - svec[i,k]
            Γ[j, k] += Δs^2 * polaron_factor
            Γ[k, j] += Δs^2 * polaron_factor
        end
    end
    for j = 1:N, k = j+1:N
        H[j, k] *= exp(-1.0 / 2 * Γ[j, k])
        H[k, j] *= exp(-1.0 / 2 * Γ[k, j])
    end
    H
end

end
