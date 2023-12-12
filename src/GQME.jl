module GQME

using HDF5
using ..Propagators, ..TTM, ..SpectralDensities, ..Utilities

function propagate(; Hamiltonian::AbstractMatrix{<:Complex}, Jw::Vector{T}, β::Real, ρ0::AbstractMatrix{<:Complex}, dt::Real, ntimes::Int, rmax::Int, kmax::Union{Int,Nothing}=nothing, path_integral_routine, extraargs::Utilities.ExtraArgs, svec=[1.0 -1.0], reference_prop=false, verbose::Bool=false, external_fields::Union{Nothing,Vector{Utilities.ExternalField}}=nothing, output::Union{Nothing,HDF5.Group}=nothing) where {T<:SpectralDensities.SpectralDensity}
    fbU = Propagators.calculate_bare_propagators(; Hamiltonian, dt, ntimes=rmax)
    U0e_within_r = path_integral_routine(; fbU, Jw, β, dt, ntimes=rmax, kmax, extraargs, svec, verbose, reference_prop, output)
    T0e = TTM.get_Ts(U0e_within_r)
    K = TTM.get_memory_kernel(T0e, fbU, dt)

    sdim = size(ρ0, 1)
    sdim2 = sdim^2
    ρs = zeros(eltype(ρ0), ntimes + 1, sdim, sdim)
    @inbounds ρs[1, :, :] = ρ0
    ρvec = collect(Iterators.flatten(transpose(ρ0)))
    ρsvec = zeros(eltype(ρs), ntimes + 1, sdim2)
    ρsvec[1, :] = ρvec

    rmax = size(T0e, 1)

    fbU = Propagators.calculate_bare_propagators(; Hamiltonian, dt, ntimes, external_fields)
    for j = 2:ntimes+1
        dρ = zero(ρvec)
        for i = 1:min(rmax, j - 1)
            dρ += K[i, :, :] * ρsvec[j-i, :]
        end
        dρ *= dt^2
        ρsvec[j, :] = fbU[j-1, :, :] * ρsvec[j-1, :] + dρ
        ρs[j, :, :] = transpose(reshape(ρsvec[j, :], (sdim, sdim)))
    end
    0:dt:ntimes*dt, ρs
end

end