module GQME

using HDF5
using ..Propagators, ..TTM, ..SpectralDensities, ..Utilities

"""
    propagate(; Hamiltonian::AbstractMatrix{<:Complex}, Jw::Vector{T}, β::Real, ρ0::AbstractMatrix{<:Complex}, dt::Real, ntimes::Int, rmax::Int, kmax::Union{Int,Nothing}=nothing, path_integral_routine, extraargs::Utilities.ExtraArgs, svec=[1.0 -1.0], reference_prop=false, verbose::Bool=false, L::Union{Nothing,Vector{Matrix{ComplexF64}}}=nothing, output::Union{Nothing,HDF5.Group}=nothing) where {T<:SpectralDensities.SpectralDensity}
Propagate an initial density matrix `ρ0` under a `Hamiltonian`, environments described by the spectral densities `Jw` held at an inverse temperature `β`. The memory kernel is obtained using the transfer tensor method.
"""
function propagate(; Hamiltonian::AbstractMatrix{<:Complex}, Jw::Vector{T}, β::Real, ρ0::AbstractMatrix{<:Complex}, dt::Real, ntimes::Int, rmax::Int, kmax::Union{Int,Nothing}=nothing, path_integral_routine, extraargs::Utilities.ExtraArgs, svec=[1.0 -1.0], reference_prop=false, verbose::Bool=false, L::Union{Nothing,Vector{Matrix{ComplexF64}}}=nothing, output::Union{Nothing,HDF5.Group}=nothing) where {T<:SpectralDensities.SpectralDensity}
    fbU = Propagators.calculate_bare_propagators(; Hamiltonian, dt, ntimes=rmax)
    U0e_within_r = path_integral_routine(; fbU, Jw, β, dt, ntimes=rmax, kmax, extraargs, svec, verbose, reference_prop, output)
    T0e = TTM.get_Ts(U0e_within_r)
    K = TTM.get_memory_kernel(T0e, fbU[1, :, :], dt)

    propagate_with_memory_kernel(; K, fbU, ρ0, dt, ntimes, L)
end

"""
    propagate_with_memory_kernel(; K::AbstractArray{<:Complex,3}, fbU::AbstractMatrix{<:Complex}, ρ0::AbstractMatrix{<:Complex}, dt::Real, ntimes::Int, L::Union{Nothing,Vector{Matrix{ComplexF64}}}=nothing)
Given a memory kernel `K`, and the bare forward-backward propagator `fbU`, propagate the initial density matrix `ρ0`. Can additionally involve empirical losses through Lindblad jump operators `L`.
"""
function propagate_with_memory_kernel(; K::AbstractArray{<:Complex,3}, fbU::AbstractMatrix{<:Complex}, ρ0::AbstractMatrix{<:Complex}, dt::Real, ntimes::Int, L::Union{Nothing,Vector{Matrix{ComplexF64}}}=nothing)
    sdim = size(ρ0, 1)
    sdim2 = sdim^2
    ρs = zeros(eltype(ρ0), ntimes + 1, sdim, sdim)
    @inbounds ρs[1, :, :] = ρ0
    ρvec = Utilities.density_matrix_to_vector(ρ0)
    ρsvec = zeros(eltype(ρs), ntimes + 1, sdim2)
    ρsvec[1, :] = ρvec

    rmax = size(K, 1)

    for j = 2:ntimes+1
        dρ = zero(ρvec)
        for i = 1:min(rmax, j - 1)
            dρ += K[i, :, :] * ρsvec[j-i, :]
        end
        dρ *= dt^2
        ρsvec[j, :] = fbU * ρsvec[j-1, :] + dρ
        ρs[j, :, :] = Utilities.density_matrix_vector_to_matrix(ρsvec[j, :])
        if !isnothing(L)
            for l in L
                ρs[j, :, :] .+= (l * ρs[j-1, :, :] * l' .- 0.5 .* l' * l * ρs[j-1, :, :] - 0.5 .* ρs[j-1, :, :] * l' * l) * dt
            end
            ρsvec[j, :] .= Utilities.density_matrix_to_vector(ρs[j, :, :])
        end
    end
    0:dt:ntimes*dt, ρs
end

end