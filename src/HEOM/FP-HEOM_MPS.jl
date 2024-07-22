module FPHEOM

using ITensors, ITensorMPS
using ..SpectralDensities, ..Utilities

# """
#     propagate_MPO(; Hamiltonian::AbstractMatrix{ComplexF64}, ρ0::AbstractMatrix{ComplexF64}, β::Real, Jw::AbstractVector{SpectralDensities.DrudeLorentz}, sys_ops::Vector{Matrix{ComplexF64}}, num_modes::Int, Lmax::Int, dt::Real, ntimes::Int, threshold::Float64=0.0, external_fields::Union{Nothing,Vector{Utilities.ExternalField}}=nothing, extraargs::Utilities.DiffEqArgs=Utilities.DiffEqArgs())
# 
# Uses HEOM to propagate the initial reduced density matrix, `ρ0`, under the given `Hamiltonian`, and set of spectral densities, `Jw`, interacting with the system through `sys_ops`.
# 
# `ρ0`: initial reduced density matrix
# `Hamiltonian`: system Hamiltonian
# `external_fields`: either `nothing` or a vector of external time-dependent fields
# `Jw`: array of spectral densities
# `sys_ops`: system operators through which the corresponding baths interact
# 
# `num_modes`: number of Matsubara modes to be considered
# `Lmax`: cutoff for maximum number of levels
# `dt`: time-step for recording the density matrices
# `ntimes`: number of time steps of simulation
# `threshold`: filtration threshold
# `extraargs`: extra arguments for the differential equation solver
# """
function propagate_MPO(; Hamiltonian::AbstractMatrix{ComplexF64}, ρ0::AbstractMatrix{ComplexF64}, β::Real, Jw::AbstractVector{SpectralDensities.DrudeLorentz}, sys_ops::Vector{Matrix{ComplexF64}}, num_modes::Int, Lmax::Int, dt::Real, ntimes::Int, threshold::Float64=0.0, scaled::Bool=true, external_fields::Union{Nothing,Vector{Utilities.ExternalField}}=nothing, extraargs::Utilities.DiffEqArgs=Utilities.DiffEqArgs())
    nbaths = length(Jw)
    num_sites = 2 + nbaths * (num_modes + 1)
    start_end = zeros(Int64, nbaths, 2)
    for j = 1:nbaths
        start_end[j, 1] = j + 2 + (num_modes + 1) * (j - 1)
        start_end[j, 2] = j + 2 + (num_modes + 1) * j - 1
    end

    γ = zeros(nbaths, num_modes + 1)
    c = zeros(ComplexF64, nbaths, num_modes + 1)
    Δk = zeros(nbaths)
    for (i, jw) in enumerate(Jw)
        γj, cj = SpectralDensities.matsubara_decomposition(jw, num_modes, β)
        @inbounds γ[i, :] .= γj
        @inbounds c[i, :] .= cj
        Δk[i] = (2 * jw.λ / (jw.Δs^2 * jw.γ * β) - real(sum(cj ./ γj))) # residual sum used to truncate the hierarchy
    end

    sdim = size(Hamiltonian, 1)
    sites = Vector{Index}(undef, num_sites)
    sites[1] = Index(sdim, "i")
    sites[2] = Index(sdim, "j")
    for j = 3:num_sites
        sites[j] = Index(Lmax + 1, "n$(j-2)")
    end

    H = ITensor(Hamiltonian, [sites[1]', sites[1]])
    L = -1im * Utilities.calculate_Liouvillian(H, [sites[1], sites[2]])
    LL, LR = factorize(L, [sites[1], sites[1]']; tags="Link,l=1")

    linds_lmpo = [Index(1, "Link,l=$(j)") for j = 2:num_sites]
    linds_lmpo = [commonind(LL, LR); linds_lmpo]
    lmpo_vec = Vector{ITensor}(undef, num_sites)
    lmpo_vec[1] = LL
    lmpo_vec[2] = LR * onehot(linds_lmpo[2] => 1)
    for j = 3:num_sites-1
        lmpo_vec[j] = onehot(linds_lmpo[j-1] => 1) * delta(sites[j], sites[j]') * onehot(linds_lmpo[j] => 1)
    end
    lmpo_vec[end] = onehot(linds_lmpo[end] => 1) * delta(sites[end], sites[end]')
    lmpo = MPO(lmpo_vec)

    lmpo_vec[1] = delta(sites[1], sites[1]') * onehot(linds_lmpo[1] => 1)
    lmpo_vec[2] = onehot(linds_lmpo[1] => 1) * delta(sites[2], sites[2]') * onehot(linds_lmpo[2] => 1)
    for j = 3:num_sites
        tmp = ITensor(sites[j], sites[j]')
        for s = 1:dim(sites[j])
            tmp[sites[j]=>s, sites[j]'=>s] = (s - 1) * γ[j-2]
        end
        lmp_vec[j] = onehot(linds_lmpo[j-1] => 1) * tmp * onehot(linds_lmpo[j] => 1)
    end

    sites, MPO(lmpo)
end

end