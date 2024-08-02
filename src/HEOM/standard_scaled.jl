module HEOM

using OrdinaryDiffEq
using ..SpectralDensities, ..Utilities

"""
    get_vecs(len::Int, L::Int)

Get a vector of vectors of length `len`, where the sum is L.
"""
function get_vecs(len::Int, L::Int)
    len == 1 && return [L]
    ans = Vector{Vector{typeof(L)}}()
    for j = 0:L
        rest = get_vecs(len - 1, L - j)
        curr = [cat(j, r; dims=1) for r in rest]
        append!(ans, curr)
    end
    ans
end

"""
    setup_simulation(num_baths::Int, num_modes::Int, Lmax::Int)

Sets up the simulation parameters for a problem with `num_baths` baths, `num_modes` extra matsubara modes, and a hierarchy `Lmax` levels deep.

Returns a tuple of:
- `nveclist`: List of the possible subscripts, `n`, in HEOM. Each element in the list is a represented as a matrix. Every row corresponds to a bath.
- `npluslocs[b,m,l]`: Given the `l`th nvector, returns the location of the nvector if the `b`th bath's `m`th Matsubara mode is increased by one.
- `nminuslocs[b,m,l]`: Given the `l`th nvector, returns the location of the nvector if the `b`th bath's `m`th Matsubara mode is decreased by one.
"""
function setup_simulation(num_baths::Int, num_modes::Int, Lmax::Int)
    nveclist = Vector{Matrix{typeof(Lmax)}}()
    len = num_baths * (num_modes + 1) # for each bath there are num_modes + 1 matsubara modes in total. (+1 coming from the zero mode)
    num = 1
    for L = 0:Lmax
        vecs = get_vecs(len, L)
        for v in vecs
            push!(nveclist, reshape(v, num_baths, num_modes + 1))
            num += 1
        end
    end
    npluslocs = zeros(Int, num_baths, num_modes + 1, length(nveclist))
    nminuslocs = zeros(Int, num_baths, num_modes + 1, length(nveclist))
    for (j, nvec) in enumerate(nveclist)
        for m = 1:num_baths
            for k = 1:num_modes+1
                nvecnew = deepcopy(nvec)
                nvecnew[m, k] += 1
                for (i, vec) in enumerate(nveclist)
                    if nvecnew == vec
                        npluslocs[m, k, j] = i
                        break
                    end
                end
                nvecnew = deepcopy(nvec)
                nvecnew[m, k] -= 1
                for (i, vec) in enumerate(nveclist)
                    if nvecnew == vec
                        nminuslocs[m, k, j] = i
                        break
                    end
                end
            end
        end
    end
    nveclist, npluslocs, nminuslocs
end

struct HEOMParams
    H::Matrix{ComplexF64}
    external_fields::Union{Nothing,Vector{Utilities.ExternalField}}
    Jw::Vector{SpectralDensities.DrudeLorentz}
    coupl::Vector{Matrix{ComplexF64}}
    nveclist
    npluslocs
    nminuslocs
    γ
    c
    Δk
    β
    threshold
end

function scaled_HEOM_RHS!(dρ, ρ, params, t)
    @inbounds begin
        H = deepcopy(params.H)
        if !isnothing(params.external_fields)
            for ef in params.external_fields
                H .+= ef.V(t) * ef.coupling_op
            end
        end
        for n in axes(ρ, 3)
            if maximum(abs.(ρ[:, :, n])) ≤ params.threshold
                dρ[:, :, n] .= 0
            else
                dρ[:, :, n] .= -1im * Utilities.nh_commutator(H, ρ[:, :, n])
                dρ[:, :, n] .-= sum(params.nveclist[n] .* params.γ) .* ρ[:, :, n]
                for (Δk, co) in zip(params.Δk, params.coupl)
                    dρ[:, :, n] .-= Δk .* Utilities.commutator(co, Utilities.commutator(co, ρ[:, :, n]))
                end

                nvec = params.nveclist[n]
                npluslocs = params.npluslocs[:, :, n]
                nminuslocs = params.nminuslocs[:, :, n]
                for (m, co) in enumerate(params.coupl)
                    ρplus = zeros(ComplexF64, size(params.H, 1), size(params.H, 2))
                    for k in axes(npluslocs, 2)
                        if npluslocs[m, k] > 0
                            ρplus .+= sqrt((nvec[m, k] + 1) * abs(params.c[m, k])) * ρ[:, :, npluslocs[m, k]]
                        end
                        if nminuslocs[m, k] > 0
                            dρ[:, :, n] .+= -1im * sqrt(nvec[m, k] / abs(params.c[m, k])) * (params.c[m, k] * co * ρ[:, :, nminuslocs[m, k]] .- conj(params.c[m, k]) * ρ[:, :, nminuslocs[m, k]] * co)
                        end
                    end
                    dρ[:, :, n] .+= -1im * Utilities.commutator(co, ρplus)
                end
            end
        end
    end
    nothing
end

function unscaled_HEOM_RHS!(dρ, ρ, params, t)
    @inbounds begin
        H = deepcopy(params.H)
        if !isnothing(params.external_fields)
            for ef in params.external_fields
                H .+= ef.V(t) * ef.coupling_op
            end
        end
        for n in axes(ρ, 3)
            dρ[:, :, n] .= -1im * Utilities.nh_commutator(H, ρ[:, :, n])
            dρ[:, :, n] .-= sum(params.nveclist[n] .* params.γ) .* ρ[:, :, n]
            for (Δk, co) in zip(params.Δk, params.coupl)
                dρ[:, :, n] .-= Δk .* Utilities.commutator(co, Utilities.commutator(co, ρ[:, :, n]))
            end
        end

        for n in axes(ρ, 3)
            nvec = params.nveclist[n]
            npluslocs = params.npluslocs[:, :, n]
            nminuslocs = params.nminuslocs[:, :, n]
            for (m, co) in enumerate(params.coupl)
                ρplus = zeros(ComplexF64, size(params.H, 1), size(params.H, 2))
                for k in axes(npluslocs, 2)
                    if npluslocs[m, k] > 0
                        ρplus .+= ρ[:, :, npluslocs[m, k]]
                    end
                    if nminuslocs[m, k] > 0
                        dρ[:, :, n] .+= -1im * nvec[m, k] * (params.c[m, k] * co * ρ[:, :, nminuslocs[m, k]] - conj(params.c[m, k]) * ρ[:, :, nminuslocs[m, k]] * co)
                    end
                end
                dρ[:, :, n] .+= -1im * Utilities.commutator(co, ρplus)
            end
        end
    end
    nothing
end

"""
    propagate(; Hamiltonian::AbstractMatrix{ComplexF64}, ρ0::AbstractMatrix{ComplexF64}, β::Real, Jw::AbstractVector{SpectralDensities.DrudeLorentz}, sys_ops::Vector{Matrix{ComplexF64}}, num_modes::Int, Lmax::Int, dt::Real, ntimes::Int, threshold::Float64=0.0, external_fields::Union{Nothing,Vector{Utilities.ExternalField}}=nothing, extraargs::Utilities.DiffEqArgs=Utilities.DiffEqArgs())

Uses HEOM to propagate the initial reduced density matrix, `ρ0`, under the given `Hamiltonian`, and set of spectral densities, `Jw`, interacting with the system through `sys_ops`.

- `ρ0`: initial reduced density matrix
- `Hamiltonian`: system Hamiltonian
- `external_fields`: either `nothing` or a vector of external time-dependent fields
- `Jw`: array of spectral densities
- `sys_ops`: system operators through which the corresponding baths interact

- `num_modes`: number of Matsubara modes to be considered
- `Lmax`: cutoff for maximum number of levels
- `dt`: time-step for recording the density matrices
- `ntimes`: number of time steps of simulation
- `threshold`: filtration threshold
- `extraargs`: extra arguments for the differential equation solver
"""
function propagate(; Hamiltonian::AbstractMatrix{ComplexF64}, ρ0::AbstractMatrix{ComplexF64}, β::Real, Jw::AbstractVector{SpectralDensities.DrudeLorentz}, sys_ops::Vector{Matrix{ComplexF64}}, num_modes::Int, Lmax::Int, dt::Real, ntimes::Int, threshold::Float64=0.0, scaled::Bool=true, external_fields::Union{Nothing,Vector{Utilities.ExternalField}}=nothing, extraargs::Utilities.DiffEqArgs=Utilities.DiffEqArgs())
    γ = zeros(length(Jw), num_modes + 1)
    c = zeros(ComplexF64, length(Jw), num_modes + 1)
    Δk = zeros(length(Jw))
    for (i, jw) in enumerate(Jw)
        γj, cj = SpectralDensities.matsubara_decomposition(jw, num_modes, β)
        @inbounds γ[i, :] .= γj
        @inbounds c[i, :] .= cj
        Δk[i] = (2 * jw.λ / (jw.Δs^2 * jw.γ * β) - real(sum(cj ./ γj))) # residual sum used to truncate the hierarchy
    end
    nveclist, npluslocs, nminuslocs = setup_simulation(length(Jw), num_modes, Lmax)

    params = HEOMParams(Hamiltonian, external_fields, Jw, sys_ops, nveclist, npluslocs, nminuslocs, γ, c, Δk, β, threshold)
    tspan = (0.0, dt * ntimes)
    sdim = size(ρ0, 1)
    ρ0_expanded = zeros(ComplexF64, sdim, sdim, length(nveclist))
    for j in axes(ρ0_expanded, 3)
        ρ0_expanded[:, :, j] .= ρ0
    end
    prob = scaled ? ODEProblem{true}(scaled_HEOM_RHS!, ρ0_expanded, tspan, params) : ODEProblem{true}(unscaled_HEOM_RHS!, ρ0_expanded, tspan, params)
    sol = solve(prob, extraargs.solver, reltol=extraargs.reltol, abstol=extraargs.abstol, saveat=dt)
    ρs = zeros(ComplexF64, length(sol.t), sdim, sdim)
    for j = 1:length(sol.t)
        @inbounds ρs[j, :, :] .= sol.u[j][:, :, 1]
    end
    sol.t, ρs
end

end