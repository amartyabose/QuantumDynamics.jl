module HEOM

using OrdinaryDiffEq
using FLoops
using ..SpectralDensities, ..Utilities

const references = """
- Y. Tanimura and R. Kubo, Time Evolution of a Quantum System in Contact with a Nearly Gaussian-Markoffian Noise Bath, Journal of the Physical Society of Japan 58, 101 (1989).
- Q. Shi, L. Chen, G. Nan, R.-X. Xu, and Y. Yan, Efficient hierarchical Liouville space propagator to quantum dissipative dynamics, J. Chem. Phys. 130, 084105 (2009)."""

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
    nveclist = Vector{Matrix{Int}}()
    len = num_baths * (num_modes + 1)
    for L = 0:Lmax
        vecs = HEOM.get_vecs(len, L)
        for v in vecs
            push!(nveclist, reshape(v, num_baths, num_modes + 1))
        end
    end
    Nh = length(nveclist)
    index = Dict{NTuple{len,Int}, Int}()
    for (i, v) in enumerate(nveclist)
        index[Tuple(vec(v))] = i
    end

    npluslocs  = zeros(Int, num_baths, num_modes + 1, Nh)
    nminuslocs = zeros(Int, num_baths, num_modes + 1, Nh)
    for (j, nvec) in enumerate(nveclist)
        base_key = Tuple(vec(nvec))
        for m = 1:num_baths
            for k = 1:(num_modes + 1)
                nvec_plus = copy(nvec)
                nvec_plus[m, k] += 1
                npluslocs[m, k, j] = get(index, Tuple(vec(nvec_plus)), 0)
                if nvec[m, k] > 0
                    nvec_minus = copy(nvec)
                    nvec_minus[m, k] -= 1
                    nminuslocs[m, k, j] = get(index, Tuple(vec(nvec_minus)), 0)
                else
                    nminuslocs[m, k, j] = 0
                end
            end
        end
    end

    nveclist, npluslocs, nminuslocs
end

struct HEOMParams
    H::Matrix{ComplexF64}
    L::Union{Nothing,Vector{Matrix{ComplexF64}}}
    LdagL::Union{Nothing, Vector{Matrix{ComplexF64}}}
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
    decay::Vector{Float64}
    workspace::Matrix{ComplexF64}
    tmp1::Matrix{ComplexF64}
end

function scaled_HEOM_RHS!(dρ, ρ, params, t)
    @inbounds begin
        params.tmp1 .= params.H
        if !isnothing(params.external_fields)
            for ef in params.external_fields
                params.tmp1 .+= ef.V(t) * ef.coupling_op
            end
        end
        for n in axes(ρ, 3)
#             @init begin
#                 ρplus = similar(params.workspace)
#             end

            if maximum(abs, ρ[:, :, n]) ≤ params.threshold
                dρ[:, :, n] .= 0
            else
                dρ[:, :, n] .= -1im * Utilities.nh_commutator(params.tmp1, ρ[:, :, n])
                if !isnothing(params.L)
                    for (L, LdagL) in zip(params.L, params.LdagL)
                        dρ[:, :, n] .+= L * ρ[:, :, n] * L' .- 0.5 .* LdagL * ρ[:, :, n] .- 0.5 .* ρ[:, :, n] * LdagL
                    end
                end
                @. dρ[:, :, n] -= params.decay[n] * ρ[:, :, n]
                for (Δk, co) in zip(params.Δk, params.coupl)
                    dρ[:, :, n] .-= Δk .* Utilities.commutator(co, Utilities.commutator(co, ρ[:, :, n]))
                end

                @views begin
                    nvec = params.nveclist[n]
                    npluslocs = params.npluslocs[:, :, n]
                    nminuslocs = params.nminuslocs[:, :, n]
                    ρplus = params.workspace
                end
                for (m, co) in enumerate(params.coupl)
                    fill!(ρplus, 0.0)
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
            if !isnothing(params.L)
                for L in params.L
                    dρ[:, :, n] .+= L * ρ[:, :, n] * L' .- 0.5 .* L' * L * ρ[:, :, n] .- 0.5 .* ρ[:, :, n] * L' * L
                end
            end
            @. dρ[:, :, n] -= params.decay * ρ[:, :, n]
            for (Δk, co) in zip(params.Δk, params.coupl)
                dρ[:, :, n] .-= Δk .* Utilities.commutator(co, Utilities.commutator(co, ρ[:, :, n]))
            end
        end

        for n in axes(ρ, 3)
            nvec = params.nveclist[n]
            npluslocs = params.npluslocs[:, :, n]
            nminuslocs = params.nminuslocs[:, :, n]
            for (m, co) in enumerate(params.coupl)
                ρplus = params.workspace
                fill!(ρplus, 0.0)
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
    propagate(; Hamiltonian::AbstractMatrix{ComplexF64}, ρ0::AbstractMatrix{ComplexF64}, β::Real, Jw::AbstractVector{SpectralDensities.SpectralDensity}, sys_ops::Vector{Matrix{ComplexF64}}, num_modes::Int, Lmax::Int, dt::Real, ntimes::Int, threshold::Float64=0.0, scaled::Bool=true, L::Union{Nothing,Vector{Matrix{ComplexF64}}}=nothing, external_fields::Union{Nothing,Vector{Utilities.ExternalField}}=nothing, extraargs::Utilities.DiffEqArgs=Utilities.DiffEqArgs())
Uses HEOM to propagate the initial reduced density matrix, `ρ0`, under the given `Hamiltonian`, and set of spectral densities, `Jw`, interacting with the system through `sys_ops`.

- `ρ0`: initial reduced density matrix
- `Hamiltonian`: system Hamiltonian
- `external_fields`: either `nothing` or a vector of external time-dependent fields
- `Jw`: array of spectral densities
- `sys_ops`: system operators through which the corresponding baths interact
- `L`: vector of Lindblad jump operators

- `num_modes`: number of Matsubara modes to be considered
- `Lmax`: cutoff for maximum number of levels
- `dt`: time-step for recording the density matrices
- `ntimes`: number of time steps of simulation
- `threshold`: filtration threshold
- `extraargs`: extra arguments for the differential equation solver
"""
function propagate(; Hamiltonian::AbstractMatrix{ComplexF64}, ρ0::AbstractMatrix{ComplexF64}, β::Real, Jw::AbstractVector{SpectralDensities.SpectralDensity}, sys_ops::Vector{Matrix{ComplexF64}}, num_modes::Int, Lmax::Int, dt::Real, ntimes::Int, threshold::Float64=0.0, scaled::Bool=true, L::Union{Nothing,Vector{Matrix{ComplexF64}}}=nothing, external_fields::Union{Nothing,Vector{Utilities.ExternalField}}=nothing, extraargs::Utilities.DiffEqArgs=Utilities.DiffEqArgs(), decomposition::String, verbose=false)
    γ = zeros(length(Jw), num_modes + 1)
    c = zeros(ComplexF64, length(Jw), num_modes + 1)
    Δk = zeros(length(Jw))
    Δk_imag = zeros(length(Jw))
    for (i, jw) in enumerate(Jw)
        @assert typeof(jw) == SpectralDensities.DrudeLorentz "HEOM has only been implemented for the Drude-Lorentz spectral density."
        γj, cj = decomposition == "matsubara" ? SpectralDensities.matsubara_decomposition(jw, num_modes, β) : SpectralDensities.pade_decomposition(jw, num_modes, β)
        @inbounds γ[i, :] .= γj
        @inbounds c[i, :] .= cj
        tmp = sum(cj ./ γj)
        Δk[i] = (2 * jw.λ / (jw.Δs^2 * jw.γ * β) - real(tmp)) # residual sum used to truncate the hierarchy
        Δk_imag[i] = (-jw.λ - imag(tmp))
        verbose && @info "Decomposed bath number $i."
    end
    nveclist, npluslocs, nminuslocs = setup_simulation(length(Jw), num_modes, Lmax)
    verbose && @info "Setup complete. Starting run"

    H = deepcopy(Hamiltonian)
    for (Δi, co) in zip(Δk_imag, sys_ops)
        H .+= Δi * (co * co)
    end

    Nh = length(nveclist)
    sdim = size(ρ0, 1)
    workspace = zeros(ComplexF64, sdim, sdim)
    tmp1 = zeros(ComplexF64, sdim, sdim)
    tmp2 = zeros(ComplexF64, sdim, sdim)

    LdagL = if isnothing(L)
        nothing
    else
        [l' * l for l in L]
    end
    decay = zeros(Float64, length(nveclist))
    for (i, nvec) in enumerate(nveclist)
        decay[i] = sum(nvec .* γ)
    end
    params = HEOMParams(H, L, LdagL, external_fields, Jw, sys_ops, nveclist, npluslocs, nminuslocs, γ, c, Δk, β, threshold, decay, workspace, tmp1)
    tspan = (0.0, dt * ntimes)
    sdim = size(ρ0, 1)
    ρ0_expanded = zeros(ComplexF64, sdim, sdim, length(nveclist))
    for j in axes(ρ0_expanded, 3)
        ρ0_expanded[:, :, j] .= ρ0
    end
    prob = scaled ? ODEProblem{true}(scaled_HEOM_RHS!, ρ0_expanded, tspan, params) : ODEProblem{true}(unscaled_HEOM_RHS!, ρ0_expanded, tspan, params)
    sol = solve(prob, extraargs.solver, reltol=extraargs.reltol, abstol=extraargs.abstol, saveat=dt, progress=verbose)
    ρs = zeros(ComplexF64, length(sol.t), sdim, sdim)
    for j = 1:length(sol.t)
        @inbounds ρs[j, :, :] .= sol.u[j][:, :, 1]
    end
    sol.t, ρs
end

end
