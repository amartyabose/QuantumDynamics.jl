module HEOM

using SparseArrays
using DifferentialEquations
using ..SpectralDensities, ..Utilities

function get_vecs(len::Int, L::Int)
    len==1 && return [L]
    ans = []
    for j = 0:L
        rest = get_vecs(len-1, L-j)
        curr = [cat(j, r; dims=1) for r in rest]
        append!(ans, curr)
    end
    ans
end

function setup_simulation(num_baths::Int, num_modes::Int, Lmax::Int)
    nveclist = []
    len = num_baths * (num_modes + 1)
    num = 1
    for L = 0:Lmax
        vecs = get_vecs(len, L)
        for v in vecs
            push!(nveclist, reshape(v, num_baths, num_modes+1))
            num += 1
        end
    end
    npluslocs = zeros(Int, num_baths, num_modes+1, length(nveclist))
    nminuslocs = zeros(Int, num_baths, num_modes+1, length(nveclist))
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
    H :: Matrix{ComplexF64}
    Jw :: Vector{SpectralDensities.DrudeLorentzCutoff}
    coupl :: Vector{Matrix{ComplexF64}}
    nveclist
    npluslocs
    nminuslocs
    γ
    c
    Δk
    β
end

function HEOM_RHS!(dρ, ρ, params, t)
    @inbounds begin
        for n = 1:size(ρ, 3)
            dρ[:,:,n] .= -1im * Utilities.commutator(params.H, ρ[:,:,n])
            dρ[:,:,n] .-= sum(params.nveclist[n] .* params.γ) .* ρ[:,:,n]
            for (Δk, co) in zip(params.Δk, params.coupl)
                dρ[:,:,n] .-= Δk .* Utilities.commutator(co, Utilities.commutator(co, ρ[:,:,n]))
            end
        end

        for n = 1:size(ρ, 3)
            nvec = params.nveclist[n]
            npluslocs = params.npluslocs[:,:,n]
            nminuslocs = params.nminuslocs[:,:,n]
            for (m, co) in enumerate(params.coupl)
                ρplus = zeros(ComplexF64, size(params.H, 1), size(params.H, 2))
                for k = 1:size(npluslocs, 2)
                    if npluslocs[m, k] > 0
                        ρplus .+= sqrt((nvec[m, k] + 1) * abs(params.c[m, k])) * ρ[:,:,npluslocs[m, k]]
                    end
                    if nminuslocs[m, k] > 0
                        dρ[:,:,n] .+= -1im * sqrt(nvec[m, k] / abs(params.c[m, k])) * (params.c[m, k] * co * ρ[:,:,nminuslocs[m, k]] - conj(params.c[m,k]) * ρ[:,:,nminuslocs[m,k]] * co)
                    end
                end
                dρ[:,:,n] .+= -1im * Utilities.commutator(co, ρplus)
            end
        end
    end
end

"""
    propagate(; Hamiltonian::Matrix{ComplexF64}, ρ0::Matrix{ComplexF64}, β::Real, Jw::Vector{SpectralDensities.DrudeLorentzCutoff}, sys_ops::Vector{Matrix{ComplexF64}}, num_modes::Int, Lmax::Int, dt::Real, ntimes::Int, extraargs::Utilities.DiffEqArgs=Utilities.DiffEqArgs())

Uses HEOM to propagate the initial reduced density matrix, `ρ0`, under the given `Hamiltonian`, and set of spectral densities, `Jw`, interacting with the system through `sys_ops`.

    `ρ0`: initial reduced density matrix
    `Hamiltonian`: system Hamiltonian
    `Jw`: array of spectral densities
    `sys_ops`: system operators through which the corresponding baths interact

    `num_modes`: number of Matsubara modes to be considered
    `Lmax`: cutoff for maximum number of levels
    `dt`: time-step for recording the density matrices
    `ntimes`: number of time steps of simulation
    `extraargs`: extra arguments for the differential equation solver
"""
function propagate(; Hamiltonian::Matrix{ComplexF64}, ρ0::Matrix{ComplexF64}, β::Real, Jw::Vector{SpectralDensities.DrudeLorentzCutoff}, sys_ops::Vector{Matrix{ComplexF64}}, num_modes::Int, Lmax::Int, dt::Real, ntimes::Int, extraargs::Utilities.DiffEqArgs=Utilities.DiffEqArgs())
    γ = zeros(length(Jw), num_modes+1)
    c = zeros(ComplexF64, length(Jw), num_modes+1)
    Δk = zeros(length(Jw))
    for (i,jw) in enumerate(Jw)
        γj, cj = SpectralDensities.matsubara_decomposition(jw, num_modes, β)
        @inbounds γ[i,:] .= γj
        @inbounds c[i,:] .= cj
        Δk[i] = (2 * jw.λ / (jw.Δs^2 * jw.γ * β) - real(sum(cj ./ γj)))
    end
    nveclist, npluslocs, nminuslocs = setup_simulation(length(Jw), num_modes, Lmax)
    params = HEOMParams(Hamiltonian, Jw, sys_ops, nveclist, npluslocs, nminuslocs, γ, c, Δk, β)
    tspan = (0.0, dt*ntimes)
    sdim = size(ρ0, 1)
    ρ0_expanded = zeros(ComplexF64, sdim, sdim, length(nveclist))
    ρ0_expanded[:,:,1] .= ρ0
    prob = ODEProblem{true}(HEOM_RHS!, ρ0_expanded, tspan, params)
    sol = solve(prob, extraargs.solver, reltol=extraargs.reltol, abstol=extraargs.abstol, saveat=dt)
    ρs = zeros(ComplexF64, length(sol.t), sdim, sdim)
    for j = 1:length(sol.t)
        @inbounds ρs[j, :, :] .= sol.u[j][:,:,1]
    end
    sol.t, ρs
end

end