module HEOM

using SparseArrays
using DifferentialEquations
using ..SpectralDensities, ..Utilities

nvec2string(xs) = string(["$(x)" for x in xs]...)

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
    for n = 1:size(ρ, 3)
        dρ[:,:,n] .= -1im * (params.H * ρ[:,:,n] - ρ[:,:,n] * params.H)
        dρ[:,:,n] .-= sum(params.nveclist[n] .* params.γ) .* ρ[:,:,n]
        for (i, (jw, co)) in enumerate(zip(params.Jw, params.coupl))
            # dρ[:,:,n] .-= (2 * jw.λ / (jw.Δs^2 * jw.γ * params.β) - real(sum(params.c[i,:] ./ params.γ[i,:]))) .* Utilities.commutator(co, Utilities.commutator(co, ρ[:,:,n]))
            dρ[:,:,n] .-= params.Δk .* Utilities.commutator(co, Utilities.commutator(co, ρ[:,:,n]))
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
                    dρ[:,:,n] .+= -1im * sqrt(nvec[m, k] / abs(params.c[m, k])) * (params.c[m, k] * co[m] * ρ[:,:,nminuslocs[m, k]] - conj(params.c[m,k]) * ρ[:,:,nminuslocs[m,k]] * co[m])
                end
            end
            dρ[:,:,n] .+= -1im * (co[m] * ρplus - ρplus * co[m])
        end
    end
end

function propagate(; Hamiltonian::Matrix{ComplexF64}, ρ0::Matrix{ComplexF64}, β::Real, Jw::Vector{SpectralDensities.DrudeLorentzCutoff}, sys_ops::Vector{Matrix{ComplexF64}}, num_modes::Int, Lmax::Int, dt::Real, ntimes::Int, extraargs::Utilities.DiffEqArgs=Utilities.DiffEqArgs())
    γ = zeros(length(Jw), num_modes+1)
    c = zeros(ComplexF64, length(Jw), num_modes+1)
    Δk = zeros(length(Jw))
    for (i,jw) in enumerate(Jw)
        γj, cj, δk = SpectralDensities.matsubara_decomposition(jw, num_modes, β)
        @show γj, cj, 2 * jw.λ / (jw.Δs^2 * jw.γ * β) - real(sum(cj ./ γj)), δk
        γ[i,:] .= γj
        c[i,:] .= cj
        Δk[i] = δk
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