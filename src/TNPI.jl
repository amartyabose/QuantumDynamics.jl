module TNPI

using ITensors
using ..EtaCoefficients, ..Propagators, ..SpectralDensities, ..Blip, ..Utilities

struct TNPIArgs <: Utilities.ExtraArgs
    cutoff :: Float64
    maxdim :: Int
    method :: String
end
TNPIArgs(; cutoff=1e-8, maxdim=50, method="naive") = TNPIArgs(cutoff, maxdim, method)

function build_path_amplitude_mps(fbU, sites)
    fbUtens = ITensor(fbU, sites)
    U, V = factorize(fbUtens, (sites[1]); ortho="left", which_decomp="svd")
    ans = MPS(2)
    ans[1] = U
    ans[2] = V
    ans
end

function extend_path_amplitude_mps(pamps, fbU, sites)
    additional_part = build_path_amplitude_mps(fbU, sites)
    ans = MPS(length(pamps)+1)
    for (i, pa) in enumerate(pamps)
        ans[i] = deepcopy(pa)
    end
    ans[end-1] = ITensor(unioninds(pamps[end], additional_part[1]))
    old_last_link = linkinds(pamps)[end]
    site_ind = siteinds(pamps)[end]
    new_last_link = linkinds(additional_part)[1]
    for ol = 1:dim(old_last_link)
        for s = 1:dim(site_ind)
            for nl = 1:dim(new_last_link)
                ans[end-1][old_last_link=>ol, site_ind=>s, new_last_link=>nl] = pamps[end][old_last_link=>ol, site_ind=>s] * additional_part[1][site_ind=>s, new_last_link=>nl]
            end
        end
    end 
    ans[end] = additional_part[2]
    ans
end

function extend_path_amplitude_mps_beyond_memory(pamps, fbU, sites)
    old_sites = siteinds(pamps)[1:3]
    trace_op = ITensor(old_sites[2])
    for iv in eachindval(old_sites[2])
        trace_op[iv] = 1
    end
    tensor13 = pamps[1] * pamps[2] * pamps[3] * trace_op
    U, V = factorize(tensor13, (old_sites[1]); ortho="left", which_decomp="svd")
    additional_part = build_path_amplitude_mps(fbU, sites)
    ans = MPS(length(pamps))
    ans[1] = U
    ans[2] = V
    for (i, pa) in enumerate(pamps[4:end])
        ans[i+2] = deepcopy(pa)
    end
    ans[end-1] = ITensor(unioninds(pamps[end], additional_part[1]))
    old_last_link = linkinds(pamps)[end]
    site_ind = siteinds(pamps)[end]
    new_last_link = linkinds(additional_part)[1]
    for ol = 1:dim(old_last_link)
        for s = 1:dim(site_ind)
            for nl = 1:dim(new_last_link)
                ans[end-1][old_last_link=>ol, site_ind=>s, new_last_link=>nl] = pamps[end][old_last_link=>ol, site_ind=>s] * additional_part[1][site_ind=>s, new_last_link=>nl]
            end
        end
    end 
    ans[end] = additional_part[2]
    ans
end

function path_amplitude_to_propagator(pamps)
    ans = pamps[1]
    sinds = siteinds(pamps)
    curr_site = sinds[1]
    trace_op = ITensor(curr_site)
    for iv in eachindval(curr_site)
        trace_op[iv] = 1
    end
    for (i, ps) in enumerate(pamps[2:end-1])
        swapinds!(trace_op, [sinds[i]], [sinds[i+1]])
        ans *= ps * trace_op
    end
    noprime(ans * pamps[end])
end

function apply_contract_propagator(pamps, ifmpo)
    ans = pamps[1] * ifmpo[1]
    sinds = siteinds(pamps)
    curr_site = sinds[1]
    trace_op = ITensor(curr_site)
    for iv in eachindval(curr_site)
        trace_op[iv] = 1
    end
    for (i, ps) in enumerate(pamps[2:end-1])
        swapinds!(trace_op, [sinds[i]], [sinds[i+1]])
        ans *= ps * ifmpo[i+1] * trace_op'
    end
    noprime(ans * pamps[end] * ifmpo[end])
end

function build_ifmpo(; ηs::Vector{EtaCoefficients.EtaCoeffs}, group_Δs, Δs, sbar, sites)
    nsites = length(sites)
    sitedim = dim(sites[1])
    linkdim = length(group_Δs)
    term_ifmpo = MPO(nsites)
    cont_ifmpo = MPO(nsites)
    links = [Index(linkdim, "Link") for j in 1:nsites-1]

    tensor = ITensor(sites[1], sites[1]', links[1])
    for s = 1:sitedim
        for β = 1:linkdim
            tensor[sites[1]=>s, sites[1]'=>s, links[1]=>β] = exp(sum([-group_Δs[bn, β] * (real(η.η0e[nsites-1]) * Δs[bn, s] + 2im * imag(η.η0e[nsites-1]) * sbar[bn, s]) for (bn, η) in enumerate(ηs)]))
        end
    end
    term_ifmpo[1] = tensor
    for (i, site) in enumerate(sites[2:end-1])
        tensor = ITensor(site, site', links[i], links[i+1])
        for s = 1:sitedim
            for β = 1:linkdim
                tensor[site=>s, site'=>s, links[i]=>β, links[i+1]=>β] = exp(sum([-group_Δs[bn, β] * (real(η.η0m[nsites-1-i]) * Δs[bn, s] + 2im * imag(η.η0m[nsites-1-i]) * sbar[bn, s]) for (bn, η) in enumerate(ηs)]))
            end
        end
        term_ifmpo[i+1] = tensor
    end
    tensor = ITensor(sites[end], sites[end]', links[end])
    for s = 1:sitedim
        for β = 1:linkdim
            tensor[sites[end]=>s, sites[end]'=>s, links[end]=>β] = Δs[:, s]==group_Δs[:, β] ? exp(sum([-group_Δs[bn, β] * (real(η.η00) * Δs[bn, s] + 2im * imag(η.η00) * sbar[bn, s]) for (bn, η) in enumerate(ηs)])) : 0
        end
    end
    term_ifmpo[end] = tensor

    tensor = ITensor(sites[1], sites[1]', links[1])
    for s = 1:sitedim
        for β = 1:linkdim
            tensor[sites[1]=>s, sites[1]'=>s, links[1]=>β] = exp(sum([-group_Δs[bn, β] * (real(η.η0m[nsites-1]) * Δs[bn, s] + 2im * imag(η.η0m[nsites-1]) * sbar[bn, s]) for (bn, η) in enumerate(ηs)]))
        end
    end
    cont_ifmpo[1] = tensor
    for (i, site) in enumerate(sites[2:end-1])
        tensor = ITensor(site, site', links[i], links[i+1])
        for s = 1:sitedim
            for β = 1:linkdim
                tensor[site=>s, site'=>s, links[i]=>β, links[i+1]=>β] = exp(sum([-group_Δs[bn, β] * (real(η.ηmn[nsites-1-i]) * Δs[bn, s] + 2im * imag(η.ηmn[nsites-1-i]) * sbar[bn, s]) for (bn, η) in enumerate(ηs)]))
            end
        end
        cont_ifmpo[i+1] = tensor
    end
    tensor = ITensor(sites[end], sites[end]', links[end])
    for s = 1:sitedim
        for β = 1:linkdim
            tensor[sites[end]=>s, sites[end]'=>s, links[end]=>β] = Δs[:, s]==group_Δs[:, β] ? exp(sum([-group_Δs[bn, β] * (real(η.ηmm) * Δs[bn, s] + 2im * imag(η.ηmm) * sbar[bn, s]) for (bn, η) in enumerate(ηs)])) : 0
        end
    end
    cont_ifmpo[end] = tensor

    cont_ifmpo, term_ifmpo
end

function extend_ifmpo(; ηs::Vector{EtaCoefficients.EtaCoeffs}, group_Δs, Δs, sbar, sites, old_cont_ifmpo, old_term_ifmpo)
    nsites = length(sites)
    sitedim = dim(sites[1])
    linkdim = length(group_Δs)
    term_ifmpo = MPO(nsites)
    cont_ifmpo = MPO(nsites)
    for j = 3:nsites
        term_ifmpo[j] = swapinds(old_term_ifmpo[j-1], [sites[j-1], sites[j-1]'], [sites[j], sites[j]'])
        cont_ifmpo[j] = swapinds(old_cont_ifmpo[j-1], [sites[j-1], sites[j-1]'], [sites[j], sites[j]'])
    end

    link12 = Index(linkdim, "Link")
    tensor = ITensor(sites[1], sites[1]', link12)
    for s = 1:sitedim
        for β = 1:linkdim
            tensor[sites[1]=>s, sites[1]'=>s, link12=>β] = exp(sum([-group_Δs[bn, β] * (real(η.η0e[nsites-1]) * Δs[bn, s] + 2im * imag(η.η0e[nsites-1]) * sbar[bn, s]) for (bn, η) in enumerate(ηs)]))
        end
    end
    term_ifmpo[1] = tensor

    tensor = ITensor(sites[1], sites[1]', link12)
    for s = 1:sitedim
        for β = 1:linkdim
            tensor[sites[1]=>s, sites[1]'=>s, link12=>β] = exp(sum([-group_Δs[bn, β] * (real(η.η0m[nsites-1]) * Δs[bn, s] + 2im * imag(η.η0m[nsites-1]) * sbar[bn, s]) for (bn, η) in enumerate(ηs)]))
        end
    end
    cont_ifmpo[1] = tensor

    oldlink = linkinds(old_term_ifmpo)[1]
    tensor = ITensor(sites[2], sites[2]', link12, oldlink)
    for s = 1:sitedim
        for β = 1:linkdim
            tensor[sites[2]=>s, sites[2]'=>s, link12=>β, oldlink=>β] = exp(sum([-group_Δs[bn, β] * (real(η.η0m[nsites-2]) * Δs[bn, s] + 2im * imag(η.η0m[nsites-2]) * sbar[bn, s]) for (bn, η) in enumerate(ηs)]))
        end
    end
    term_ifmpo[2] = tensor

    oldlink = linkinds(old_cont_ifmpo)[1]
    tensor = ITensor(sites[2], sites[2]', link12, oldlink)
    for s = 1:sitedim
        for β = 1:linkdim
            tensor[sites[2]=>s, sites[2]'=>s, link12=>β, oldlink=>β] = exp(sum([-group_Δs[bn, β] * (real(η.ηmn[nsites-2]) * Δs[bn, s] + 2im * imag(η.ηmn[nsites-2]) * sbar[bn, s]) for (bn, η) in enumerate(ηs)]))
        end
    end
    cont_ifmpo[2] = tensor
    cont_ifmpo, term_ifmpo
end

function extend_ifmpo_kmax_plus_1(; ηs::Vector{EtaCoefficients.EtaCoeffs}, group_Δs, Δs, sbar, sites, old_cont_ifmpo, old_term_ifmpo)
    nsites = length(sites)
    sitedim = dim(sites[1])
    linkdim = length(group_Δs)
    term_ifmpo = MPO(nsites)
    cont_ifmpo = MPO(nsites)
    for j = 3:nsites
        term_ifmpo[j] = swapinds(old_term_ifmpo[j-1], [sites[j-1], sites[j-1]'], [sites[j], sites[j]'])
        cont_ifmpo[j] = swapinds(old_cont_ifmpo[j-1], [sites[j-1], sites[j-1]'], [sites[j], sites[j]'])
    end

    link12 = Index(1, "Link")
    tensor = ITensor(sites[1], sites[1]', link12)
    for s = 1:sitedim
        tensor[sites[1]=>s, sites[1]'=>s, link12=>1] = 1.0
    end
    term_ifmpo[1] = tensor

    tensor = ITensor(sites[1], sites[1]', link12)
    for s = 1:sitedim
        tensor[sites[1]=>s, sites[1]'=>s, link12=>1] = 1.0
    end
    cont_ifmpo[1] = tensor

    oldlink = linkinds(old_term_ifmpo)[1]
    tensor = ITensor(sites[2], sites[2]', link12, oldlink)
    for s = 1:sitedim
        for β = 1:linkdim
            tensor[sites[2]=>s, sites[2]'=>s, link12=>1, oldlink=>β] = exp(sum([-group_Δs[bn, β] * (real(η.η0m[nsites-2]) * Δs[bn, s] + 2im * imag(η.η0m[nsites-2]) * sbar[bn, s]) for (bn, η) in enumerate(ηs)]))
        end
    end
    term_ifmpo[2] = tensor

    oldlink = linkinds(old_cont_ifmpo)[1]
    tensor = ITensor(sites[2], sites[2]', link12, oldlink)
    for s = 1:sitedim
        for β = 1:linkdim
            tensor[sites[2]=>s, sites[2]'=>s, link12=>1, oldlink=>β] = exp(sum([-group_Δs[bn, β] * (real(η.ηmn[nsites-2]) * Δs[bn, s] + 2im * imag(η.ηmn[nsites-2]) * sbar[bn, s]) for (bn, η) in enumerate(ηs)]))
        end
    end
    cont_ifmpo[2] = tensor
    cont_ifmpo, term_ifmpo
end

function extend_ifmpo_beyond_memory(; sites, old_cont_ifmpo, old_term_ifmpo, count)
    nsites = length(old_cont_ifmpo)
    term_ifmpo = MPO(nsites)
    cont_ifmpo = MPO(nsites)
    for j = 2:nsites
        term_ifmpo[j] = swapinds(old_term_ifmpo[j], [sites[count+j-1], sites[count+j-1]'], [sites[count+j], sites[count+j]'])
        cont_ifmpo[j] = swapinds(old_cont_ifmpo[j], [sites[count+j-1], sites[count+j-1]'], [sites[count+j], sites[count+j]'])
    end
    term_ifmpo[1] = old_term_ifmpo[1]
    cont_ifmpo[1] = old_cont_ifmpo[1]
    cont_ifmpo, term_ifmpo
end

function convert_to_matrix(prop_tens, sinit, sterm)
    prop = zeros(ComplexF64, dim(sterm), dim(sinit))
    for j = 1:dim(sterm)
        for k = 1:dim(sinit)
            prop[j, k] = prop_tens[sinit=>k, sterm=>j]
        end
    end
    prop
end

function build_augmented_propagator(; fbU::Array{ComplexF64, 3}, Jw::Vector{T}, β::Real, dt::Real, ntimes::Int, kmax::Int, svec=[1.0 -1.0], reference_prop=false, extraargs::TNPIArgs, end_prop=false) where {T<:SpectralDensities.SpectralDensity}
    @assert length(Jw) == size(svec, 1)
    ηs = [EtaCoefficients.calculate_η(jw; β, dt, kmax=min(kmax, ntimes), imaginary_only=reference_prop) for jw in Jw]
    sdim2 = size(fbU, 2)
    _, group_Δs, sbar, Δs = Blip.setup_simulation(svec)

    sites = siteinds(sdim2, ntimes+1)

    fbU1 = deepcopy(fbU[1, :, :])
    infl = zeros(ComplexF64, sdim2)
    for (b, η) in enumerate(ηs)
        infl .+= -Δs[b, :] .* (real(η.η00) .* Δs[b, :] .+ 2im .* imag(η.η00) .* sbar[b, :])
    end
    infl .= exp.(infl)
    for j in eachcol(fbU1)
        j .*= infl
    end
    pamps = build_path_amplitude_mps(fbU1, sites[1:2])

    U0e = zeros(ComplexF64, ntimes, sdim2, sdim2)
    cont_ifmpo, term_ifmpo = build_ifmpo(; ηs, group_Δs, Δs, sbar, sites=sites[1:2])
    U0e[1, :, :] .= convert_to_matrix(apply_contract_propagator(pamps, term_ifmpo), sites[1], sites[2])
    for j = 2:min(kmax, ntimes)
        pamps_cont = apply(cont_ifmpo, pamps; cutoff=extraargs.cutoff, maxdim=extraargs.maxdim, method=extraargs.method)
        pamps = extend_path_amplitude_mps(pamps_cont, fbU[j, :, :], sites[j:j+1])
        cont_ifmpo, term_ifmpo = extend_ifmpo(; ηs, group_Δs, Δs, sbar, sites=sites[1:j+1], old_cont_ifmpo=cont_ifmpo, old_term_ifmpo=term_ifmpo)
        U0e[j, :, :] .= convert_to_matrix(apply_contract_propagator(pamps, term_ifmpo), sites[1], sites[j+1])
    end

    if ntimes > kmax
        pamps_cont = apply(cont_ifmpo, pamps; cutoff=extraargs.cutoff, maxdim=extraargs.maxdim, method=extraargs.method)
        pamps = extend_path_amplitude_mps(pamps_cont, fbU[kmax+1, :, :], sites[kmax+1:kmax+2])
        cont_ifmpo, term_ifmpo = extend_ifmpo_kmax_plus_1(; ηs, group_Δs, Δs, sbar, sites=sites[1:kmax+2], old_cont_ifmpo=cont_ifmpo, old_term_ifmpo=term_ifmpo)
        U0e[kmax+1, :, :] .= convert_to_matrix(apply_contract_propagator(pamps, term_ifmpo), sites[1], sites[kmax+2])

        count = 1
        for j = kmax+2 : ntimes
            pamps_cont = apply(cont_ifmpo, pamps; cutoff=extraargs.cutoff, maxdim=extraargs.maxdim, method=extraargs.method)
            pamps = extend_path_amplitude_mps_beyond_memory(pamps_cont, fbU[j, :, :], sites[j:j+1])
            cont_ifmpo, term_ifmpo = extend_ifmpo_beyond_memory(; sites=sites[1:j+1], old_cont_ifmpo=cont_ifmpo, old_term_ifmpo=term_ifmpo, count)
            U0e[j, :, :] .= convert_to_matrix(apply_contract_propagator(pamps, term_ifmpo), sites[1], sites[j+1])
            count += 1
        end
    end

    U0e
end

"""
    propagate(; fbU::Matrix{ComplexF64}, Jw::Vector{T}, β::Real, ρ0, dt::Real, ntimes::Int, kmax::Int, extraargs::QuAPIArgs=QuAPIArgs(), svec=[1.0 -1.0], reference_prop=false, verbose::Bool=false) where {T<:SpectralDensities.SpectralDensity}

Given a time-series of system forward-backward propagators, `fbU`, the spectral densities describing the solvent, `Jw`, and an inverse temperature, this uses TNPI to propagate the input initial reduced density matrix, ρ0, with a time-step of `dt` for `ntimes` time steps. A non-Markovian memory of `kmax` steps is used in this simulation. The i^th bath, described by `Jw[i]`, interacts with the system through the diagonal operator with the values of `svec[j,:]`.

    `ρ0`: initial reduced density matrix
    `fbU`: time-series of forward-backward propagators
    `Jw`: array of spectral densities
    `svec`: diagonal elements of system operators through which the corresponding baths interact. QuAPI currently only works for baths with diagonal coupling to the system.

    `dt`: time-step for recording the density matrices
    `ntimes`: number of time steps of simulation
    `kmax`: number of steps within memory
    `extraargs`: extra arguments for the TNPI algorithm. Contains the `cutoff` threshold for SVD filtration, the maximum bond dimension, `maxdim`, and the `method` of applying an MPO to an MPS.
"""
function propagate(; fbU::Array{ComplexF64, 3}, Jw::Vector{T}, β::Real, ρ0::Matrix{ComplexF64}, dt::Real, ntimes::Int, kmax::Int, extraargs::TNPIArgs = TNPIArgs(), svec=[1.0 -1.0], reference_prop=false, verbose::Bool=false) where {T<:SpectralDensities.SpectralDensity}
    U0e = build_augmented_propagator(; fbU, Jw, β, dt, ntimes, kmax, extraargs, svec, reference_prop)
    Utilities.apply_propagator(; propagators=U0e, ρ0, ntimes, dt)
end

end