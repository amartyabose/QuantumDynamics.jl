module TEMPO

using LinearAlgebra
using HDF5
using FLoops
using ITensors, ITensorMPS
using ..EtaCoefficients, ..Propagators, ..SpectralDensities, ..Blip, ..Utilities, ..TTM

const references = """
- Strathearn, A.; Kirton, P.; Kilda, D.; Keeling, J.; Lovett, B. W. Efficient Non-Markovian Quantum Dynamics Using Time-Evolving Matrix Product Operators. Nature Communications 2018, 9, 3322. https://doi.org/10.1038/s41467-018-05617-3.
- Bose, A.; Walters, P. L. A Tensor Network Representation of Path Integrals: Implementation and Analysis. arXiv pre-print server arXiv:2106.12523 2021."""

const TEMPOArgs = Utilities.TensorNetworkArgs

function build_ifmpo(; ηs::Vector{EtaCoefficients.EtaCoeffs}, group_Δs, Δs, sbar, sites)
    nsites = length(sites)
    sitedim = dim(sites[1])
    linkdim = size(group_Δs, 2)
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
            tensor[sites[end]=>s, sites[end]'=>s, links[end]=>β] = Δs[:, s] == group_Δs[:, β] ? exp(sum([-group_Δs[bn, β] * (real(η.η00) * Δs[bn, s] + 2im * imag(η.η00) * sbar[bn, s]) for (bn, η) in enumerate(ηs)])) : 0
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
            tensor[sites[end]=>s, sites[end]'=>s, links[end]=>β] = Δs[:, s] == group_Δs[:, β] ? exp(sum([-group_Δs[bn, β] * (real(η.ηmm) * Δs[bn, s] + 2im * imag(η.ηmm) * sbar[bn, s]) for (bn, η) in enumerate(ηs)])) : 0
        end
    end
    cont_ifmpo[end] = tensor

    cont_ifmpo, term_ifmpo
end

function extend_ifmpo(; ηs::Vector{EtaCoefficients.EtaCoeffs}, group_Δs, Δs, sbar, sites, old_cont_ifmpo, old_term_ifmpo)
    nsites = length(sites)
    sitedim = dim(sites[1])
    linkdim = size(group_Δs, 2)
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
    linkdim = size(group_Δs, 2)
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

"""
    build_augmented_propagator(; fbU::Matrix{ComplexF64}, Jw::Vector{T}, β::Real, dt::Real, ntimes::Int, kmax::Union{Int, Nothing}=nothing, extraargs::TEMPOArgs=TEMPOArgs(), svec=[1.0 -1.0], reference_prop=false, verbose::Bool=false) where {T<:SpectralDensities.SpectralDensity}
Builds the propagators, augmented with the influence of the harmonic baths defined by the spectral densities `Jw`,  upto `ntimes` time-steps using the **TEMPO scheme**. If `kmax` is specified, the full memory simulation is only done for `kmax` steps, else it is done for all `ntimes` steps. The j^th bath, described by `Jw[j]`, interacts with the system through the diagonal operator with the values of `svec[j,:]`.

Relevant references:
$(references)
"""
function build_augmented_propagator(; fbU::Array{<:Complex,3}, Jw::Vector{T}, β::Real, dt::Real, ntimes::Int, kmax::Union{Int,Nothing}=nothing, svec=[1.0 -1.0], reference_prop=false, extraargs::TEMPOArgs=TEMPOArgs(), verbose::Bool=false, output::Union{Nothing,HDF5.Group}=nothing, from_TTM::Bool=false, exec=ThreadedEx()) where {T<:SpectralDensities.SpectralDensity}
    @assert isnothing(kmax) || kmax > 1
    @assert length(Jw) == size(svec, 1)
    nmem = isnothing(kmax) ? ntimes : min(kmax, ntimes)
    ηs = [EtaCoefficients.calculate_η(jw; β, dt, kmax=nmem, imaginary_only=reference_prop) for jw in Jw]
    sdim2 = size(fbU, 2)
    _, _, group_Δs, sbar, Δs = Blip.setup_simulation(svec)

    sites = siteinds(sdim2, ntimes + 1)

    fbU1 = deepcopy(fbU[1, :, :])
    infl = zeros(ComplexF64, sdim2)
    for (b, η) in enumerate(ηs)
        infl .+= -Δs[b, :] .* (real(η.η00) .* Δs[b, :] .+ 2im .* imag(η.η00) .* sbar[b, :])
    end
    infl .= exp.(infl)
    for j in eachcol(fbU1)
        j .*= infl
    end
    pamps = Utilities.build_path_amplitude_mps(fbU1, sites[1:2])

    if verbose
        @info "Starting propagation within memory"
    end
    _, time_taken, memory_allocated, gc_time, _ = @timed begin
        U0e = zeros(ComplexF64, ntimes, sdim2, sdim2)
        T0e = zero(U0e)
        cont_ifmpo, term_ifmpo = build_ifmpo(; ηs, group_Δs, Δs, sbar, sites=sites[1:2])
        U0e[1, :, :] .= Utilities.convert_ITensor_to_matrix(Utilities.apply_contract_propagator(pamps, term_ifmpo), sites[1], sites[2])
        TTM.update_Ts!(T0e, U0e, 1)
    end
    ldims = linkdims(pamps)
    maxldim = maximum(ldims)
    avgldim = sum(ldims) / length(ldims)
    if verbose
        @info "Step = 1; max bond dimension = $(maxldim); avg bond dimension = $(round(avgldim; digits=3)); norm(Tmat) = $(round(norm(T0e[1, :, :]); digits=3)); time = $(round(time_taken; digits=3)) sec; memory allocated = $(round(memory_allocated / 1e6; digits=3)) GB; gc time = $(round(gc_time; digits=3)) sec"
    end
    if !isnothing(output)
        if !from_TTM
            Utilities.check_or_insert_value(output, "U0e", U0e)
        end
        Utilities.check_or_insert_value(output, "T0e", T0e)
        Utilities.check_or_insert_value(output, "maxbonddim", zeros(Int64, ntimes))
        Utilities.check_or_insert_value(output, "avgbonddim", zeros(Float64, ntimes))
        Utilities.check_or_insert_value(output, "time_taken", zeros(Float64, ntimes))
        output["U0e"][1, :, :] = U0e[1, :, :]
        output["T0e"][1, :, :] = T0e[1, :, :]
        output["time_taken"][1] = time_taken
        output["maxbonddim"][1] = maxldim
        output["avgbonddim"][1] = avgldim
        flush(output)
    end

    for j = 2:nmem
        _, time_taken, memory_allocated, gc_time, _ = @timed begin
            pamps = extraargs.algorithm!="fit" ? Utilities.extend_path_amplitude_mps(apply(cont_ifmpo, pamps; cutoff=extraargs.cutoff, maxdim=extraargs.maxdim, alg=extraargs.algorithm), fbU[j, :, :], sites[j:j+1]) : Utilities.extend_path_amplitude_mps(apply(cont_ifmpo, pamps; cutoff=extraargs.cutoff, maxdim=extraargs.maxdim, alg=extraargs.algorithm, nsweeps=1), fbU[j, :, :], sites[j:j+1])
            cont_ifmpo, term_ifmpo = extend_ifmpo(; ηs, group_Δs, Δs, sbar, sites=sites[1:j+1], old_cont_ifmpo=cont_ifmpo, old_term_ifmpo=term_ifmpo)
            U0e[j, :, :] .= Utilities.convert_ITensor_to_matrix(Utilities.apply_contract_propagator(pamps, term_ifmpo), sites[1], sites[j+1])
            TTM.update_Ts!(T0e, U0e, j)
            GC.gc()
        end
        ldims = linkdims(pamps)
        maxldim = maximum(ldims)
        avgldim = sum(ldims) / length(ldims)
        if verbose
            @info "Step = $(j); max bond dimension = $(maxldim); avg bond dimension = $(round(avgldim; digits=3)); norm(Tmat) = $(round(norm(T0e[j, :, :]); digits=3)); time = $(round(time_taken; digits=3)) sec; memory allocated = $(round(memory_allocated / 1e6; digits=3)) GB; gc time = $(round(gc_time; digits=3)) sec"
        end
        if !isnothing(output)
            output["U0e"][j, :, :] = U0e[j, :, :]
            output["time_taken"][j] = time_taken
            output["maxbonddim"][j] = maxldim
            output["avgbonddim"][j] = avgldim
            flush(output)
        end
    end

    if !isnothing(kmax) && ntimes > kmax
        if verbose
            @info "Starting iteration"
        end
        _, time_taken, memory_allocated, gc_time, _ = @timed begin
            pamps = extraargs.algorithm!="fit" ? Utilities.extend_path_amplitude_mps(apply(cont_ifmpo, pamps; cutoff=extraargs.cutoff, maxdim=extraargs.maxdim, alg=extraargs.algorithm), fbU[kmax+1, :, :], sites[kmax+1:kmax+2]) : Utilities.extend_path_amplitude_mps(apply(cont_ifmpo, pamps; cutoff=extraargs.cutoff, maxdim=extraargs.maxdim, alg=extraargs.algorithm, nsweeps=1), fbU[kmax+1, :, :], sites[kmax+1:kmax+2])
            cont_ifmpo, term_ifmpo = extend_ifmpo_kmax_plus_1(; ηs, group_Δs, Δs, sbar, sites=sites[1:kmax+2], old_cont_ifmpo=cont_ifmpo, old_term_ifmpo=term_ifmpo)
            U0e[kmax+1, :, :] .= Utilities.convert_ITensor_to_matrix(Utilities.apply_contract_propagator(pamps, term_ifmpo), sites[1], sites[kmax+2])
            TTM.update_Ts!(T0e, U0e, kmax+1)
            GC.gc()
        end
        ldims = linkdims(pamps)
        maxldim = maximum(ldims)
        avgldim = sum(ldims) / length(ldims)
        if verbose
            @info "Step = $(kmax+1); max bond dimension = $(maxldim); avg bond dimension = $(round(avgldim; digits=3)); norm(Tmat) = $(round(norm(T0e[kmax+1, :, :]); digits=3)); time = $(round(time_taken; digits=3)) sec; memory allocated = $(round(memory_allocated / 1e6; digits=3)) GB; gc time = $(round(gc_time; digits=3)) sec"
        end
        if !isnothing(output)
            output["U0e"][kmax+1, :, :] = U0e[kmax+1, :, :]
            output["T0e"][kmax+1, :, :] = T0e[kmax+1, :, :]
            output["time_taken"][kmax+1] = time_taken
            output["maxbonddim"][kmax+1] = maxldim
            output["avgbonddim"][kmax+1] = avgldim
            flush(output)
        end

        count = 1
        for j = kmax+2:ntimes
            _, time_taken, memory_allocated, gc_time, _ = @timed begin
                pamps = extraargs.algorithm!="fit" ? Utilities.extend_path_amplitude_mps_beyond_memory(apply(cont_ifmpo, pamps; cutoff=extraargs.cutoff, maxdim=extraargs.maxdim, alg=extraargs.algorithm), fbU[j, :, :], sites[j:j+1]) : Utilities.extend_path_amplitude_mps_beyond_memory(apply(cont_ifmpo, pamps; cutoff=extraargs.cutoff, maxdim=extraargs.maxdim, alg=extraargs.algorithm, nsweeps=1), fbU[j, :, :], sites[j:j+1])
                cont_ifmpo, term_ifmpo = extend_ifmpo_beyond_memory(; sites=sites[1:j+1], old_cont_ifmpo=cont_ifmpo, old_term_ifmpo=term_ifmpo, count)
                U0e[j, :, :] .= Utilities.convert_ITensor_to_matrix(Utilities.apply_contract_propagator(pamps, term_ifmpo), sites[1], sites[j+1])
                TTM.update_Ts!(T0e, U0e, j)
                GC.gc()
            end
            ldims = linkdims(pamps)
            maxldim = maximum(ldims)
            avgldim = sum(ldims) / length(ldims)
            if verbose
                @info "Step = $(j); max bond dimension = $(maxldim); avg bond dimension = $(round(avgldim; digits=3)); norm(Tmat) = $(round(norm(T0e[j, :, :]); digits=3)); time = $(round(time_taken; digits=3)) sec; memory allocated = $(round(memory_allocated / 1e6; digits=3)) GB; gc time = $(round(gc_time; digits=3)) sec"
            end
            if !isnothing(output)
                output["U0e"][j, :, :] = U0e[j, :, :]
                output["T0e"][j, :, :] = T0e[j, :, :]
                output["time_taken"][j] = time_taken
                output["maxbonddim"][j] = maxldim
                output["avgbonddim"][j] = avgldim
                flush(output)
            end
            count += 1
        end
    end

    U0e
end

"""
    propagate(; fbU::Matrix{ComplexF64}, Jw::Vector{T}, β::Real, ρ0, dt::Real, ntimes::Int, kmax::Int, extraargs::TEMPOArgs=TEMPOArgs(), svec=[1.0 -1.0], reference_prop=false, verbose::Bool=false) where {T<:SpectralDensities.SpectralDensity}

Given a time-series of system forward-backward propagators, `fbU`, the spectral densities describing the solvent, `Jw`, and an inverse temperature, this uses TEMPO to propagate the input initial reduced density matrix, ρ0, with a time-step of `dt` for `ntimes` time steps. A non-Markovian memory of `kmax` steps is used in this simulation. The i^th bath, described by `Jw[i]`, interacts with the system through the diagonal operator with the values of `svec[j,:]`.

Relevant references:
$(references)

Arguments:
- `ρ0`: initial reduced density matrix
- `fbU`: time-series of forward-backward propagators
- `Jw`: array of spectral densities
- `svec`: diagonal elements of system operators through which the corresponding baths interact. QuAPI currently only works for baths with diagonal coupling to the system.
- `β`: inverse temperature
- `dt`: time-step for recording the density matrices
- `ntimes`: number of time steps of simulation
- `kmax`: number of steps within memory
- `extraargs`: extra arguments for the TEMPO algorithm. Contains the `cutoff` threshold for SVD filtration, the maximum bond dimension, `maxdim`, and the `algorithm` of applying an MPO to an MPS.
"""
function propagate(; fbU::AbstractArray{ComplexF64,3}, Jw::AbstractVector{T}, β::Real, ρ0::AbstractMatrix{ComplexF64}, dt::Real, ntimes::Int, kmax::Int, extraargs::TEMPOArgs=TEMPOArgs(), svec=[1.0 -1.0], reference_prop=false, verbose::Bool=false, exec=ThreadedEx()) where {T<:SpectralDensities.SpectralDensity}
    U0e = build_augmented_propagator(; fbU, Jw, β, dt, ntimes, kmax, extraargs, svec, reference_prop, verbose, exec)
    Utilities.apply_propagator(; propagators=U0e, ρ0, ntimes, dt)
end

end