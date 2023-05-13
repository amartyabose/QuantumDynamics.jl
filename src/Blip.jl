module Blip

using ..EtaCoefficients, ..Propagators, ..SpectralDensities, ..Utilities

function setup_simulation(svec)
    nbaths = size(svec, 1)
    nstates = size(svec, 2)
    Δs = zeros(nbaths, nstates^2)
    sbar = zeros(nbaths, nstates^2)
    for (bn, sv) in enumerate(eachrow(svec))
        count = 1
        for sf in sv
            for sb in sv
                Δs[bn, count] = sf - sb
                sbar[bn, count] = (sf + sb) / 2.0
                count += 1
            end
        end
    end
    checked = zeros(Bool, nstates^2)
    group_states = Vector{Vector{UInt8}}()
    group_Δs = Vector{Vector{Float64}}()
    count = 1
    for (i, dels) in enumerate(eachcol(Δs))
        if checked[i]
            continue
        end
        push!(group_Δs, dels)
        states = [i]
        checked[i] = true
        for j = i+1:size(Δs, 2)
            if Δs[:, j] == dels
                push!(states, j)
                checked[j] = true
            end
        end
        push!(group_states, states)
    end
    group_Δs_final = zeros(nbaths, length(group_Δs))
    for s = 1:length(group_Δs)
        for b = 1:nbaths
            group_Δs_final[b, s] = group_Δs[s][b]
        end
    end
    state_to_blip_map = Vector{UInt8}()
    for j = 1:nstates^2
        for (s, gs) in enumerate(group_states)
            if !isnothing(findfirst(x -> x == j, gs))
                push!(state_to_blip_map, s)
                break
            end
        end
    end
    group_states, state_to_blip_map, group_Δs_final, sbar, Δs
end

function get_total_amplitude(; tmpprops, path, group_Δs, sbar, η, propagator_type, nsteps, sdim2, val1, valend, valjkp)
    @inbounds begin
        # val1 = zeros(ComplexF64, sdim2)
        # valend = zeros(ComplexF64, sdim2)
        # valjkp = zeros(ComplexF64, nsteps, sdim2)
        for (bn, bη) in enumerate(η)
            ηee = propagator_type == "0e" || propagator_type == "me" ? bη.η00 : zero(ComplexF64)
            η00 = propagator_type == "0e" || propagator_type == "0m" ? bη.η00 : bη.ηmm
            η0m = propagator_type == "0e" || propagator_type == "0m" ? bη.η0m : bη.ηmn
            ηme = propagator_type == "0e" || propagator_type == "me" ? bη.η0m : bη.ηmn
            η0e = propagator_type == "0e" ? bη.η0e : (propagator_type == "mn" ? bη.ηmn : bη.η0m)
            for j = 1:sdim2
                val1[j] += -group_Δs[bn, path[1]] * (real(ηee) * group_Δs[bn, path[1]] + 2im * imag(ηee) * sbar[bn, j])
                exponent = -group_Δs[bn, path[end]] * (real(η00) * group_Δs[bn, path[end]] + 2im * imag(η00) * sbar[bn, j]) - group_Δs[bn, path[1]] * (real(η0e[nsteps]) * group_Δs[bn, path[end]] + 2im * imag(η0e[nsteps]) * sbar[bn, j])
                for k = 2:nsteps
                    exponent += -group_Δs[bn, path[k]] * (real(η0m[nsteps+1-k]) * group_Δs[bn, path[end]] + 2im * imag(η0m[nsteps+1-k]) * sbar[bn, j])
                end
                valend[j] += exponent
            end
            for kp = 2:nsteps
                for j = 1:sdim2
                    exponent = -group_Δs[bn, path[1]] * (real(ηme[kp-1]) * group_Δs[bn, path[kp]] + 2im * imag(ηme[kp-1]) * sbar[bn, j]) - group_Δs[bn, path[kp]] * (real(bη.ηmm) * group_Δs[bn, path[kp]] + 2im * imag(bη.ηmm) * sbar[bn, j])
                    for k = 2:kp-1
                        exponent += -group_Δs[bn, path[k]] * (real(bη.ηmn[kp-k]) * group_Δs[bn, path[kp]] + 2im * imag(bη.ηmn[kp-k]) * sbar[bn, j])
                    end
                    valjkp[kp, j] += exponent
                end
            end
        end
        for j = 1:sdim2
            val1[j] = exp(val1[j])
            valend[j] = exp(valend[j])
            for k = 1:sdim2
                tmpprops[1, j, k] *= val1[j]
                tmpprops[nsteps, k, j] *= valend[j]
            end
        end
        valjkp .= exp.(valjkp)
        for k = 1:sdim2
            for j = 1:sdim2
                for kp = 2:nsteps
                    tmpprops[kp, j, k] *= valjkp[kp, j]
                end
            end
        end
        tmpprop = tmpprops[1, :, :]
        for j = 2:nsteps
            tmpprop *= tmpprops[j, :, :]
        end
    end
    tmpprop
end

"""
Filtration parameters for blips. Currently has the maximum number of blips allowed which by default is -1 (implying all blips are allowed).
"""
struct BlipArgs <: Utilities.ExtraArgs
    max_blips::Int
end
BlipArgs(; max_blips::Int=-1) = BlipArgs(max_blips)

"""
    build_augmented_propagator(; fbU::Matrix{ComplexF64}, Jw::Vector{T}, β::Real, dt::Real, ntimes::Int, extraargs::BlipArgs=BlipArgs(), svec=[1.0 -1.0], reference_prop=false, verbose::Bool=false) where {T<:SpectralDensities.SpectralDensity}
Builds the propagators, augmented with the influence of the harmonic baths defined by the spectral densities `Jw`,  upto `ntimes` time-steps without iteration using the **blip decomposition**. The paths are, consequently, generated in the space of unique blips and not stored. So, while the space requirement is minimal and constant, the time complexity for each time-step grows by an additional factor of ``b``, where ``b`` is the number of unique blip-values. The i^th bath, described by `Jw[i]`, interacts with the system through the diagonal operator with the values of `svec[j,:]`.
"""
function build_augmented_propagator(; fbU::AbstractArray{ComplexF64,3}, Jw::Vector{T}, β::Real, dt::Real, ntimes::Int, extraargs::BlipArgs=BlipArgs(), svec=[1.0 -1.0], reference_prop=false, verbose::Bool=false) where {T<:SpectralDensities.SpectralDensity}
    @assert length(Jw) == size(svec, 1)
    cutoff = extraargs.max_blips == -1 ? ntimes + 1 : extraargs.max_blips
    η = [EtaCoefficients.calculate_η(jw; β, dt, kmax=ntimes, imaginary_only=reference_prop) for jw in Jw]
    sdim2 = size(fbU, 2)
    group_states, _, group_Δs, sbar, _ = setup_simulation(svec)

    ndim = length(group_states)
    U0e = zeros(ComplexF64, ntimes, sdim2, sdim2)
    propagators = zeros(ComplexF64, ntimes, sdim2, sdim2)
    val1 = zeros(ComplexF64, sdim2)
    valend = zeros(ComplexF64, sdim2)
    @inbounds begin
        for i = 1:ntimes
            valjkp = zeros(ComplexF64, i, sdim2)
            if verbose
                @info "Starting time step $(i)."
            end
            num_paths = 0
            for b = 0:cutoff
                path_list = Utilities.unhash_path_blips(i, ndim, b)
                for path in path_list
                    num_paths += 1
                    val1 .= 0.0 + 0.0im
                    valend .= 0.0 + 0.0im
                    valjkp .= 0.0 + 0.0im
                    for (j, (sf, si)) in enumerate(zip(path, path[2:end]))
                        propagators[j, :, :] .= 0.0 + 0.0im
                        propagators[j, group_states[sf], group_states[si]] .= fbU[i, group_states[sf], group_states[si]]
                    end
                    U0e[i, :, :] .+= get_total_amplitude(; tmpprops=propagators, path, group_Δs, sbar, η, propagator_type="0e", nsteps=i, sdim2, val1, valend, valjkp)
                end
            end
            if verbose
                @info "Done time step $(i). # paths = $(num_paths)."
            end
        end
    end
    U0e
end

"""
    build_augmented_propagator_QuAPI_TTM(; fbU::Matrix{ComplexF64}, Jw::Vector{T}, β::Real, dt::Real, ntimes::Int, extraargs::BlipArgs=BlipArgs(), svec=[1.0 -1.0], reference_prop=false, verbose::Bool=false) where {T<:SpectralDensities.SpectralDensity}
Builds the propagators, augmented with the influence of the harmonic baths defined by the spectral densities `Jw`,  upto `ntimes` time-steps without iteration using the **blip decomposition**. The paths are, consequently, generated in the space of unique blips and not stored. So, while the space requirement is minimal and constant, the time complexity for each time-step grows by an additional factor of ``b``, where ``b`` is the number of unique blip-values. The i^th bath, described by `Jw[i]`, interacts with the system through the diagonal operator with the values of `svec[j,:]`. In this version, multiple ``types'' of propagators are calculated. These are required to make the TTM scheme consistent with QuAPI splitting.
"""
function build_augmented_propagator_QuAPI_TTM(; fbU::AbstractArray{ComplexF64,3}, Jw::Vector{T}, β::Real, dt::Real, ntimes::Int, extraargs::BlipArgs=BlipArgs(), svec=[1.0 -1.0], reference_prop=false, verbose::Bool=false) where {T<:SpectralDensities.SpectralDensity}
    @assert length(Jw) == size(svec, 1)
    cutoff = extraargs.max_blips == -1 ? ntimes + 1 : extraargs.max_blips
    η = [EtaCoefficients.calculate_η(jw; β, dt, kmax=ntimes, imaginary_only=reference_prop) for jw in Jw]
    sdim2 = size(fbU, 2)
    group_states, _, group_Δs, sbar, _ = setup_simulation(svec)

    ndim = length(group_states)
    U0e = zeros(ComplexF64, ntimes, sdim2, sdim2)
    U0m = zeros(ComplexF64, ntimes, sdim2, sdim2)
    Ume = zeros(ComplexF64, ntimes, sdim2, sdim2)
    Umn = zeros(ComplexF64, ntimes, sdim2, sdim2)
    propagators = zeros(ComplexF64, ntimes, sdim2, sdim2)
    tmpprops = zeros(ComplexF64, ntimes, sdim2, sdim2)
    val1 = zeros(ComplexF64, sdim2)
    valend = zeros(ComplexF64, sdim2)
    @show cutoff
    @inbounds begin
        for i = 1:ntimes
            valjkp = zeros(ComplexF64, i, sdim2)
            if verbose
                @info "Starting time step $(i)."
            end
            num_paths = 0
            for b = 0:cutoff
                path_list = Utilities.unhash_path_blips(i, ndim, b)
                for path in path_list
                    num_paths += 1
                    val1 .= 0.0 + 0.0im
                    valend .= 0.0 + 0.0im
                    valjkp .= 0.0 + 0.0im
                    for (j, (sf, si)) in enumerate(zip(path, path[2:end]))
                        propagators[j, :, :] .= 0.0 + 0.0im
                        propagators[j, group_states[sf], group_states[si]] .= fbU[i, group_states[sf], group_states[si]]
                    end
                    tmpprops .= propagators
                    U0e[i, :, :] .+= get_total_amplitude(; tmpprops, path, group_Δs, sbar, η, propagator_type="0e", nsteps=i, sdim2, val1, valend, valjkp)
                    tmpprops .= propagators
                    U0m[i, :, :] .+= get_total_amplitude(; tmpprops, path, group_Δs, sbar, η, propagator_type="0m", nsteps=i, sdim2, val1, valend, valjkp)
                    tmpprops .= propagators
                    Ume[i, :, :] .+= get_total_amplitude(; tmpprops, path, group_Δs, sbar, η, propagator_type="me", nsteps=i, sdim2, val1, valend, valjkp)
                    tmpprops .= propagators
                    Umn[i, :, :] .+= get_total_amplitude(; tmpprops, path, group_Δs, sbar, η, propagator_type="mn", nsteps=i, sdim2, val1, valend, valjkp)
                end
            end
            if verbose
                @info "Done time step $(i). # paths = $(num_paths)."
            end
        end
    end
    U0e, U0m, Ume, Umn
end

end
