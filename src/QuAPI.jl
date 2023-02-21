module QuAPI

using ..EtaCoefficients, ..SpectralDensities, ..Utilities

struct States
    sbar::Matrix{Float64}
    Δs::Matrix{Float64}
    forward_ind::Vector{UInt8}
    backward_ind::Vector{UInt8}
end

struct Path{T}
    states::Vector{UInt8}
    amplitude::T
    sdim2::UInt8
end
function next_path(p::Path{T}) where {T<:Number}
    states = p.states
    sdim2 = p.sdim2
    for j = 1:length(states)
        if states[j] < sdim2
            states[j] += 1
            break
        elseif states[j] == sdim2
            states[j] = 1
        end
    end
    return Path{T}(states, zero(T), sdim2), states != repeat([sdim2], length(states))
end

function hash_path(states::Vector{UInt8}, sdim)
    factor = 1
    number = 0
    for s in states
        number += (s - 1) * factor
        factor *= sdim
    end
    number + 1
end

@inbounds function setup_simulation(ρ0, η, svec, cutoff)
    sdim = size(ρ0, 1)
    sdim_square = sdim^2
    nbaths = size(svec, 1)
    sbar = zeros(nbaths, sdim_square)
    Δs = zeros(nbaths, sdim_square)
    forward_ind = zeros(UInt8, sdim_square)
    backward_ind = zeros(UInt8, sdim_square)
    amplitudes = zeros(ComplexF64, sdim_square)
    count = 1
    for forind = 1:size(svec, 2)
        for backind = 1:size(svec, 2)
            forward_ind[count] = forind
            backward_ind[count] = backind
            amplitudes[count] = ρ0[forind, backind]
            count += 1
        end
    end
    for (bath_num, bath_svec) in enumerate(eachrow(svec))
        count = 1
        for sf in bath_svec
            for sb in bath_svec
                sbar[bath_num, count] = (sf + sb) / 2.0
                Δs[bath_num, count] = sf - sb
                count += 1
            end
        end
    end
    state_values = States(sbar, Δs, forward_ind, backward_ind)
    paths = Vector{Path{ComplexF64}}()
    for count = 1:sdim_square
        states = [count]
        if abs(amplitudes[count]) > cutoff.cutoff
            amp = amplitudes[count]
            for (bn, bη) in enumerate(η)
                amp *= exp(-state_values.Δs[bn, count] * (real(bη.η00) * state_values.Δs[bn, count] + 2im * imag(bη.η00) * state_values.sbar[bn, count]))
            end
            push!(paths, Path{ComplexF64}(states, amp, sdim_square))
        end
    end
    state_values, paths
end

function get_influence(η::EtaCoefficients.EtaCoeffs, bath_number::Int, state_values::States, path::Vector{UInt8}, terminal::Bool, in_memory::Bool)
    @inbounds begin
        Δsfinal = state_values.Δs[bath_number, path[end]]
        if Δsfinal == 0
            return 1
        end
        num_time = length(path)
        if in_memory
            if terminal
                influence_exponent = real(η.η00) * Δsfinal + 2im * imag(η.η00) * state_values.sbar[bath_number, path[end]]
                for i = 2:num_time-1
                    influence_exponent += real(η.η0m[num_time-i]) * state_values.Δs[bath_number, path[i]] + 2im * imag(η.η0m[num_time-i]) * state_values.sbar[bath_number, path[i]]
                end
                influence_exponent += real(η.η0e[num_time-1]) * state_values.Δs[bath_number, path[1]] + 2im * imag(η.η0e[num_time-1]) * state_values.sbar[bath_number, path[1]]
                return exp(-Δsfinal * influence_exponent)
            else
                influence_exponent = real(η.ηmm) * Δsfinal + 2im * imag(η.ηmm) * state_values.sbar[bath_number, path[end]]
                for i = 2:num_time-1
                    influence_exponent += real(η.ηmn[num_time-i]) * state_values.Δs[bath_number, path[i]] + 2im * imag(η.ηmn[num_time-i]) * state_values.sbar[bath_number, path[i]]
                end
                influence_exponent += real(η.η0m[num_time-1]) * state_values.Δs[bath_number, path[1]] + 2im * imag(η.η0m[num_time-1]) * state_values.sbar[bath_number, path[1]]
                return exp(-Δsfinal * influence_exponent)
            end
        else
            if terminal
                influence_exponent = real(η.η00) * Δsfinal + 2im * imag(η.η00) * state_values.sbar[bath_number, path[end]]
                for i = 1:num_time-1
                    influence_exponent += real(η.η0m[num_time-i]) * state_values.Δs[bath_number, path[i]] + 2im * imag(η.η0m[num_time-i]) * state_values.sbar[bath_number, path[i]]
                end
                return exp(-Δsfinal * influence_exponent)
            else
                influence_exponent = real(η.ηmm) * Δsfinal + 2im * imag(η.ηmm) * state_values.sbar[bath_number, path[end]]
                for i = 1:num_time-1
                    influence_exponent += real(η.ηmn[num_time-i]) * state_values.Δs[bath_number, path[i]] + 2im * imag(η.ηmn[num_time-i]) * state_values.sbar[bath_number, path[i]]
                end
                return exp(-Δsfinal * influence_exponent)
            end
        end
    end
end

get_influence(η::Vector{EtaCoefficients.EtaCoeffs}, state_values::States, path::Vector{UInt8}, terminal::Bool, in_memory::Bool) = prod([get_influence(bη, bn, state_values, path, terminal, in_memory) for (bn, bη) in enumerate(η)])

"""
Filtration parameters for QuAPI. Currently has a threshold for magnitude-based filtering, with a default value of `cutoff=0` (no filtering).
"""
struct QuAPIArgs <: Utilities.ExtraArgs
    cutoff::Float64
end
QuAPIArgs(; cutoff=0.0) = QuAPIArgs(cutoff)

"""
    propagate(; fbU::Matrix{ComplexF64}, Jw::Vector{T}, β::Real, ρ0, dt::Real, ntimes::Int, kmax::Int, extraargs::QuAPIArgs=QuAPIArgs(), svec=[1.0 -1.0], reference_prop=false, verbose::Bool=false) where {T<:SpectralDensities.SpectralDensity}

Given a time-series of system forward-backward propagators, `fbU`, the spectral densities describing the solvent, `Jw`, and an inverse temperature, this uses QuAPI to propagate the input initial reduced density matrix, ρ0, with a time-step of `dt` for `ntimes` time steps. A non-Markovian memory of `kmax` steps is used in this simulation. The i^th bath, described by `Jw[i]`, interacts with the system through the diagonal operator with the values of `svec[j,:]`.

`ρ0`: initial reduced density matrix
`fbU`: time-series of forward-backward propagators
`Jw`: array of spectral densities
`svec`: diagonal elements of system operators through which the corresponding baths interact. QuAPI currently only works for baths with diagonal coupling to the system

`dt`: time-step for recording the density matrices
`ntimes`: number of time steps of simulation
`kmax`: number of steps within memory
`extraargs`: extra arguments for the QuAPI algorithm. Contains the filtration cutoff threshold
"""
function propagate(; fbU::AbstractArray{ComplexF64,3}, Jw::AbstractVector{T}, β::Real, ρ0::AbstractMatrix{ComplexF64}, dt::Real, ntimes::Int, kmax::Int, extraargs::QuAPIArgs=QuAPIArgs(), svec=[1.0 -1.0], reference_prop=false, verbose::Bool=false) where {T<:SpectralDensities.SpectralDensity}
    @assert length(Jw) == size(svec, 1)
    if kmax > ntimes
        kmax = ntimes + 2
    end
    η = [EtaCoefficients.calculate_η(jw; β, dt, kmax, imaginary_only=reference_prop) for jw in Jw]
    sdim = size(ρ0, 1)
    ρs = zeros(ComplexF64, ntimes + 1, sdim, sdim)
    ρs[1, :, :] = ρ0
    state_values, paths = setup_simulation(ρ0, η, svec, extraargs)
    eacp = false
    if kmax == 0
        kmax += 1
        eacp = true
    end

    sdim2 = sdim^2
    if verbose
        @info "Starting propagation within memory"
    end
    tmprho = zeros(ComplexF64, sdim^2)
    for i = 1:min(kmax, ntimes)
        if verbose
            @info "Step = $(i), #paths = $(length(paths))"
        end
        new_paths = Vector{Path{ComplexF64}}()
        for path in paths
            states = path.states
            tmprho .= zero(ComplexF64)
            @inbounds tmprho[states[end]] = 1.0
            tmprho = fbU[i, :, :] * tmprho * path.amplitude
            for (s, amp) in enumerate(tmprho)
                tmpstates = deepcopy(states)
                push!(tmpstates, s)
                amplitude = amp * get_influence(η, state_values, tmpstates, false, true)
                if eacp
                    @inbounds ρs[i+1, state_values.forward_ind[s], state_values.backward_ind[s]] += amp
                    amplitude = amp
                else
                    @inbounds ρs[i+1, state_values.forward_ind[s], state_values.backward_ind[s]] += amp * get_influence(η, state_values, tmpstates, true, true)
                end
                if abs(amplitude) > extraargs.cutoff
                    push!(new_paths, Path{ComplexF64}(tmpstates, amplitude, sdim2))
                end
            end
        end
        paths = new_paths
        new_paths = nothing
    end

    if verbose
        @info "Starting iteration"
    end
    path_max = repeat(Vector{UInt8}([sdim2]), length(paths[1].states) - 1)
    max_length = hash_path(path_max, sdim2)
    checked = repeat([false], max_length)
    tmat = zeros(ComplexF64, max_length)
    for i = kmax+1:ntimes
        if verbose
            @info "Step = $(i), #paths = $(length(paths))"
        end
        tmat .= zero(ComplexF64)
        for path in paths
            trunc_path = path.states[2:end]
            @inbounds tmat[hash_path(trunc_path, sdim2)] += path.amplitude
        end
        new_paths = Vector{Path{ComplexF64}}()
        checked .= false
        for path in paths
            @inbounds states = path.states[2:end]
            hash_val = hash_path(states, sdim2)
            if checked[hash_val]
                continue
            else
                tmprho .= zero(ComplexF64)
                @inbounds tmprho[states[end]] = 1.0
                tmprho = fbU[i, :, :] * tmprho * tmat[hash_val]
                for (s, amp) in enumerate(tmprho)
                    tmpstates = deepcopy(states)
                    push!(tmpstates, s)
                    amplitude = amp * get_influence(η, state_values, tmpstates, false, true)
                    if eacp
                        @inbounds ρs[i+1, state_values.forward_ind[s], state_values.backward_ind[s]] += amp
                        amplitude = amp
                    else
                        @inbounds ρs[i+1, state_values.forward_ind[s], state_values.backward_ind[s]] += amp * get_influence(η, state_values, tmpstates, true, true)
                    end
                    if abs(amplitude) > extraargs.cutoff
                        push!(new_paths, Path{ComplexF64}(tmpstates, amplitude, sdim2))
                    end
                end
                @inbounds checked[hash_val] = true
            end
        end
        paths = new_paths
        new_paths = nothing
    end
    0:dt:ntimes*dt, ρs
end

function get_path_influence(η::EtaCoefficients.EtaCoeffs, bath_number::Int, state_values, path)
    @inbounds begin
        interm_influence = zero(ComplexF64)
        i = length(path) - 1
        for sk = 2:i
            val = real(η.ηmm) * state_values.Δs[bath_number, path[sk]] + 2im * imag(η.ηmm) * state_values.sbar[bath_number, path[sk]]
            for skp = 2:sk-1
                val += real(η.ηmn[sk-skp]) * state_values.Δs[bath_number, path[skp]] + 2im * imag(η.ηmn[sk-skp]) * state_values.sbar[bath_number, path[skp]]
            end
            interm_influence += -state_values.Δs[bath_number, path[sk]] * val
        end
        amplitude = exp(interm_influence)

        init_influence_0 = -state_values.Δs[bath_number, path[1]] * (real(η.η00) * state_values.Δs[bath_number, path[1]] + 2im * imag(η.η00) * state_values.sbar[bath_number, path[1]])
        init_influence_m = -state_values.Δs[bath_number, path[1]] * (real(η.ηmm) * state_values.Δs[bath_number, path[1]] + 2im * imag(η.ηmm) * state_values.sbar[bath_number, path[1]])
        for j = 2:i
            init_influence_0 += -state_values.Δs[bath_number, path[j]] * (real(η.η0m[j-1] * state_values.Δs[bath_number, path[1]]) + 2im * imag(η.η0m[j-1]) * state_values.sbar[bath_number, path[1]])
            init_influence_m += -state_values.Δs[bath_number, path[j]] * (real(η.ηmn[j-1] * state_values.Δs[bath_number, path[1]]) + 2im * imag(η.ηmn[j-1]) * state_values.sbar[bath_number, path[1]])
        end
        init_influence_0 = exp(init_influence_0)
        init_influence_m = exp(init_influence_m)

        final_influence_0 = -state_values.Δs[bath_number, path[end]] * (real(η.η00) * state_values.Δs[bath_number, path[end]] + 2im * imag(η.η00) * state_values.sbar[bath_number, path[end]])
        final_influence_m = -state_values.Δs[bath_number, path[end]] * (real(η.ηmm) * state_values.Δs[bath_number, path[end]] + 2im * imag(η.ηmm) * state_values.sbar[bath_number, path[end]])
        for j = 2:i
            final_influence_0 += -state_values.Δs[bath_number, path[end]] * (real(η.η0m[i+1-j] * state_values.Δs[bath_number, path[j]]) + 2im * imag(η.η0m[i+1-j]) * state_values.sbar[bath_number, path[j]])
            final_influence_m += -state_values.Δs[bath_number, path[end]] * (real(η.ηmn[i+1-j] * state_values.Δs[bath_number, path[j]]) + 2im * imag(η.ηmn[i+1-j]) * state_values.sbar[bath_number, path[j]])
        end
        final_influence_0 = exp(final_influence_0)
        final_influence_m = exp(final_influence_m)

        term_influence_0e = exp(-state_values.Δs[bath_number, path[end]] * (real(η.η0e[i]) * state_values.Δs[bath_number, path[1]] + 2im * imag(η.η0e[i]) * state_values.sbar[bath_number, path[1]]))
        term_influence_me = exp(-state_values.Δs[bath_number, path[end]] * (real(η.η0m[i]) * state_values.Δs[bath_number, path[1]] + 2im * imag(η.η0m[i]) * state_values.sbar[bath_number, path[1]]))
        term_influence_0m = exp(-state_values.Δs[bath_number, path[end]] * (real(η.η0m[i]) * state_values.Δs[bath_number, path[1]] + 2im * imag(η.η0m[i]) * state_values.sbar[bath_number, path[1]]))
        term_influence_mn = exp(-state_values.Δs[bath_number, path[end]] * (real(η.ηmn[i]) * state_values.Δs[bath_number, path[1]] + 2im * imag(η.ηmn[i]) * state_values.sbar[bath_number, path[1]]))

        amplitude * init_influence_0 * final_influence_0 * term_influence_0e, amplitude * init_influence_0 * term_influence_0m, amplitude * init_influence_m * final_influence_0 * term_influence_me, amplitude * init_influence_m * term_influence_mn
    end
end

"""
    build_augmented_propagator(; fbU::Matrix{ComplexF64}, Jw::Vector{T}, β::Real, dt::Real, ntimes::Int, extraargs=QuAPIArgs(), svec=[1.0 -1.0], reference_prop=false, verbose::Bool=false) where {T<:SpectralDensities.SpectralDensity}
Builds the propagators, augmented with the influence of the harmonic baths defined by the spectral densities `Jw`,  upto `ntimes` time-steps without iteration. The paths are generated in full forward-backward space but not stored. So, while the space requirement is minimal and constant, the time complexity for each time-step grows by an additional factor of ``d^2``, where ``d`` is the dimensionality of the system. This i^th bath, described by `Jw[i]`, interacts with the system through the diagonal operator with the values of `svec[j,:]`.
"""
function build_augmented_propagator(; fbU::AbstractArray{ComplexF64,3}, Jw::Vector{T}, β::Real, dt::Real, ntimes::Int, extraargs::QuAPIArgs=QuAPIArgs(), svec=[1.0 -1.0], reference_prop=false, verbose::Bool=false, end_prop=false) where {T<:SpectralDensities.SpectralDensity}
    @assert length(Jw) == size(svec, 1)
    η = [EtaCoefficients.calculate_η(jw; β, dt, kmax=ntimes, imaginary_only=reference_prop) for jw in Jw]
    sdim2 = size(fbU, 2)
    sdim = trunc(Int, sqrt(sdim2))
    state_values, _ = setup_simulation(ones(sdim, sdim), η, svec, extraargs)

    if verbose
        @info "Starting propagation within memory"
    end
    U0e = zeros(ComplexF64, ntimes, sdim2, sdim2)
    U0m = zeros(ComplexF64, ntimes, sdim2, sdim2)
    Ume = zeros(ComplexF64, ntimes, sdim2, sdim2)
    Umn = zeros(ComplexF64, ntimes, sdim2, sdim2)
    for i = 1:ntimes
        if verbose
            @info "Step = $(i)"
        end
        num_paths = 0
        for path_num = 1:sdim2^(i+1)
            states = Utilities.unhash_path(path_num, i, sdim2)
            @inbounds state_pairs = collect(zip(states, states[2:end]))
            @inbounds bare_amplitude = prod([fbU[i, s[1], s[2]] for s in state_pairs])
            if abs(bare_amplitude) < extraargs.cutoff
                continue
            end
            num_paths += 1
            amplitudes = [bare_amplitude, bare_amplitude, bare_amplitude, bare_amplitude]
            for (bn, bη) in enumerate(η)
                influence = get_path_influence(bη, bn, state_values, states)
                amplitudes .*= influence
            end
            @inbounds U0e[i, states[end], states[1]] += amplitudes[1]
            @inbounds U0m[i, states[end], states[1]] += amplitudes[2]
            @inbounds Ume[i, states[end], states[1]] += amplitudes[3]
            @inbounds Umn[i, states[end], states[1]] += amplitudes[4]
        end
        if verbose
            @info "Done time step $(i). # paths = $(num_paths)."
        end
    end
    end_prop ? U0e : U0e, U0m, Ume, Umn
end

end