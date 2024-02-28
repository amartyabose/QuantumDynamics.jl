module QuAPI

using HDF5
using FLoops
using ..EtaCoefficients, ..SpectralDensities, ..Utilities

const references = """
- Makri, N.; Makarov, D. E. Tensor Propagator for Iterative Quantum Time Evolution of Reduced Density Matrices. II. Numerical Methodology. The Journal of Chemical Physics 1995, 102 (11), 4611–4618. https://doi.org/10.1063/1.469508.
- Makri, N.; Makarov, D. E. Tensor Propagator for Iterative Quantum Time Evolution of Reduced Density Matrices. I. Theory. The Journal of Chemical Physics 1995, 102 (11), 4600–4610. https://doi.org/10.1063/1.469509."""

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

Given a time-series of system forward-backward propagators, `fbU`, the spectral densities describing the solvent, `Jw`, and an inverse temperature, this uses QuAPI to propagate the input initial reduced density matrix, ρ0, with a time-step of `dt` for `ntimes` time steps. A non-Markovian memory of `kmax` steps is used in this simulation. The j^th bath, described by `Jw[j]`, interacts with the system through the diagonal operator with the values of `svec[j,:]`.

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
    max_length = Utilities.hash_path(path_max, sdim2)
    checked = repeat([false], max_length)
    tmat = zeros(ComplexF64, max_length)
    for i = kmax+1:ntimes
        if verbose
            @info "Step = $(i), #paths = $(length(paths))"
        end
        tmat .= zero(ComplexF64)
        for path in paths
            trunc_path = path.states[2:end]
            @inbounds tmat[Utilities.hash_path(trunc_path, sdim2)] += path.amplitude
        end
        new_paths = Vector{Path{ComplexF64}}()
        checked .= false
        for path in paths
            @inbounds states = path.states[2:end]
            hash_val = Utilities.hash_path(states, sdim2)
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

function get_path_influence(η::EtaCoefficients.EtaCoeffs, bath_number::Int, state_values, path, quapi=true)
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

        quapi ? (amplitude * init_influence_0 * final_influence_0 * term_influence_0e, amplitude * init_influence_0 * term_influence_0m, amplitude * init_influence_m * final_influence_0 * term_influence_me, amplitude * init_influence_m * term_influence_mn) : amplitude * init_influence_0 * final_influence_0 * term_influence_0e
    end
end

"""
    build_augmented_propagator(; fbU::Matrix{ComplexF64}, Jw::Vector{T}, β::Real, dt::Real, ntimes::Int, kmax::Union{Int,Nothing}=nothing, extraargs=QuAPIArgs(), svec=[1.0 -1.0], reference_prop=false, verbose::Bool=false) where {T<:SpectralDensities.SpectralDensity}
Builds the propagators, augmented with the influence of the harmonic baths defined by the spectral densities `Jw`,  upto `ntimes` time-steps without iteration. The paths are generated in full forward-backward space but not stored. So, while the space requirement is minimal and constant, the time complexity for each time-step grows by an additional factor of ``d^2``, where ``d`` is the dimensionality of the system. This j^th bath, described by `Jw[j]`, interacts with the system through the diagonal operator with the values of `svec[j,:]`.
"""
function build_augmented_propagator(; fbU::AbstractArray{ComplexF64,3}, Jw::Vector{T}, β::Real, dt::Real, ntimes::Int, kmax::Union{Int,Nothing}=nothing, extraargs::QuAPIArgs=QuAPIArgs(), svec=[1.0 -1.0], reference_prop=false, verbose::Bool=false, output::Union{Nothing,HDF5.Group}=nothing, from_TTM::Bool=false) where {T<:SpectralDensities.SpectralDensity}
    @assert length(Jw) == size(svec, 1)
    η = [EtaCoefficients.calculate_η(jw; β, dt, kmax=ntimes, imaginary_only=reference_prop) for jw in Jw]
    sdim2 = size(fbU, 2)
    sdim = trunc(Int, sqrt(sdim2))
    state_values, _ = setup_simulation(ones(sdim, sdim), η, svec, extraargs)

    if verbose
        @info "Starting propagation within memory"
    end
    U0e = zeros(ComplexF64, ntimes, sdim2, sdim2)
    if !isnothing(output)
        if !from_TTM
            Utilities.check_or_insert_value(output, "U0e", U0e)
        end
        Utilities.check_or_insert_value(output, "time_taken", zeros(Float64, ntimes))
        Utilities.check_or_insert_value(output, "num_paths", zeros(Int64, ntimes))
    end
    for i = 1:ntimes
        if verbose
            @info "Step = $(i)"
        end
        num_paths = 0
        states = ones(UInt8, i+1)
        _, time_taken, memory_allocated, gc_time, _ = @timed begin
            for path_num = 1:sdim2^(i+1)
                Utilities.unhash_path(path_num, states, sdim2)
                bare_amplitude = one(ComplexF64)
                @inbounds for (j, (s1, s2)) in enumerate(zip(states, states[2:end]))
                    @inbounds bare_amplitude *= fbU[j, s1, s2]
                end
                if abs(bare_amplitude) < extraargs.cutoff
                    continue
                end
                num_paths += 1
                for (bn, bη) in enumerate(η)
                    @inbounds bare_amplitude *= get_path_influence(bη, bn, state_values, states, false)
                end
                @inbounds U0e[i, states[end], states[1]] += bare_amplitude
            end
        end
        if !isnothing(output)
            output["U0e"][i, :, :] = U0e[i, :, :]
            output["time_taken"][i] = time_taken
            output["num_paths"][i] = num_paths
            flush(output)
        end
        if verbose
            @info "Done time step $(i); # paths = $(sum(num_paths)); time = $(round(time_taken; digits=3)) sec; memory allocated = $(round(memory_allocated / 1e6; digits=3)) GB; gc time = $(round(gc_time; digits=3)) sec"
        end
    end
    U0e
end

"""
    build_augmented_propagator_parallel(; fbU::Matrix{ComplexF64}, Jw::Vector{T}, β::Real, dt::Real, ntimes::Int, kmax::Union{Int,Nothing}=nothing, extraargs=QuAPIArgs(), svec=[1.0 -1.0], reference_prop=false, verbose::Bool=false) where {T<:SpectralDensities.SpectralDensity}
Builds the propagators, augmented with the influence of the harmonic baths defined by the spectral densities `Jw`,  upto `ntimes` time-steps without iteration using shared memory parallelism. The paths are generated in full forward-backward space but not stored. So, while the space requirement is minimal and constant, the time complexity for each time-step grows by an additional factor of ``d^2``, where ``d`` is the dimensionality of the system. This j^th bath, described by `Jw[j]`, interacts with the system through the diagonal operator with the values of `svec[j,:]`.
"""
function build_augmented_propagator_parallel(; fbU::AbstractArray{ComplexF64,3}, Jw::Vector{T}, β::Real, dt::Real, ntimes::Int, kmax::Union{Int,Nothing}=nothing, extraargs::QuAPIArgs=QuAPIArgs(), svec=[1.0 -1.0], reference_prop=false, verbose::Bool=false, output::Union{Nothing,HDF5.Group}=nothing, from_TTM::Bool=false) where {T<:SpectralDensities.SpectralDensity}
    @assert length(Jw) == size(svec, 1)
    η = [EtaCoefficients.calculate_η(jw; β, dt, kmax=ntimes, imaginary_only=reference_prop) for jw in Jw]
    sdim2 = size(fbU, 2)
    sdim = trunc(Int, sqrt(sdim2))
    state_values, _ = setup_simulation(ones(sdim, sdim), η, svec, extraargs)

    if verbose
        @info "Starting propagation within memory"
    end
    U0e = Array{ComplexF64}(undef, ntimes, sdim2, sdim2)
    if !isnothing(output)
        if !from_TTM
            Utilities.check_or_insert_value(output, "U0e", U0e)
        end
        Utilities.check_or_insert_value(output, "time_taken", zeros(Float64, ntimes))
        Utilities.check_or_insert_value(output, "num_paths", zeros(Int64, ntimes))
    end
    for i = 1:ntimes
        if verbose
            @info "Step = $(i)"
        end
        _, time_taken, memory_allocated, gc_time, _ = @timed begin
            @floop for path_num = 1:sdim2^(i+1)
                @init states = ones(UInt8, i+1)
                Utilities.unhash_path(path_num, states, sdim2)
                bare_amplitude = one(ComplexF64)
                @inbounds for (j, (s1, s2)) in enumerate(zip(states, states[2:end]))
                    @inbounds bare_amplitude *= fbU[j, s1, s2]
                end
                if abs(bare_amplitude) < extraargs.cutoff
                    continue
                end
                @reduce num_paths = 0 + 1
                for (bn, bη) in enumerate(η)
                    @inbounds bare_amplitude *= get_path_influence(bη, bn, state_values, states, false)
                end
                @init tmpval = zeros(ComplexF64, sdim2, sdim2)
                tmpval .= 0.0
                @inbounds tmpval[states[end], states[1]] = bare_amplitude
                @reduce tmpU0e = zeros(ComplexF64, sdim2, sdim2) .+ tmpval
            end
            @inbounds U0e[i, :, :] .= tmpU0e
        end
        if !isnothing(output)
            output["U0e"][i, :, :] = U0e[i, :, :]
            output["time_taken"][i] = time_taken
            output["num_paths"][i] = num_paths
            flush(output)
        end
        if verbose
            @info "Done time step $(i); # paths = $(sum(num_paths)); time = $(round(time_taken; digits=3)) sec; memory allocated = $(round(memory_allocated / 1e6; digits=3)) GB; gc time = $(round(gc_time; digits=3)) sec"
        end
    end
    U0e
end

"""
    build_augmented_propagator_QuAPI_TTM(; fbU::Matrix{ComplexF64}, Jw::Vector{T}, β::Real, dt::Real, ntimes::Int, kmax::Union{Int,Nothing}=nothing, extraargs=QuAPIArgs(), svec=[1.0 -1.0], reference_prop=false, verbose::Bool=false) where {T<:SpectralDensities.SpectralDensity}
Builds the propagators, augmented with the influence of the harmonic baths defined by the spectral densities `Jw`,  upto `ntimes` time-steps without iteration. The paths are generated in full forward-backward space but not stored. So, while the space requirement is minimal and constant, the time complexity for each time-step grows by an additional factor of ``d^2``, where ``d`` is the dimensionality of the system. This j^th bath, described by `Jw[j]`, interacts with the system through the diagonal operator with the values of `svec[j,:]`. In this version, multiple ``types'' of propagators are calculated. These are required to make the TTM scheme consistent with QuAPI splitting.
"""
function build_augmented_propagator_QuAPI_TTM(; fbU::AbstractArray{ComplexF64,3}, Jw::Vector{T}, β::Real, dt::Real, ntimes::Int, kmax::Union{Int,Nothing}=nothing, extraargs::QuAPIArgs=QuAPIArgs(), svec=[1.0 -1.0], reference_prop=false, verbose::Bool=false, output::Union{Nothing,HDF5.Group}=nothing, from_TTM::Bool=false) where {T<:SpectralDensities.SpectralDensity}
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
    if !isnothing(output)
        if !from_TTM
            Utilities.check_or_insert_value(output, "U0e", U0e)
        end
        Utilities.check_or_insert_value(output, "time_taken", zeros(Float64, ntimes))
        Utilities.check_or_insert_value(output, "num_paths", zeros(Int64, ntimes))
    end
    for i = 1:ntimes
        if verbose
            @info "Step = $(i)"
        end
        _, time_taken, memory_allocated, gc_time, _ = @timed begin
            num_paths = 0
            states = ones(UInt8, i+1)
            for path_num = 1:sdim2^(i+1)
                Utilities.unhash_path(path_num, states, sdim2)
                bare_amplitude = one(ComplexF64)
                for (j, (s1, s2)) in enumerate(zip(states, states[2:end]))
                    @inbounds bare_amplitude *= fbU[j, s1, s2]
                end
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
        end
        if !isnothing(output)
            output["U0e"][i, :, :] = U0e[i, :, :]
            output["time_taken"][i] = time_taken
            output["num_paths"][i] = num_paths
            flush(output)
        end
        if verbose
            @info "Done time step $(i); # paths = $(sum(num_paths)); time = $(round(time_taken; digits=3)) sec; memory allocated = $(round(memory_allocated / 1e6; digits=3)) GB; gc time = $(round(gc_time; digits=3)) sec"
        end
    end
    U0e, U0m, Ume, Umn
end

"""
    build_augmented_propagator_QuAPI_TTM_parallel(; fbU::Matrix{ComplexF64}, Jw::Vector{T}, β::Real, dt::Real, ntimes::Int, kmax::Union{Int,Nothing}=nothing, extraargs=QuAPIArgs(), svec=[1.0 -1.0], reference_prop=false, verbose::Bool=false) where {T<:SpectralDensities.SpectralDensity}
Builds the propagators, augmented with the influence of the harmonic baths defined by the spectral densities `Jw`,  upto `ntimes` time-steps without iteration using shared memory parallelism. The paths are generated in full forward-backward space but not stored. So, while the space requirement is minimal and constant, the time complexity for each time-step grows by an additional factor of ``d^2``, where ``d`` is the dimensionality of the system. This j^th bath, described by `Jw[j]`, interacts with the system through the diagonal operator with the values of `svec[j,:]`. In this version, multiple ``types'' of propagators are calculated. These are required to make the TTM scheme consistent with QuAPI splitting.
"""
function build_augmented_propagator_QuAPI_TTM_parallel(; fbU::AbstractArray{ComplexF64,3}, Jw::Vector{T}, β::Real, dt::Real, ntimes::Int, kmax::Union{Int,Nothing}=nothing, extraargs::QuAPIArgs=QuAPIArgs(), svec=[1.0 -1.0], reference_prop=false, verbose::Bool=false, output::Union{Nothing,HDF5.Group}=nothing, from_TTM::Bool=false) where {T<:SpectralDensities.SpectralDensity}
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
    if !isnothing(output)
        if !from_TTM
            Utilities.check_or_insert_value(output, "U0e", U0e)
        end
        Utilities.check_or_insert_value(output, "time_taken", zeros(Float64, ntimes))
        Utilities.check_or_insert_value(output, "num_paths", zeros(Int64, ntimes))
    end
    for i = 1:ntimes
        if verbose
            @info "Step = $(i)"
        end
        _, time_taken, memory_allocated, gc_time, _ = @timed begin
            @floop for path_num = 1:sdim2^(i+1)
                @init states = ones(UInt8, i+1)
                Utilities.unhash_path(path_num, states, sdim2)
                bare_amplitude = one(ComplexF64)
                for (j, (s1, s2)) in enumerate(zip(states, states[2:end]))
                    @inbounds bare_amplitude *= fbU[j, s1, s2]
                end
                if abs(bare_amplitude) < extraargs.cutoff
                    continue
                end
                @reduce num_paths = 0 + 1
                @init amplitudes = zeros(ComplexF64, 4)
                amplitudes .= bare_amplitude
                for (bn, bη) in enumerate(η)
                    amplitudes .*= get_path_influence(bη, bn, state_values, states)
                end

                @init tmpval = zeros(ComplexF64, sdim2, sdim2)
                tmpval .= 0.0
                tmpval[states[end], states[1]] = amplitudes[1]
                @reduce tmpU0e = zeros(ComplexF64, sdim2, sdim2) .+ tmpval
                tmpval .= 0.0
                tmpval[states[end], states[1]] = amplitudes[2]
                @reduce tmpU0m = zeros(ComplexF64, sdim2, sdim2) .+ tmpval
                tmpval .= 0.0
                tmpval[states[end], states[1]] = amplitudes[3]
                @reduce tmpUme = zeros(ComplexF64, sdim2, sdim2) .+ tmpval
                tmpval .= 0.0
                tmpval[states[end], states[1]] = amplitudes[4]
                @reduce tmpUmn = zeros(ComplexF64, sdim2, sdim2) .+ tmpval
            end
            @inbounds U0e[i, :, :] .= tmpU0e
            @inbounds U0m[i, :, :] .= tmpU0m
            @inbounds Ume[i, :, :] .= tmpUme
            @inbounds Umn[i, :, :] .= tmpUmn
        end
        if !isnothing(output)
            output["U0e"][i, :, :] = U0e[i, :, :]
            output["time_taken"][i] = time_taken
            output["num_paths"][i] = num_paths
            flush(output)
        end
        if verbose
            @info "Done time step $(i); # paths = $(sum(num_paths)); time = $(round(time_taken; digits=3)) sec; memory allocated = $(round(memory_allocated / 1e6; digits=3)) GB; gc time = $(round(gc_time; digits=3)) sec"
        end
    end
    U0e, U0m, Ume, Umn
end

end